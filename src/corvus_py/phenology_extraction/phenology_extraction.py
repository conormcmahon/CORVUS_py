import ee 



# Extract smoothed, evenly-spaced phenology values from image collection
#    input_collection        image collection where each image contains bands DOY (day-of-year) and BAND_NAME
#    advance_window_width    width of search window before target date (in days)
#    following_window_width  width of search window after target date (in days)
#    num_time_steps          integer number of steps at which to fit phenology (e.g. 12 for monthly, 52 for weekly)
#    band_name               string containing name of target band (e.g. NDVI, SR_1, T1, etc.)
#    min_date                specify minimum date for range of fitted phenology values (note - imagery which occurs before that date
#                              won't be sampled for output datestamps, but can participate in the sampling window around a target date
#    max_date                specify maximum date for range of fitted phenology values (note - see above)
#    min_value               minimum possible value for parameter of interest (e.g. -1.0 for NDVI, 0 for Kelvin, etc.). Defaults to -1.0
#    max_value               maximum possible value for parameter of interest (e.g. 1.0 for NDVI). Defaults to 1.0
def fit_phenology(input_collection, advance_window_width, following_window_width, num_time_steps, band_name, min_date, max_date, min_value=-1.0, max_value=1.0, wrap_data=False):
    input_collection = input_collection.select(['time', band_name])
    collection_median = input_collection.median()
    
    # List of target dates (based on input number of timestamps to model)
    timing_interval = (max_date.difference(min_date, 'day')).divide(num_time_steps).floor()
    def sampleDate(index):
        return min_date.advance(ee.Number(index).multiply(timing_interval), 'day')
    target_dates = ee.List.sequence(0, num_time_steps).map(sampleDate)
    
    # Generate and return phenology ImageCollection
    # Contains NUM_TIME_STEPS images 
    # Each image contains the following bands:
    #    median:                 median of neighborhood values
    #    time:                   timestamp of image (in milliseconds since Jan 1, 1970)
    #    prediction_filtered:    predicted values at DOY 
    def fitAllDays(target_date):
        return(get_window_fit(input_collection, advance_window_width, following_window_width, ee.Date(target_date), band_name, min_value, max_value, collection_median, wrap_data))
    
    return ee.ImageCollection(target_dates.map(fitAllDays))



# Retrieves prediction at a given day-of-year
def get_window_fit(input_collection, advance_window_width, following_window_width, target_date, band_name, min_value=-1.0, max_value=1.0, collection_median=None, wrap_data=False):
    # Generate collection median, if unprovided
    if collection_median is None: 
        collection_median = input_collection.median()
    
    target_time = (ee.Number(target_date.difference(ee.Date("1970-01-01"), 'second'))
                                          .divide(365.25*24*3600).add(1970))
    
    # ************************ Temporal Window of Reference Scenes ************************
    # Set up first and last day within window 
    window_start_date = target_date.advance(-advance_window_width, 'day');
    window_end_date = target_date.advance(following_window_width, 'day');
    
    # Normal advance/following windows just collect data from before and after the target date
    advance_window = input_collection.filterDate(window_start_date, target_date)
    following_window = input_collection.filterDate(target_date, window_end_date)
    
    # "Wrapped" windows also collect data from the end or start of the previous or following year, for dates close to the start/end of the year
    target_year = target_date.get('year')
    target_year_start = ee.Date.fromYMD(target_year,1,1)
    target_year_end = ee.Date.fromYMD(target_year,12,31)
    def addYear(img):
        return  img.select(band_name).addBands((img.select('time').add(1)).toFloat())
    def subtractYear(img):
        return  img.select(band_name).addBands((img.select('time').subtract(1)).toFloat())
    following_window_wrapped = (input_collection.filterDate(target_year_start,
                                                            target_date.advance(-1,'year').advance(following_window_width,'day'))
                                                .map(addYear))
    advance_window_wrapped = (input_collection.filterDate(target_date.advance(1,'year').advance(-advance_window_width,'day'), 
                                                           target_year_end)
                                                .map(subtractYear))
  
    # Wrapped data are only included if an option is specified at function input
    window = ee.Algorithms.If(wrap_data,
                     advance_window.merge(following_window).merge(advance_window_wrapped).merge(following_window_wrapped),
                     advance_window.merge(following_window))
    window = ee.ImageCollection(window)
    
    # If there are NO images in the entire window, then some logic below fails. 
    # To prevent this, augment the collection with a single image which is 0 and masked everywhere
    window = window.merge(ee.ImageCollection(collection_median.multiply(ee.Image(0.0))
                                                            .mask(ee.Image(0)))
                                                            .cast(collection_median.bandTypes(), ['time', band_name]))
  
  
    # ************************ Regression Fitting and Prediction ************************
    # Add squared day-of-year variable (for quadratic fit)
    def addTimeSquared(img):
        return img.addBands(img.select('time').multiply(img.select('time')).rename('time_sq'))
    # Add constant term (for all regression fits)
    def addConstant(img):
        return img.addBands(ee.Image(1))
    # Generate regression data - for each point, image with: response variable, constant 1, DOY, and DOY^2
    regression_data = window.map(addTimeSquared).map(addConstant)
    
    # Fit quadratic regression to scenes in window 
    #   "robustLinearRegression() uses a cost function based on regression residuals to iteratively 
    #      de-weight outliers in the data (O’Leary, 1990)."
    quadratic_regression = regression_data.select(['constant','time','time_sq',band_name]).reduce(
            ee.Reducer.robustLinearRegression(numX=3, numY=1)
    )
    # Fit linear regression to scenes in window
    linear_regression = regression_data.select(['constant','time',band_name]).reduce(
        ee.Reducer.robustLinearRegression(numX=2, numY=1)
    )
    
    # Band names for variables in quadratic and linear fits
    quadratic_bands = [['constant','time','time_sq'],[band_name]]
    linear_bands = [['constant','time'],[band_name]]
    
    # The 3x1 regression outputs a 3x1 matrix, which has to be 'flattened' into a 3-band image
    quadratic_coeffs = quadratic_regression.select('coefficients').arrayFlatten(quadratic_bands)
    linear_coeffs = linear_regression.select('coefficients').arrayFlatten(linear_bands)
    
    # Apply both regression models to predict the result at the target date
    quadratic_prediction = (quadratic_coeffs.select('constant_'+band_name)
                                     .add(quadratic_coeffs.select('time_'+band_name).multiply(ee.Image(target_time)))
                                     .add(quadratic_coeffs.select('time_sq_'+band_name).multiply(ee.Image(target_time)).multiply(ee.Image(target_time))))
    linear_prediction = (linear_coeffs.select('constant_'+band_name)
                                     .add(linear_coeffs.select('time_'+band_name).multiply(ee.Image(target_time))))

    # ************************ Choose Fit by Number of Scenes ************************
    # Get mask for cases when there is too little data, so regressions will be unreliable
    #   Only use quadratic regression if count(scenes) > 6
    #   Only use linear regression if count(scense) > 3
    # This prevents overfitting to noise in cases with just a few values
    too_little_data_for_linear_mask = regression_data.select(band_name).count().lt(3).unmask();
    too_little_data_for_quadratic_mask = regression_data.select(band_name).count().lt(6).unmask();
    
    # At locations where there is too little data (too few uncloudy scenes) for the model to work, 
    #   switch from quadratic to linear regression or the median
    # Applies masks evaluated above
    model_prediction_filtered = ((quadratic_prediction.multiply(too_little_data_for_quadratic_mask.eq(0))).unmask()
                                                    .add(too_little_data_for_quadratic_mask.multiply(linear_prediction)))
    model_prediction_filtered = ((model_prediction_filtered.multiply(too_little_data_for_linear_mask.eq(0))).unmask()
                                                         .add(too_little_data_for_linear_mask.multiply(window.select(band_name).median())))
    
    # ************************ Reject Implausible Values ************************
    # At locations where the fit values are implausible, use the medians instead 
    # NOTE NDVI is valued between -1 and 1, so any values outside that range are physically unreasonable
    #   This can be a problem if there are only a few points and they set up a best-fit line that's vertical
    #   Optionally the user can provide a min_value and max_value which constrain the physical limits on the variable being fit
    too_high_mask = model_prediction_filtered.gt(max_value).unmask()
    too_low_mask = model_prediction_filtered.lt(min_value).unmask()
    
    # ************************ Reject Statistically Unlikely Values ************************
    # Also check if values are statistically unreasonable when compared to other values in neighborhood
    neighborhood_median = window.select(band_name).median()
    neighborhood_min = window.select(band_name).min()
    neighborhood_max = window.select(band_name).max()
    neighborhood_stdev = window.select(band_name).reduce(ee.Reducer.stdDev())
    #    Here, using a Z-score threshold of 1.5 to remove predictions which are outliers relative to
    #       temporal neighborhood window. If ((prediction - mean) > 1.5 * stdev), remove it (see below)
    statistical_outlier_high = model_prediction_filtered.gt(neighborhood_median.add(neighborhood_stdev.multiply(ee.Number(1.5)))).unmask()
    statistical_outlier_low = model_prediction_filtered.lt(neighborhood_median.subtract(neighborhood_stdev.multiply(ee.Number(1.5)))).unmask()
    
    # If the point is too high or too low, replace it with:
    #    median(neighborhood) +- 1.5 stdev(neighborhood) 
    # Direction of offset based on whether prediction was high or low relative to median
    bad_data_high = too_high_mask.add(statistical_outlier_high).gt(0)
    bad_data_low = too_low_mask.add(statistical_outlier_low).gt(0)
    model_prediction_filtered = (model_prediction_filtered.multiply(bad_data_high.eq(0))).unmask().add(bad_data_high.multiply(neighborhood_median).add(neighborhood_stdev))
    model_prediction_filtered = (model_prediction_filtered.multiply(bad_data_low.eq(0))).unmask().add(bad_data_low.multiply(neighborhood_median).subtract(neighborhood_stdev))
    
    # ************************ Return result ************************
    # Output results with median, day-of-year, and predicted value
    def unmask_image(img):
        return img.unmask()
    return (ee.Image(neighborhood_median)
            .addBands(ee.Image(target_time).float())
            .addBands(model_prediction_filtered)
            .addBands(window.select(band_name).count())
            .addBands(window.select(band_name).map(unmask_image).count())
            .addBands(neighborhood_min)
            .addBands(neighborhood_max)
            .addBands(neighborhood_stdev)
            .addBands(linear_coeffs.select('constant_'+band_name))
            .addBands(quadratic_coeffs.select('constant_'+band_name))
            .rename(['median','time','prediction_filtered','clear_images','window_size','min','max','stdev','linear_coef','quadratic_coef']))

  
  
  
  



'''


# Extract smoothed, evenly-spaced phenology values from image collection
#    input_collection   image collection where each image contains bands DOY (day-of-year) and BAND_NAME
#    window_radius      half-width of time window in which to perform fitting (before/after target date)
#    num_timestamps     integer number of steps at which to fit phenology (e.g. 12 for monthly, 52 for weekly)
#    band_name          string containing name of target band (e.g. NDVI, SR_1, T1, etc.)
#    min_value          minimum possible value for parameter of interest (e.g. -1.0 for NDVI, 0 for Kelvin, etc.). Defaults to -1.0
#    max_value          maximum possible value for parameter of interest (e.g. 1.0 for NDVI). Defaults to 1.0
def fit_phenology(input_collection, window_radius, num_timestamps, band_name, min_value=-1.0, max_value=1.0):
    input_collection = input_collection.select(['DOY', band_name])
    collection_median = input_collection.median()

    # List of target dates (based on input number of timestamps to model)
    days_of_year = ee.List.sequence(start = 0,
                                    step = ee.Number(365.25).divide(num_timestamps).floor(),
                                    count = num_timestamps)

    # Generate and return phenology ImageCollection
    # Contains NUM_TIMESTEPS images
    # Each image contains the following bands:
    #    median:                 median of window values
    #    DOY:                    day-of-year
    #    prediction_filtered:    predicted values at DOY
    def fit_all_days(target_doy):
            return(get_window_fit(input_collection, window_radius, band_name, target_doy, collection_median, min_value, max_value))
    
    return( ee.ImageCollection(days_of_year.map(fit_all_days)))




# Retrieves prediction at a given day-of-year
def get_window_fit(input_collection, window_radius, band_name, target_doy, collection_median, min_value=-1.0, max_value=1.0):
    # ************************ Temporal Window of Reference Scenes ************************
    # Convert input, client-side value to Earth Engine object
    target_doy = ee.Number(target_doy)
    # Set up first and last day within window
    start_doy = target_doy.subtract(window_radius)
    end_doy = target_doy.add(window_radius)
    # Roll over values that go past the start/end of a year
    #    If the start date is < 0, convert to a late-season value (e.g. -10 is 355)
    start_doy = start_doy.add(start_doy.lt(0).multiply(365)) 
    #    If the end date is > 365, convert to an early-season value (e.g. 375 is 10)
    end_doy = end_doy.subtract(end_doy.gt(365).multiply(365))
    # Extract relevant scenes from collection
    #    Remember, DOY for start of window can be negative if period overlaps 0...
    window_early = input_collection.filter(ee.Filter.dayOfYear(0,end_doy))
    #    OR the end DOY can be > 365, on the other extreme...
    window_late = input_collection.filter(ee.Filter.dayOfYear(start_doy,366))
    #    OR it can be in the middle, with neither end clipped:
    window_mid = input_collection.filter(ee.Filter.dayOfYear(start_doy,end_doy))
    # Now block the three windows together
    #    Again we need a bit of extra logic to handle windows that overlap ends of the year:
    window = ee.Algorithms.If(target_doy.subtract(window_radius).lt(ee.Number(0)),
                              window_early.merge(window_late),
                              window_mid)
    # Finally, our image collection with all relevant scenes:
    window = ee.ImageCollection(window)
    # If there are NO images in the entire window, then some logic below fails.
    # To prevent this, buffer the collection with a single image which is 0 and masked everywhere
    window = window.merge(ee.ImageCollection(collection_median.multiply(ee.Image(0.0))
                                                              .mask(ee.Image(0)))
                            .cast(collection_median.bandTypes(), ['DOY', band_name]))
    window_median = window.select(band_name).median()

    # ************************ Regression Fitting and Prediction ************************
    # Add squared day-of-year variable (for quadratic fit)
    def add_time_squared(img):
            return img.addBands(img.select('DOY').multiply(img.select('DOY')).rename('DOY_sq'))
    # Add constant term (for all regression fits)
    def add_constant(img):
            return img.addBands(ee.Image(1))
    # Generate regression data - for each point, image with: response variable, constant 1, DOY, and DOY^2
    regression_data = (window.map(add_time_squared)
                             .map(add_constant))

    # Fit quadratic regression to scenes in window
    #   "robustLinearRegression() uses a cost function based on regression residuals to iteratively
    #      de-weight outliers in the data (O’Leary, 1990)."
    quadratic_regression = regression_data.select(['constant','DOY','DOY_sq',band_name]).reduce(
        ee.Reducer.robustLinearRegression(numX = 3, 
                                          numY = 1)
    )
    # Fit linear regression to scenes in window
    linear_regression = regression_data.select(['constant','DOY',band_name]).reduce(
        ee.Reducer.robustLinearRegression(numX = 2, 
                                          numY = 1)
    )

    # Band names for variables in quadratic and linear fits
    quadratic_bands = [['constant','DOY','DOY_sq'],[band_name]]
    linear_bands = [['constant','DOY'],[band_name]]

    # The 3x1 regression outputs a 3x1 matrix, which has to be 'flattened' into a 3-band image
    quadratic_coeffs = quadratic_regression.select('coefficients').arrayFlatten(quadratic_bands)
    linear_coeffs = linear_regression.select('coefficients').arrayFlatten(linear_bands)

    # Apply both regression models to predict the result at the target date
    target_doy_img = ee.Image(target_doy)
    target_doy_squared = target_doy_img.multiply(target_doy_img)
    quadratic_prediction = (quadratic_coeffs.select('constant_'+band_name)
                                            .add(quadratic_coeffs.select('DOY_'+band_name).multiply(target_doy_img))
                                            .add(quadratic_coeffs.select('DOY_sq_'+band_name).multiply(target_doy_squared)))
    linear_prediction = (linear_coeffs.select('constant_'+band_name)
                                      .add(linear_coeffs.select('DOY_'+band_name).multiply(target_doy_img)))

    # ************************ Choose Fit by Number of Scenes ************************
    # Get mask for cases when there is too little data, so regressions will be unreliable
    #   Only use quadratic regression if count(scenes) > 6
    #   Only use linear regression if count(scense) > 3
    # This prevents overfitting to noise in cases with just a few values
    too_little_data_for_linear_mask = regression_data.select(band_name).count().lt(3).unmask()
    too_little_data_for_quadratic_mask = regression_data.select(band_name).count().lt(6).unmask()

    # At locations where there is too little data (too few uncloudy scenes) for the model to work,
    #   switch from quadratic to linear regression or the median
    # Use quadratic prediction, but fall back to linear if there's too little data:
    model_prediction_filtered = ((quadratic_prediction.multiply(too_little_data_for_quadratic_mask.eq(0))).unmask()
                                                      .add(too_little_data_for_quadratic_mask.multiply(linear_prediction)))
    # And if there's too little even for that, use the window median 
    model_prediction_filtered = ((model_prediction_filtered.multiply(too_little_data_for_linear_mask.eq(0))).unmask()
                                                           .add(too_little_data_for_linear_mask.multiply(window_median)))

    # ************************ Reject Implausible Values ************************
    # At locations where the fit values are implausible, use the medians instead
    # NOTE NDVI is valued between -1 and 1, so any values outside that range are physically unreasonable
    #   This can be a problem if there are only a few points and they set up a best-fit line that's vertical
    #   Optionally the user can provide a min_value and max_value which constrain the physical limits on the variable being fit
    too_high_mask = model_prediction_filtered.gt(max_value).unmask()
    too_low_mask = model_prediction_filtered.lt(min_value).unmask()

    # ************************ Reject Statistically Unlikely Values ************************
    # Also check if values are statistically unreasonable when compared to other values in window
    #window_min = window.select(band_name).min()
    #window_max = window.select(band_name).max()
    window_stdev = window.select(band_name).reduce(ee.Reducer.stdDev())
    #    Here, using a Z-score threshold of 1.5 to remove predictions which are outliers relative to
    #       temporal window. If ((prediction - mean) > 1.5 * stdev), remove it (see below)
    statistical_outlier_high = model_prediction_filtered.gt(window_median.add(window_stdev.multiply(ee.Number(1.5)))).unmask()
    statistical_outlier_low = model_prediction_filtered.lt(window_median.subtract(window_stdev.multiply(ee.Number(1.5)))).unmask()

    # If the point is too high or too low, replace it with:
    #    median(window) +- 1.5 stdev(window)
    # Direction of offset based on whether prediction was high or low relative to median
    bad_data_high = too_high_mask.add(statistical_outlier_high).gt(0)
    bad_data_low = too_low_mask.add(statistical_outlier_low).gt(0)
    model_prediction_filtered = ((model_prediction_filtered.multiply(bad_data_high.eq(0))).unmask()
                                                           .add(bad_data_high.multiply(window_median).add(window_stdev)))
    model_prediction_filtered = ((model_prediction_filtered.multiply(bad_data_low.eq(0))).unmask()
                                                           .add(bad_data_low.multiply(window_median).subtract(window_stdev)))

    # ************************ Return result ************************
    # Output results with median, day-of-year, and predicted value
    return (ee.Image(window_median)
              .addBands(ee.Image(target_doy)
              .addBands(model_prediction_filtered))
              .rename(['median','DOY','prediction_filtered']))


'''