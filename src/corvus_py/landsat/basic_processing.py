
import ee 

# *******************************************************************************************************
# ******************************************** CLOUD MASKING ********************************************
# *******************************************************************************************************
# Landsat 4,5,7 Cloud Mask - see QA_PIXEL bit designations at:
#  Landsat 4: https://developers.google.com/earth-engine/datasets/catalog/LANDSAT_LT04_C02_T1_L2
#  Landsat 5: https://developers.google.com/earth-engine/datasets/catalog/LANDSAT_LT05_C02_T1_L2
#  Landsat 7: https://developers.google.com/earth-engine/datasets/catalog/LANDSAT_LE07_C02_T1_L2
# NOTE - the difference between the 4/5/7 and 8/9 series is addition of Cirrus information in 8/9
# if mask_mask is set to ee.Image(0), then the mask will not change
def landsat_cloud_mask_457(image, mask_mask=1):
    # Bit-shift to get bits which need to be zero-valued to avoid masking
    dilated_cloud = (1 << 1)         
    cloud = (1 << 3)                 
    cloud_shadow = (1 << 4)          
    qa = image.select('QA_PIXEL');  
    mask = (qa.bitwiseAnd(dilated_cloud)
              .Or(qa.bitwiseAnd(cloud_shadow))
              .Or(qa.bitwiseAnd(cloud))
              .multiply(mask_mask))
    return image.updateMask(mask.Not())

# Landsat 8,9 Cloud Mask - see QA_PIXEL bit designations at:
#  Landsat 8: https://developers.google.com/earth-engine/datasets/catalog/LANDSAT_LC08_C02_T1_L2
#  Landsat 9: https://developers.google.com/earth-engine/datasets/catalog/LANDSAT_LC09_C02_T1_L2
# NOTE - the difference between the 4/5/7 and 8/9 series is addition of Cirrus information in 8/9
# if maskMask is set to ee.Image(0), then the mask will not change
def landsat_cloud_mask_89(image, mask_mask=1):
    # Bit-shift to get bits which need to be zero-valued to avoid masking
    dilated_cloud = (1 << 1)         
    cirrus = (1 << 2)
    cloud = (1 << 3)                 
    cloud_shadow = (1 << 4)          
    qa = image.select('QA_PIXEL');  
    mask = (qa.bitwiseAnd(dilated_cloud)
              .Or(qa.bitwiseAnd(cirrus))
              .Or(qa.bitwiseAnd(cloud_shadow))
              .Or(qa.bitwiseAnd(cloud))
              .multiply(mask_mask))
    return image.updateMask(mask.Not())


# *******************************************************************************************************
# ********************************** Surface Reflectance Scale Factors **********************************
# *******************************************************************************************************
def landsat_rescale(image):
    image_sr_offset = image.select('SR_B.*').multiply(0.0000275).add(-0.2);      
    image_st_offset = image.select('ST_B.*').multiply(0.00341802).add(149);      
    image_sr_qa = (image_sr_offset.addBands(image_st_offset)
                                  .addBands(image.select('QA_PIXEL')))
    return (image_sr_qa.set('system:index', image.get('system:index'))
                       .set('system:time_start', image.get('system:time_start')))


# *******************************************************************************************************
# ******************************************** NDVI Function ********************************************
# *******************************************************************************************************

# For Landsat 4,5,7, there is one band numbering convention (NIR->4, Red->3)
def get_ndvi_457(image):
    return image.addBands(image.normalizedDifference(["SR_B4","SR_B3"]).rename('NDVI'))

# For Landsat 8, there is a different band numbering convention (NIR->5, Red->4)
def get_ndvi_89(image):
    return image.addBands(image.normalizedDifference(["SR_B5","SR_B4"]).rename('NDVI'))


# *******************************************************************************************************
# ******************************************** EVI Function ********************************************
# *******************************************************************************************************

# For Landsat 4,5,7, there is one band numbering convention (NIR->4, Red->3, etc.)
def get_evi_457(image):
    # EVI = G * (NIR - R) / (NIR + 6*R - 7.5*B + 1)
    return (image.addBands(ee.Image(2.5).multiply((image.select('SR_B4').subtract(image.select('SR_B3')))
                 .divide((image.select('SR_B4')
                 .add(ee.Image(6).multiply(image.select('SR_B3')))
                 .subtract(ee.Image(7.5).multiply(image.select('SR_B1')))
                 .add(ee.Image(1)))))
                 .rename('EVI')))

# For Landsat 8, there is a different band numbering convention (NIR->5, Red->4, etc.)
def get_evi_89(image):
    return (image.addBands(ee.Image(2.5).multiply((image.select('SR_B5').subtract(image.select('SR_B4')))
                 .divide((image.select('SR_B5')
                 .add(ee.Image(6).multiply(image.select('SR_B4')))
                 .subtract(ee.Image(7.5).multiply(image.select('SR_B2')))
                 .add(ee.Image(1)))))
                 .rename('EVI')))

# *******************************************************************************************************
# ******************************************** NIRv Function ********************************************
# *******************************************************************************************************

# For Landsat 4,5,7, there is one band numbering convention (NIR->4, Red->3)
def get_nirv_457(image):
    return image.addBands(image.normalizedDifference(["SR_B4","SR_B3"]).multiply(image.select("SR_B4")).rename('NIRv'))

# For Landsat 8, there is a different band numbering convention (NIR->5, Red->4)
def get_nirv_89(image):
    return image.addBands(image.normalizedDifference(["SR_B5","SR_B4"]).multiply(image.select("SR_B5")).rename('NIRv'))

# *******************************************************************************************************
# ******************************************** SAVI Function ********************************************
# *******************************************************************************************************

# For Landsat 4,5,7, there is one band numbering convention (NIR->4, Red->3)
def get_savi_457(image):
    return image.addBands(ee.Image(1.5) \
    .multiply(image.select('SR_B4').subtract(image.select('SR_B3'))) \
    .divide(ee.Image(0.5).add(image.select('SR_B4')).add(image.select('SR_B3'))) \
    .rename('SAVI'))

# For Landsat 8, there is a different band numbering convention (NIR->5, Red->4)
def get_savi_89(image):
    return image.addBands(ee.Image(1.5) \
    .multiply(image.select('SR_B5').subtract(image.select('SR_B4'))) \
    .divide(ee.Image(0.5).add(image.select('SR_B5')).add(image.select('SR_B4'))) \
    .rename('SAVI'))


# *******************************************************************************************************
# ******************************************** NDSVI Function ********************************************
# *******************************************************************************************************

# For Landsat 4,5,7, there is one band numbering convention (NIR->4, Red->3)
def get_ndsvi_457(image):
        return (image.addBands((ee.Image(1).subtract(image.select("SR_B7").divide(image.select("SR_B5")))
                 .multiply(image.select("SR_B3")).divide(image.select("SR_B4"))).rename("NDSVI")))

# For Landsat 8, there is a different band numbering convention (NIR->5, Red->4)
def get_ndsvi_457(image):
    return (image.addBands((ee.Image(1).subtract(image.select("SR_B7").divide(image.select("SR_B6")))
                 .multiply(image.select("SR_B4")).divide(image.select("SR_B5"))).rename("NDSVI")))

# *******************************************************************************************************
# ****************************** Temperature Difference (Surface vs. Air) *******************************
# *******************************************************************************************************

# We implement two versions here. Using Daymet is preferred, but the timeseries is less extensive, so NOA NCEI can be used in cases without Daymet coverage

# Get difference between surface temperature and local air temperature
def get_temp_difference(img):
    # Load DAYMET climatology data to extract
    daymet_climatology = ee.ImageCollection("NASA/ORNL/DAYMET_V4")
    date = img.date()
    # Find the closest climate data to the target imagery
    # Get all climate data within 3 days of target
    nearby_climates = daymet_climatology.filterDate(date.advance(-3,'days'), date.advance(3,'days'))
    # For each climate value, get time difference vs. target image
    def func_npv(clim_img):
        return clim_img.set(
            'time_difference',
            ee.Number(clim_img.get('system:time_start')).subtract(img.get('system:time_start')).abs()
        )
    nearby_climates = nearby_climates.map(func_npv)
    # Sort by time difference (ascending)
    nearby_climates.sort('time_difference')
    # Add climate data (air temperature and relative land vs. air temperature) to input image
    new_img = (img.addBands(nearby_climates.first().select('tmax'))
                  .addBands(img.select('ST')
                  .subtract(nearby_climates.first().select('tmax')).subtract(273.15)
                  .rename('ST_air_diff')))
    return(new_img)

# Get difference between surface temperature and local air temperature
def get_temp_differnce_noancei(img):    
    # Load NCEI climatology data to extract
    noa_ncei_climatology = ee.ImageCollection("projects/climate-engine-pro/assets/noaa-ncei-nclimgrid/daily")
    date = img.date()
    # Find the closest climate data to the target imagery
    # Get all climate data within 3 days of target
    nearby_climates = noa_ncei_climatology.filterDate(date.advance(-3,'days'), date.advance(3,'days'))
    # For each climate value, get time difference vs. target image
    def func_wls(clim_img):
            return clim_img.set(
                'time_difference',
                ee.Number(clim_img.get('system:time_start')).subtract(img.get('system:time_start')).abs()
            )
    nearby_climates = nearby_climates.map(func_wls)
    # Sort by time difference (ascending)
    nearby_climates.sort('time_difference')
    # Add climate data (air temperature and relative land vs. air temperature) to input image
    newimg = img.addBands(nearby_climates.first().select('tmax')) \
    .addBands(img.select('ST') \
    .subtract(nearby_climates.first().select('tmax')).subtract(273.15) \
    .rename('ST_air_diff'))
    return(newimg)


# *******************************************************************************************************
# *************************************** Common Band Extraction ****************************************
# *******************************************************************************************************

# Pull out just the common bands among Landsat 5 - 9
# Rename SR bands to Landsat 5 basis
# Rename ST_10 and ST_6 to just ST
def get_common_band_names_89(image):
    return(image.select([
             'SR_B2',
             'SR_B3',
             'SR_B4',
             'SR_B5',
             'SR_B6',
             'SR_B7',
             'ST_B10',
             'SR_B1']
            ).rename(['SR_B1',
                      'SR_B2',
                      'SR_B3',
                      'SR_B4',
                      'SR_B5',
                      'SR_B7',
                      'ST',
                      'SR_coastal'
                     ]))

def get_common_band_names_457(image):
    return(image.select([
             'SR_B1',
             'SR_B2',
             'SR_B3',
             'SR_B4',
             'SR_B5',
             'SR_B7',
             'ST_B6']
            ).rename([
                'SR_B1',
                'SR_B2',
                'SR_B3',
                'SR_B4',
                'SR_B5',
                'SR_B7',
                'ST'
            ]))

def remove_coastal_band(image):
    return image.select(
                image.bandNames().filter(ee.Filter.stringEndsWith('item', 'SR_coastal').Not())
            )

