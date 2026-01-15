import ee 

# ********************* NOTE ON USAGE *********************
# This file contains two functions which are very similar. They require the same inputs and generate
#    the same outputs, but the implementation is slightly different. Both functions use the 'iterate'
#    utility, which requires converting the input imagery into a List format (memory-intensive).
#    The trade-off is in how much imagery is converted to list, and in how much post-processing is needed.
# cloudTemporalFilter
#    This function converts ALL the data in pheno_in into the list for the iterate function. For files with
#    lots of bands, this can be problematic and cause memory issues.
# cloudTemporalFilterLowRAM
#    This function converts only the needed bands (DOY and band_name) into a list, which can reduce memory
#    usage for inputs with lots of bands. HOWEVER, to re-apply the mask to all bands afterwards, we have to
#    include a post-processing step which maps through all images. For some collections, THIS can cause
#    memory errors.
# In my usage so far, I have found Landsat and Sentinel-2 to work better with the 'cloudTemporalFilter' and
#    MODIS to work better with 'cloudTemporalFilterLowRAM.' Users should try both to see which performs
#    better for their use case.
# A longer-term goal of mine is to rework these functions so that they no longer use the iterate() utility
#    at all, which may help get around these problems and improve scalability.

# Search for and filter out probable cloudy scenes
#   At each pixel, compares each date to its two closest unmasked neighbor dates
#   Fits a linear regression (BAND_NAME ~ Time) to the neighbor dates
#   If the target date is much different in NDVI than the neighbors, rule it a cloud
#   Uses 2 input values to determine threshold differences beyond which a pixel is masked
def cloud_temporal_filter(pheno_in, band_name, threshold_low, threshold_high, num_padding_scenes = 0):
    
    # ----------------- Buffer Timeseries -----------------
    
    # Extra scenes, from start of timeseries, to be folded onto end as a buffer
    def buffer_start(img):
        img = img.addBands(srcImg=(img.select('DOY')
                                      .add(ee.Image(365))
                                      .rename('DOY')),
                           overwrite=True)
        img = img.set('system:time_start', 
                      ee.Number(img.get('system:time_start')).add(ee.Number(31536000000)))
        return img
    first_few_scenes = pheno_in.limit(num_padding_scenes, 'system:time_start', True).map(buffer_start)

    # Extra scenes, from start of timeseries, to be folded onto end as a buffer
    def buffer_end(img):
        img = img.addBands(srcImg=(img.select('DOY')
                                      .subtract(ee.Image(365))
                                      .rename('DOY')),
                           overwrite=True)
        img = img.set('system:time_start', 
                      ee.Number(img.get('system:time_start')).subtract(ee.Number(31536000000)))
        return img
    last_few_scenes = pheno_in.limit(num_padding_scenes, 'system:time_start', False).map(buffer_end)

    # Combine leading and ending buffers with original dataset 
    def subset_bands_and_update_mask(img):
        return (img.updateMask(img.select(band_name).mask())
                   .select([band_name, 'DOY']))
    pheno_buffered = (last_few_scenes.merge(pheno_in.merge(first_few_scenes))
                                     .sort('system:time_start')
                                     .map(subset_bands_and_update_mask))

    # ----------------- Find Preceding and Following Images for Each Image -----------------    
    
    # The first scene is set to be 0
    first = ee.List([
        ee.Image(0).addBands(ee.Image(0)).set('system:time_start', 0).select([0, 1], [band_name, 'DOY']).float()
    ])

    # Function which, at each index of an input ImageCollection, retains the most recent PREVIOUS image which was not masked
    def get_previous_unmasked_scene(image, running_list):
            # Get the value which WAS the previous value, at the previous index (ie, n-2)
            previous = ee.Image(ee.List(running_list).get(-1))
            # Get a boolean image, 0 where current image is masked, 1 where it is not masked
            current_unmasked = image.select(band_name).mask()
            # Set the previous value for this scene to the previous image, OR
            #   to the value which WAS previous for that image if the previous image is masked
            return (ee.List(running_list).add(ee.Image(image.select([band_name,'DOY']).unmask(None, False).multiply(current_unmasked)
                                                .add(previous.unmask(None, False).multiply(current_unmasked.Not())))))
    # Collection where each value is the most recent scene in the collection which was not masked
    previous_unmasked_scene = ee.List(pheno_buffered.iterate(get_previous_unmasked_scene, first))
    # Get list of all indices which should be kept in final version (trimming off buffers on start/end)
    # Current length is=1 (filler) + 2*num_padding_scenes + pheno_in.size()
    indices_to_retain = ee.List.sequence(
        start=ee.Number(num_padding_scenes).add(ee.Number(1)).int(),
        end=ee.Number(previous_unmasked_scene.size()).subtract(ee.Number(num_padding_scenes)).subtract(ee.Number(1)).int()
    )
    # Remove extra buffered values off start and end
    def get_retained_previous_image(ind):
        return previous_unmasked_scene.get(ee.Number(ind).subtract(1))
    previous_unmasked_scene = indices_to_retain.map(get_retained_previous_image)
    
    # The last scene is set to be 0
    last = ee.List([
        ee.Image(0).addBands(ee.Image(365)).set('system:time_start', 0).select([0, 1], [band_name, 'DOY']).float()
    ])
    # Collection where each value is the next scene in the collection which was not masked
    next_unmasked_scene = ee.List(pheno_buffered.sort('system:time_start', False).iterate(get_previous_unmasked_scene, last))
    # Remove extra buffered values off start and end
    def get_retained_next_image(ind):
        return next_unmasked_scene.get((ee.Number(ind).multiply(ee.Number(-1))).subtract(1))
    next_unmasked_scene = indices_to_retain.map(get_retained_next_image)

    
    # ----------------- Create list of actual target images in collection ----------------- 

    # Convert input phenology data to a list
    def get_list_iterator(image, image_list):
            return ee.List(image_list).add(image)
    
    first_input = ee.List([
        pheno_in.first().mask(1).multiply(0)
    ])
    # NOTE the first value in this is empty - make sure to offset index by one
    pheno_in_list = ee.List(pheno_in.sort('system:time_start').iterate(get_list_iterator, first_input))
    
    # Now, for each target image in the collection, get the predicted value from a linear regression between the previous and next data points
    def check_value(ind):
            # Get the previous and following images
            previous_scene = ee.Image(previous_unmasked_scene.get(ee.Number(ind)))
            next_scene = ee.Image(next_unmasked_scene.get(ee.Number(ind)))
            # Get the actual observed value
            current_scene = ee.Image(pheno_in_list.get(ee.Number(ind).add(ee.Number(1))))
            # Get change between previous and following image values
            neighbor_change = next_scene.subtract(previous_scene)
            # Get the slope of change
            slope = neighbor_change.select(band_name).divide(neighbor_change.select('DOY'))
            # Get the expected value for the target date
            expected_value = ((current_scene.select('DOY').subtract(previous_scene.select('DOY')))
                                            .multiply(slope)
                                            .add(previous_scene.select(band_name)))
            # If the difference in expected and actual value is outside allowed bounds, mask the image
            difference = current_scene.select(band_name).subtract(expected_value)
            within_low_threshold = ee.Image(threshold_low).lt(difference)
            within_high_threshold = ee.Image(threshold_high).gt(difference)
            allowed = within_low_threshold.multiply(within_high_threshold)
            return (current_scene.updateMask(allowed)
                                 .addBands(difference.rename('reg_diff'))
                                 .addBands(expected_value.rename('reg_expected'))
                                 .addBands(previous_scene.select(band_name).rename('previous'))
                                 .addBands(previous_scene.select('DOY').rename('previous_DOY'))
                                 .addBands(next_scene.select(band_name).rename('next'))
                                 .addBands(next_scene.select('DOY').rename('next_DOY')))

    target_sequence = ee.List.sequence(0,ee.Number(pheno_in.size()).subtract(ee.Number(1)))

    output_collection = ee.ImageCollection(target_sequence.map(check_value))

    return output_collection


'''


# Search for and filter out probable cloudy scenes
#   At each pixel, compares each date to its two closest unmasked neighbor dates
#   Fits a linear regression (BAND_NAME ~ Time) to the neighbor dates
#   If the target date is much different in NDVI than the neighbors, rule it a cloud
#   Uses 2 input values to determine threshold differences beyond which a pixel is masked
{
    # Set default value for num_padding_scenes to 0
if (num_padding_scenes == undefined or num_padding_scenes == None) =

    # Subset input collection to hold ONLY the bands DOY and target band for filtering
    #   This prevents loading lots of unnecessary data into memory
    pheno_in_subset = pheno_in.select(['DOY',band_name])

    # Extra scenes, from start of timeseries, folded onto end as a buffer
    first_few_scenes = pheno_in_subset.limit(num_padding_scenes, 'system =time_start', True)

    def func_pjg(img) =
            img = img.addBands(
                    srcImg=img.select('DOY') \
                    .add(ee.Image(365)) \
                    .rename('DOY'),
                    overwrite=True})
            img = img.set('system =time_start', ee.Number(img.get('system =time_start')).add(ee.Number(31536000000)))
            return(img) \
    .map(func_pjg)








    # Extra scenes, from end of timeseries, folded onto start as a buffer
    last_few_scenes = pheno_in_subset.limit(num_padding_scenes, 'system =time_start', False)

    def func_omh(img) =
            img = img.addBands(
                    srcImg=img.select('DOY') \
                    .subtract(ee.Image(365)) \
                    .rename('DOY'),
                    overwrite=True})
            img = img.set('system =time_start', ee.Number(img.get('system =time_start')).subtract(ee.Number(31536000000)))
            return(img) \
    .map(func_omh)








    # Add buffers onto start and end of initial dataset
    pheno_buffered = last_few_scenes.merge(pheno_in_subset.merge(first_few_scenes)) \
    .sort('system =time_start')

    def func_rxs(img) =
            return(img.updateMask(img.select(band_name).mask()) \
            .select([band_name,'DOY'])) \
    .map(func_rxs)




    # The first scene is set to be 0
    first = ee.List([
    ee.Image(0).addBands(ee.Image(0)).set('system =time_start', 0).select([0, 1], [band_name, 'DOY']).float()
    ])

    # Function which, at each index of an input ImageCollection, retains the most recent PREVIOUS image which was not masked

    def func_mpp(image, list) =
            # Get the value which WAS the previous value, at the previous index (ie, n-2)
            previous = ee.Image(ee.List(list).get(-1))
            # Get a boolean image, 0 where current image is masked, 1 where it is not masked
            current_unmasked = image.select(band_name).mask()
            # Set the previous value for this scene to the previous image, OR
            #   to the value which WAS previous for that image if the previous image is masked
            return ee.List(list).add(ee.Image(image.select([band_name,'DOY']).unmask(None, False).multiply(current_unmasked) \
            .add(previous.unmask(None, False).multiply(current_unmasked.Not()))))

    getLastUnmaskedScene = func_mpp










    # Collection where each value is the most recent scene in the collection which was not masked
    last_unmasked_scene = ee.List(pheno_buffered.iterate(getLastUnmaskedScene, first))

    # Get list of all indices which should be kept in final version (trimming off buffers on start/end)
    # Current length is=1 (filler) + 2*num_padding_scenes + pheno_in.size()
    indices_to_retain = ee.List.sequence(
        start= ee.Number(num_padding_scenes).add(ee.Number(1)).int(),
        end=   ee.Number(last_unmasked_scene.size()).subtract(ee.Number(num_padding_scenes)).subtract(ee.Number(1)).int()
        )

    # Remove extra buffered values off start and end

    def func_nqj(ind) =
            return last_unmasked_scene.get(ee.Number(ind).subtract(1))

    last_unmasked_scene = indices_to_retain.map(func_nqj)



    # The last scene is set to be 0
    last = ee.List([
    ee.Image(0).addBands(ee.Image(365)).set('system =time_start', 0).select([0, 1], [band_name, 'DOY']).float()
    ])

    # Collection where each value is the next scene in the collection which was not masked
    next_unmasked_scene = ee.List(pheno_buffered.sort('system =time_start', False).iterate(getLastUnmaskedScene, last))
    # Remove extra buffered values off start and end

    def func_vjc(ind) =
            return next_unmasked_scene.get(ee.Number(ind).multiply(ee.Number(-1)).subtract(1))

    next_unmasked_scene = indices_to_retain.map(func_vjc)



    # Convert input phenology data to a list
    {
            return ee.List(list).add(image)
        }
    first_input = ee.List([
    pheno_in_subset.first().mask(1).multiply(0)
    ])
    # NOTE the first value in this is empty - make sure to offset index by one
    pheno_in_list = ee.List(pheno_in_subset.sort('system =time_start').iterate(getListIterator, first_input))

    # Now, for each target image in the collection, get the predicted value from a linear regression between the previous and next data points
    {
            # Get the previous and following images
            previous_scene = ee.Image(last_unmasked_scene.get(ee.Number(ind)))
            next_scene = ee.Image(next_unmasked_scene.get(ee.Number(ind)))
            # Get the actual observed value
            current_scene = ee.Image(pheno_in_list.get(ee.Number(ind).add(ee.Number(1))))
            # Get change between previous and following image values
            neighbor_change = next_scene.subtract(previous_scene)
            # Get the slope of change
            slope = neighbor_change.select(band_name).divide(neighbor_change.select('DOY'))
            # Get the expected value for the target date
            expected_value = (current_scene.select('DOY').subtract(previous_scene.select('DOY'))) \
            .multiply(slope) \
            .add(previous_scene.select(band_name))                                
            # If the difference in expected and actual value is outside allowed bounds, mask the image
            difference = current_scene.select(band_name).subtract(expected_value)
            within_low_threshold = ee.Image(threshold_low).lt(difference)
            within_high_threshold = ee.Image(threshold_high).gt(difference)
            allowed = within_low_threshold.multiply(within_high_threshold)
            return(current_scene.updateMask(allowed) \
            .addBands(difference.rename('reg_diff')) \
            .addBands(expected_value.rename('reg_expected')) \
            .addBands(previous_scene.select(band_name).rename('previous')) \
            .addBands(previous_scene.select('DOY').rename('previous_DOY')) \
            .addBands(next_scene.select(band_name).rename('next')) \
            .addBands(next_scene.select('DOY').rename('next_DOY')))
        }
    target_sequence = ee.List.sequence(0,ee.Number(pheno_in_subset.size()).subtract(ee.Number(1)))

    # List of bands in input OTHER than target band for filtering
    input_bands = pheno_in.filterBounds(ee.Geometry.Point([-115.547219,33.0889])).first().bandNames()
    non_target_bands = input_bands.remove(band_name)

    # Masked image with only target band (input parameter band_name)
    mask_collection = ee.ImageCollection(target_sequence.map(checkValue)).select(band_name)

    # Try this with a map() and filter based on scene ID...

    def func_hjd(img) =
            mask = ee.Image(mask_collection.filterBounds(img.geometry()).filter(ee.Filter.eq('system =index', img.get('system =index'))).first())
            return img.updateMask(mask.mask())

    return pheno_in.map(func_hjd)



'''