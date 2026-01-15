import ee 
from . import basic_processing

# Function to add day of year as a band to a Landsat image
def add_doy(image):
    # Overall date of image
    date = ee.Date(image.get('system:time_start'))
    # Year in which image occurred
    year = ee.Number(date.difference(ee.Date('1970-01-01'),'year')).floor()
    # Day of that year on which image occurred
    day_of_year = date.difference(ee.Date('1970-01-01'), 'day').subtract(year.multiply(ee.Number(365.25)))
    # Add DOY to image as a new band
    return image.addBands(ee.Image(day_of_year).rename('DOY').float())
    
# Function to add day of year as a band to a Landsat image
# Units are in years (e.g. 2024 is the first second of Jan 1, 2024)
def add_time(image):
  # Time in years (originally in milliseconds since Jan 1, 1970) of image
  time = ee.Number(image.get('system:time_start')).divide(365.25*24*3600*1000).add(1970)
  # Add to image and return
  return image.addBands(ee.Image.constant(time).toFloat().rename('time'))

# Returns a collection with all images from Landsat 4, 5, 7, 8, and 9 within a given time period
# Applies reflectance scaling factors and the cloud filter distributed with the product
def get_collection(date_start, date_stop, mask_clouds=1):

    def landsat_cloud_mask_457(image): 
        return basic_processing.landsat_cloud_mask_457(image, mask_clouds)
    def landsat_cloud_mask_89(image): 
        return basic_processing.landsat_cloud_mask_89(image, mask_clouds)
    
    collection_l4 = (ee.ImageCollection("LANDSAT/LT04/C02/T1_L2")
                       .filterDate(date_start, date_stop)
                       .map(basic_processing.landsat_rescale)
                       .map(landsat_cloud_mask_457)
                       .map(basic_processing.get_common_band_names_457))
    collection_l5 = (ee.ImageCollection("LANDSAT/LT05/C02/T1_L2")
                       .filterDate(date_start, date_stop)
                       .map(basic_processing.landsat_rescale)
                       .map(landsat_cloud_mask_457)
                       .map(basic_processing.get_common_band_names_457))
    collection_l7 = (ee.ImageCollection("LANDSAT/LE07/C02/T1_L2")
                       .filterDate(date_start, date_stop)
                       .map(basic_processing.landsat_rescale)
                       .map(landsat_cloud_mask_457)
                       .map(basic_processing.get_common_band_names_457))
    collection_l8 = (ee.ImageCollection("LANDSAT/LC08/C02/T1_L2")
                       .filterDate(date_start, date_stop)
                       .map(basic_processing.landsat_rescale)
                       .map(landsat_cloud_mask_89)
                       .map(basic_processing.get_common_band_names_89))
    collection_l9 = (ee.ImageCollection("LANDSAT/LC09/C02/T1_L2")
                       .filterDate(date_start, date_stop)
                       .map(basic_processing.landsat_rescale)
                       .map(landsat_cloud_mask_89)
                       .map(basic_processing.get_common_band_names_89))

    all_scenes = collection_l9.merge(collection_l8.merge(collection_l7.merge(collection_l5.merge(collection_l4))))

    return all_scenes

'''
# Returns a collection with all images from Landsat 4, 5, 8, and 9 within a given time period (NO Landsat 7)
# Applies reflectance scaling factors and the cloud filter distributed with the product
{
if (maskClouds == undefined or maskClouds == None):
    collection_L4 = ee.ImageCollection("LANDSAT/LT04/C02/T1_L2") \
    .filterDate(date_start, date_stop) \
    .map(lb.landsatReflectanceOffsets)

    def func_giw(image) return lb.cloudMaskC2L457(image, ee.Image(maskClouds))}: \
    .map(function(image) {return lb.cloudMaskC2L457(image, ee.Image(maskClouds))} \
    .map(func_giw) \
    .map(lb.getCommonBandNamesL457)
    collection_L5 = ee.ImageCollection("LANDSAT/LT05/C02/T1_L2") \
    .filterDate(date_start, date_stop) \
    .map(lb.landsatReflectanceOffsets)

    def func_cbd(image) return lb.cloudMaskC2L457(image, ee.Image(maskClouds))}: \
    .map(function(image) {return lb.cloudMaskC2L457(image, ee.Image(maskClouds))} \
    .map(func_cbd) \
    .map(lb.getCommonBandNamesL457)
    collection_L8 = ee.ImageCollection("LANDSAT/LC08/C02/T1_L2") \
    .filterDate(date_start, date_stop) \
    .map(lb.landsatReflectanceOffsets)

    def func_ooh(image) return lb.cloudMaskC2L89(image, ee.Image(maskClouds))}: \
    .map(function(image) {return lb.cloudMaskC2L89(image, ee.Image(maskClouds))} \
    .map(func_ooh) \
    .map(lb.getCommonBandNamesL89)
    collection_L9 = ee.ImageCollection("LANDSAT/LC09/C02/T1_L2") \
    .filterDate(date_start, date_stop) \
    .map(lb.landsatReflectanceOffsets)

    def func_wqu(image) return lb.cloudMaskC2L89(image, ee.Image(maskClouds))}: \
    .map(function(image) {return lb.cloudMaskC2L89(image, ee.Image(maskClouds))} \
    .map(func_wqu) \
    .map(lb.getCommonBandNamesL89)

    all_scenes = collection_L9.merge(collection_L8.merge(collection_L5.merge(collection_L4)))

    return all_scenes

'''

# Returns a collection with all images from Landsat 8 and 9 OLI within a given time period (NO Landsat 7)
# Applies reflectance scaling factors and the cloud filter distributed with the product
def get_collection_89(date_start, date_stop, mask_clouds=1):
    
    def landsat_cloud_mask_89(image): 
        return basic_processing.landsat_cloud_mask_89(image, mask_clouds)

    collection_l8 = (ee.ImageCollection("LANDSAT/LC08/C02/T1_L2")
                       .filterDate(date_start, date_stop)
                       .map(basic_processing.landsat_rescale)
                       .map(landsat_cloud_mask_89)
                       .map(basic_processing.get_common_band_names_89))
    collection_l9 = (ee.ImageCollection("LANDSAT/LC09/C02/T1_L2")
                       .filterDate(date_start, date_stop)
                       .map(basic_processing.landsat_rescale)
                       .map(landsat_cloud_mask_89)
                       .map(basic_processing.get_common_band_names_89))

    all_scenes = collection_l9.merge(collection_l8)

    return all_scenes


'''

# Returns a collection with all images from Landsat 8 and 9 OLI within a given time period (NO Landsat 7)
# Applies reflectance scaling factors and the cloud filter distributed with the product
{
if (maskClouds == undefined or maskClouds == None):
    collection_L7 = ee.ImageCollection("LANDSAT/LE07/C02/T1_L2") \
    .filterDate(date_start, date_stop) \
    .map(lb.landsatReflectanceOffsets)

    def func_rkf(image) return lb.cloudMaskC2L457(image, ee.Image(maskClouds))}: \
    .map(function(image) {return lb.cloudMaskC2L457(image, ee.Image(maskClouds))} \
    .map(func_rkf) \
    .map(lb.getCommonBandNamesL457)
    return collection_L7

'''



# Returns a collection with all images from Landsat 4, 5, and 7 within a given time period
# Applies reflectance scaling factors and the cloud filter distributed with the product
def get_collection_457(date_start, date_stop, mask_clouds=1):

    def landsat_cloud_mask_457(image): 
        return basic_processing.landsat_cloud_mask_457(image, mask_clouds)
    
    collection_l4 = (ee.ImageCollection("LANDSAT/LT04/C02/T1_L2")
                       .filterDate(date_start, date_stop)
                       .map(basic_processing.landsat_rescale)
                       .map(landsat_cloud_mask_457)
                       .map(basic_processing.get_common_band_names_457))
    collection_l5 = (ee.ImageCollection("LANDSAT/LT05/C02/T1_L2")
                       .filterDate(date_start, date_stop)
                       .map(basic_processing.landsat_rescale)
                       .map(landsat_cloud_mask_457)
                       .map(basic_processing.get_common_band_names_457))
    collection_l7 = (ee.ImageCollection("LANDSAT/LE07/C02/T1_L2")
                       .filterDate(date_start, date_stop)
                       .map(basic_processing.landsat_rescale)
                       .map(landsat_cloud_mask_457)
                       .map(basic_processing.get_common_band_names_457))

    all_scenes = collection_l7.merge(collection_l5.merge(collection_l4))

    return all_scenes
'''

# Function to add day of year as a band to a MODIS image

def func_eoq(image):
    # Overall date of image
    date = ee.Date(image.get('system:time_start'))
    # Year in which image occurred
    year = ee.Number(date.difference(ee.Date('1970-01-01'),'year')).floor()
    # Day of that year on which image occurred
    day_of_year = date.difference(ee.Date('1970-01-01'), 'day').subtract(year.multiply(ee.Number(365.25)))
    # Add DOY to image as a new band
    return image.addBands(ee.Image(day_of_year).rename('DOY').float())

addDOY = func_eoq



'''