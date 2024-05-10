def extract_CROCO(file_name, file_base = '/media/data/CHPC_SBUS_3km/',  lat_min = -36., lat_max = -28., lon_min = 15.,  lon_max = 20.):
    import numpy as np
    from netCDF4 import Dataset
    nc_file = file_base + file_name
    root = Dataset(nc_file)
    lat = root.variables['lat_rho'][:]
    lon = root.variables['lon_rho'][:]
    
    def find_nearest(array, value):
        array = np.asarray(array)
        idx = (np.abs(array - value)).argmin()
        return array[idx], idx
    
    lat_min_index = find_nearest(lat[:,0],lat_min)[1]
    lat_max_index = find_nearest(lat[:,0],lat_max)[1]
    lon_min_index = find_nearest(lon[0,:],lon_min)[1]
    lon_max_index = find_nearest(lon[0,:],lon_max)[1]
    
    
    lon_CROCO = lon[0,lon_min_index:lon_max_index]
    lat_CROCO = lat[lat_min_index:lat_max_index,0]
    sst_CROCO = root.variables['temp'][0, 59,lat_min_index:lat_max_index, lon_min_index:lon_max_index]
    sst_CROCO = np.squeeze(sst_CROCO)
    root.close()
    return sst_CROCO, lon_CROCO, lat_CROCO

def myCanny(myData, myMask, sigma = 10., lower = .8, upper = .9, use_quantiles = True):
    # because of the way masks operate,  if you read in sst using netcdf4,
    # then the mask to use is ~sst.mask
    import numpy.ma as ma
    edges, y_gradient, x_gradient, magnitude = canny1(myData, sigma = sigma, mask = myMask,low_threshold = lower, high_threshold = upper,use_quantiles = use_quantiles)
    x_gradient = ma.array(x_gradient, mask = myData.mask)
    y_gradient = ma.array(y_gradient, mask = myData.mask)
    magnitude = ma.array(magnitude, mask = myData.mask)
    return edges, x_gradient, y_gradient, magnitude

def create_canny_nc(file_year, file_month, file_day, base_dir = '/media/data/DIAGNOSTICS/',  lat_min = -36., lat_max = -28., lon_min = 15.,  lon_max = 20.):
    from netCDF4 import Dataset, num2date, date2num
    import numpy as np
    import numpy.ma as ma
    c_file_year = str(file_year)
    c_file_month = str(file_month).rjust(2,'0')
    c_file_day = str(file_day).rjust(2,'0')
    file_name = base_dir + 'Canny_Front_' + c_file_year + c_file_month + c_file_day +  '.nc'
    ncfile  = Dataset(file_name, 'w', format = 'NETCDF4')
    lat_diff = lat_max - lat_min
    latsdim = (lat_diff * 100) + 1
    lats = lat_min + (np.arange(0, latsdim) * 0.01)
    lon_diff = lon_max - lon_min
    lonsdim = (lon_diff * 100) + 1
    lons = lon_min + (np.arange(0, lonsdim) * 0.01)
    #Create Dimensions
    timedim = ncfile.createDimension('time', None)
    latdim = ncfile.createDimension('lat', latsdim)
    londim = ncfile.createDimension('lon', lonsdim)
    altdim = ncfile.createDimension('altitude', 1)
    #Create Variables
    LatLon_Projection = ncfile.createVariable('LatLon_Projection', 'i4')
    time = ncfile.createVariable('time', 'f8', ('time'), zlib = True, complevel = 2)
    altitude = ncfile.createVariable('altitude', 'f4', ('altitude'))
    latitude = ncfile.createVariable('lat', 'f4', ('lat'), zlib = True, complevel = 2)
    longitude = ncfile.createVariable('lon', 'f4', ('lon'), zlib = True, complevel = 2)
    edges = ncfile.createVariable('edges', 'f4', ('time', 'altitude', 'lat', 'lon'), fill_value = -9999.0, zlib = True, complevel = 2)
    x_gradient = ncfile.createVariable('x_gradient', 'f4', ('time', 'altitude', 'lat', 'lon'), fill_value = -9999.0, zlib = True, complevel = 2)
    y_gradient = ncfile.createVariable('y_gradient', 'f4', ('time', 'altitude', 'lat', 'lon'), fill_value = -9999.0, zlib = True, complevel = 2)
    magnitude_gradient = ncfile.createVariable('magnitude_gradient', 'f4', ('time', 'altitude', 'lat', 'lon'), fill_value = -9999.0, zlib = True, complevel = 2)
    # int LatLon_Projection ;
    LatLon_Projection.grid_mapping_name = "latitude_longitude"
    LatLon_Projection.earth_radius = 6367470.
    #float lat(lat) ;
    latitude._CoordinateAxisType = "Lat"
    junk = (lat_min, lat_max)
    latitude.actual_range = junk
    latitude.axis = "Y"
    latitude.grid_mapping = "Equidistant Cylindrical"
    latitude.ioos_category = "Location"
    latitude.long_name = "Latitude"
    latitude.reference_datum = "geographical coordinates, WGS84 projection"
    latitude.standard_name = "latitude"
    latitude.units = "degrees_north"
    latitude.valid_max = lat_max
    latitude.valid_min = lat_min
    #float lon(lon) ;
    longitude._CoordinateAxisType = "Lon"
    junk = (lon_min, lon_max)
    longitude.actual_range = junk
    longitude.axis = "X"
    longitude.grid_mapping = "Equidistant Cylindrical"
    longitude.ioos_category = "Location"
    longitude.long_name = "Longitude"
    longitude.reference_datum = "geographical coordinates, WGS84 projection"
    longitude.standard_name = "longitude"
    longitude.units = "degrees_east"
    longitude.valid_max = lon_max
    longitude.valid_min = lon_min
    #float altitude(altitude) ;
    altitude.units = "m"
    altitude.long_name = "Specified height level above ground"
    altitude.standard_name = "altitude"
    altitude.positive = "up"
    altitude.axis = "Z"
    #double time(time) ;
    time._CoordinateAxisType = "Time"
    junk = ()
    time.actual_range = junk
    time.axis = "T"
    time.calendar = "Gregorian"
    time.ioos_category = "Time"
    time.long_name = "Time"
    time.units = "Hour since 1970-01-01T00:00:00Z"
    time.standard_name = "time"
    #float edges(time, altitude, lat, lon) ;
    edges.long_name = "Frontal Edge"
    edges.missing_value = -9999.
    edges.grid_mapping = "LatLon_Projection"
    edges.coordinates = "time altitude lat lon "
    #float x_gradient(time, altitude, lat, lon) ;
    x_gradient.long_name = "East-West Gradient of SST"
    x_gradient.missing_value = -9999.
    x_gradient.grid_mapping = "LatLon_Projection"
    x_gradient.coordinates = "time altitude lat lon "
    # float y_gradient(time, altitude, lat, lon) ;
    y_gradient.long_name = "North-South Gradient of SST"
    y_gradient.missing_value = -9999.
    y_gradient.grid_mapping = "LatLon_Projection"
    y_gradient.coordinates = "time altitude lat lon "
    # float magnitude(time, altitude, lat, lon) ;
    magnitude_gradient.long_name = "Magnitude of SST Gradient"
    magnitude_gradient.missing_value = -9999.
    magnitude_gradient.grid_mapping = "LatLon_Projection"
    magnitude_gradient.coordinates = "time altitude lat lon "
    ## global
    ncfile.title = "Daily estimated MUR SST Frontal edges, x_gradient, y_gradient and gradient magnitude"
    ncfile.cdm_data_type = "Grid"
    ncfile.Conventions = "COARDS, CF-1.6, ACDD-1.3"
    ncfile.standard_name_vocabulary = "CF Standard Name Table v55"
    ncfile.creator_email = "erd.data@noaa.gov"
    ncfile.creator_name =  "NOAA NMFS SWFSC ERD"
    ncfile.creator_type =  "institution"
    ncfile.creator_url  = "https://www.pfeg.noaa.gov"
    ncfile.Easternmost_Easting = lon_max
    ncfile.Northernmost_Northing = lat_max
    ncfile.Westernmost_Easting = lon_min
    ncfile.Southernmost_Northing =  lat_max
    ncfile.geospatial_lat_max = lat_max
    ncfile.geospatial_lat_min =  lat_min
    ncfile.geospatial_lat_resolution = 0.01
    ncfile.geospatial_lat_units = "degrees_north"
    ncfile.geospatial_lon_max = lon_max
    ncfile.geospatial_lon_min = lon_min
    ncfile.geospatial_lon_resolution = 0.01
    ncfile.geospatial_lon_units = "degrees_east"
    ncfile.infoUrl = ""
    ncfile.institution = "NOAA ERD"
    ncfile.keywords = ""
    ncfile.keywords_vocabulary = "GCMD Science Keywords"
    ncfile.summary = '''Front Edges estimated from daily MUR SST files
    using the Python scikit-image canny algorithm  with sigma = 10., and
    threshold values of .8 and .9,  as well as the OpenCV algorithm findContours with a minimum length of 20.
    The SST x-gradient, y-gradient and gradient magnitude are also included
    '''
    ncfile.license = '''The data may be used and redistributed for free but is not intended
    for legal use, since it may contain inaccuracies. Neither the data
    Contributor, ERD, NOAA, nor the United States Government, nor any
    of their employees or contractors, makes any warranty, express or
    implied, including warranties of merchantability and fitness for a
    particular purpose, or assumes any legal liability for the accuracy,
    completeness, or usefulness, of this information.
    '''
    file_name1 = c_file_year + c_file_month + c_file_day + '090000-JPL-L4_GHRSST-SSTfnd-MUR-GLOB-v02.0-fv04.1.nc'
    history = 'created from MUR SST file ' + file_name1 + 'using python scikit-image canny algorithm, sigma = 10, thresholds of 0.8, 0.9 and OpenCV findContours function with minimum length 20'
    ncfile.history = history
    altitude[0] = 0.
    longitude[:] = lons[:]
    latitude[:] = lats[:]
    ncfile.close()
    return file_name

"""
canny.py - Canny Edge detector

Reference: Canny, J., A Computational Approach To Edge Detection, IEEE Trans.
    Pattern Analysis and Machine Intelligence, 8:679-714, 1986

Originally part of CellProfiler, code licensed under both GPL and BSD licenses.
Website: http://www.cellprofiler.org
Copyright (c) 2003-2009 Massachusetts Institute of Technology
Copyright (c) 2009-2011 Broad Institute
All rights reserved.
Original author: Lee Kamentsky
"""

import numpy as np
import scipy.ndimage as ndi
from scipy.ndimage import generate_binary_structure, binary_erosion, label
from skimage.filters import gaussian
from skimage import dtype_limits, img_as_float
from skimage._shared.utils import assert_nD


def smooth_with_function_and_mask(image, function, mask):
    """Smooth an image with a linear function, ignoring masked pixels

    Parameters
    ----------
    image : array
        Image you want to smooth.
    function : callable
        A function that does image smoothing.
    mask : array
        Mask with 1's for significant pixels, 0's for masked pixels.

    Notes
    ------
    This function calculates the fractional contribution of masked pixels
    by applying the function to the mask (which gets you the fraction of
    the pixel data that's due to significant points). We then mask the image
    and apply the function. The resulting values will be lower by the
    bleed-over fraction, so you can recalibrate by dividing by the function
    on the mask to recover the effect of smoothing from just the significant
    pixels.
    """
    bleed_over = function(mask.astype(float))
    masked_image = np.zeros(image.shape, image.dtype)
    masked_image[mask] = image[mask]
    smoothed_image = function(masked_image)
    output_image = smoothed_image / (bleed_over + np.finfo(float).eps)
    return output_image


def canny1(image, sigma=1., low_threshold=None, high_threshold=None, mask=None,
          use_quantiles=False):
    """Edge filter an image using the Canny algorithm.

    Parameters
    -----------
    image : 2D array
        Grayscale input image to detect edges on; can be of any dtype.
    sigma : float
        Standard deviation of the Gaussian filter.
    low_threshold : float
        Lower bound for hysteresis thresholding (linking edges).
        If None, low_threshold is set to 10% of dtype's max.
    high_threshold : float
        Upper bound for hysteresis thresholding (linking edges).
        If None, high_threshold is set to 20% of dtype's max.
    mask : array, dtype=bool, optional
        Mask to limit the application of Canny to a certain area.
    use_quantiles : bool, optional
        If True then treat low_threshold and high_threshold as quantiles of the
        edge magnitude image, rather than absolute edge magnitude values. If True
        then the thresholds must be in the range [0, 1].

    Returns
    -------
    output : 2D array (image)
        The binary edge map.

    See also
    --------
    skimage.sobel

    Notes
    -----
    The steps of the algorithm are as follows:

    * Smooth the image using a Gaussian with ``sigma`` width.

    * Apply the horizontal and vertical Sobel operators to get the gradients
      within the image. The edge strength is the norm of the gradient.

    * Thin potential edges to 1-pixel wide curves. First, find the normal
      to the edge at each point. This is done by looking at the
      signs and the relative magnitude of the X-Sobel and Y-Sobel
      to sort the points into 4 categories: horizontal, vertical,
      diagonal and antidiagonal. Then look in the normal and reverse
      directions to see if the values in either of those directions are
      greater than the point in question. Use interpolation to get a mix of
      points instead of picking the one that's the closest to the normal.

    * Perform a hysteresis thresholding: first label all points above the
      high threshold as edges. Then recursively label any point above the
      low threshold that is 8-connected to a labeled point as an edge.

    References
    -----------
    .. [1] Canny, J., A Computational Approach To Edge Detection, IEEE Trans.
           Pattern Analysis and Machine Intelligence, 8:679-714, 1986
    .. [2] William Green's Canny tutorial
           http://dasl.unlv.edu/daslDrexel/alumni/bGreen/www.pages.drexel.edu/_weg22/can_tut.html

    Examples
    --------
    >>> from skimage import feature
    >>> # Generate noisy image of a square
    >>> im = np.zeros((256, 256))
    >>> im[64:-64, 64:-64] = 1
    >>> im += 0.2 * np.random.rand(*im.shape)
    >>> # First trial with the Canny filter, with the default smoothing
    >>> edges1 = feature.canny(im)
    >>> # Increase the smoothing for better results
    >>> edges2 = feature.canny(im, sigma=3)
    """

    #
    # The steps involved:
    #
    # * Smooth using the Gaussian with sigma above.
    #
    # * Apply the horizontal and vertical Sobel operators to get the gradients
    #   within the image. The edge strength is the sum of the magnitudes
    #   of the gradients in each direction.
    #
    # * Find the normal to the edge at each point using the arctangent of the
    #   ratio of the Y sobel over the X sobel - pragmatically, we can
    #   look at the signs of X and Y and the relative magnitude of X vs Y
    #   to sort the points into 4 categories: horizontal, vertical,
    #   diagonal and antidiagonal.
    #
    # * Look in the normal and reverse directions to see if the values
    #   in either of those directions are greater than the point in question.
    #   Use interpolation to get a mix of points instead of picking the one
    #   that's the closest to the normal.
    #
    # * Label all points above the high threshold as edges.
    # * Recursively label any point above the low threshold that is 8-connected
    #   to a labeled point as an edge.
    #
    # Regarding masks, any point touching a masked point will have a gradient
    # that is "infected" by the masked point, so it's enough to erode the
    # mask by one and then mask the output. We also mask out the border points
    # because who knows what lies beyond the edge of the image?
    #
    assert_nD(image, 2)
    dtype_max = dtype_limits(image, clip_negative=False)[1]

    if low_threshold is None:
        low_threshold = 0.1 * dtype_max
    else:
        low_threshold = low_threshold / dtype_max

    if high_threshold is None:
        high_threshold = 0.2 * dtype_max
    else:
        high_threshold = high_threshold / dtype_max

    if mask is None:
        mask = np.ones(image.shape, dtype=bool)

    def fsmooth(x):
        return img_as_float(gaussian(x, sigma, mode='constant'))

    smoothed = smooth_with_function_and_mask(image, fsmooth, mask)
    jsobel = ndi.sobel(smoothed, axis=1)
    isobel = ndi.sobel(smoothed, axis=0)
    abs_isobel = np.abs(isobel)
    abs_jsobel = np.abs(jsobel)
    magnitude = np.hypot(isobel, jsobel)

    #
    # Make the eroded mask. Setting the border value to zero will wipe
    # out the image edges for us.
    #
    s = generate_binary_structure(2, 2)
    eroded_mask = binary_erosion(mask, s, border_value=0)
    eroded_mask = eroded_mask & (magnitude > 0)
    #
    #--------- Find local maxima --------------
    #
    # Assign each point to have a normal of 0-45 degrees, 45-90 degrees,
    # 90-135 degrees and 135-180 degrees.
    #
    local_maxima = np.zeros(image.shape, bool)
    #----- 0 to 45 degrees ------
    pts_plus = (isobel >= 0) & (jsobel >= 0) & (abs_isobel >= abs_jsobel)
    pts_minus = (isobel <= 0) & (jsobel <= 0) & (abs_isobel >= abs_jsobel)
    pts = pts_plus | pts_minus
    pts = eroded_mask & pts
    # Get the magnitudes shifted left to make a matrix of the points to the
    # right of pts. Similarly, shift left and down to get the points to the
    # top right of pts.
    c1 = magnitude[1:, :][pts[:-1, :]]
    c2 = magnitude[1:, 1:][pts[:-1, :-1]]
    m = magnitude[pts]
    w = abs_jsobel[pts] / abs_isobel[pts]
    c_plus = c2 * w + c1 * (1 - w) <= m
    c1 = magnitude[:-1, :][pts[1:, :]]
    c2 = magnitude[:-1, :-1][pts[1:, 1:]]
    c_minus = c2 * w + c1 * (1 - w) <= m
    local_maxima[pts] = c_plus & c_minus
    #----- 45 to 90 degrees ------
    # Mix diagonal and vertical
    #
    pts_plus = (isobel >= 0) & (jsobel >= 0) & (abs_isobel <= abs_jsobel)
    pts_minus = (isobel <= 0) & (jsobel <= 0) & (abs_isobel <= abs_jsobel)
    pts = pts_plus | pts_minus
    pts = eroded_mask & pts
    c1 = magnitude[:, 1:][pts[:, :-1]]
    c2 = magnitude[1:, 1:][pts[:-1, :-1]]
    m = magnitude[pts]
    w = abs_isobel[pts] / abs_jsobel[pts]
    c_plus = c2 * w + c1 * (1 - w) <= m
    c1 = magnitude[:, :-1][pts[:, 1:]]
    c2 = magnitude[:-1, :-1][pts[1:, 1:]]
    c_minus = c2 * w + c1 * (1 - w) <= m
    local_maxima[pts] = c_plus & c_minus
    #----- 90 to 135 degrees ------
    # Mix anti-diagonal and vertical
    #
    pts_plus = (isobel <= 0) & (jsobel >= 0) & (abs_isobel <= abs_jsobel)
    pts_minus = (isobel >= 0) & (jsobel <= 0) & (abs_isobel <= abs_jsobel)
    pts = pts_plus | pts_minus
    pts = eroded_mask & pts
    c1a = magnitude[:, 1:][pts[:, :-1]]
    c2a = magnitude[:-1, 1:][pts[1:, :-1]]
    m = magnitude[pts]
    w = abs_isobel[pts] / abs_jsobel[pts]
    c_plus = c2a * w + c1a * (1.0 - w) <= m
    c1 = magnitude[:, :-1][pts[:, 1:]]
    c2 = magnitude[1:, :-1][pts[:-1, 1:]]
    c_minus = c2 * w + c1 * (1.0 - w) <= m
    local_maxima[pts] = c_plus & c_minus
    #----- 135 to 180 degrees ------
    # Mix anti-diagonal and anti-horizontal
    #
    pts_plus = (isobel <= 0) & (jsobel >= 0) & (abs_isobel >= abs_jsobel)
    pts_minus = (isobel >= 0) & (jsobel <= 0) & (abs_isobel >= abs_jsobel)
    pts = pts_plus | pts_minus
    pts = eroded_mask & pts
    c1 = magnitude[:-1, :][pts[1:, :]]
    c2 = magnitude[:-1, 1:][pts[1:, :-1]]
    m = magnitude[pts]
    w = abs_jsobel[pts] / abs_isobel[pts]
    c_plus = c2 * w + c1 * (1 - w) <= m
    c1 = magnitude[1:, :][pts[:-1, :]]
    c2 = magnitude[1:, :-1][pts[:-1, 1:]]
    c_minus = c2 * w + c1 * (1 - w) <= m
    local_maxima[pts] = c_plus & c_minus

    #
    #---- If use_quantiles is set then calculate the thresholds to use
    #
    if use_quantiles:
        if high_threshold > 1.0 or low_threshold > 1.0:
            raise ValueError("Quantile thresholds must not be > 1.0")
        if high_threshold < 0.0 or low_threshold < 0.0:
            raise ValueError("Quantile thresholds must not be < 0.0")

        high_threshold = np.percentile(magnitude, 100.0 * high_threshold)
        low_threshold = np.percentile(magnitude, 100.0 * low_threshold)

    #
    #---- Create two masks at the two thresholds.
    #
    high_mask = local_maxima & (magnitude >= high_threshold)
    low_mask = local_maxima & (magnitude >= low_threshold)

    #
    # Segment the low-mask, then only keep low-segments that have
    # some high_mask component in them
    #
    strel = np.ones((3, 3), bool)
    labels, count = label(low_mask, strel)
    if count == 0:
        return low_mask

    sums = (np.array(ndi.sum(high_mask, labels,
                             np.arange(count, dtype=np.int32) + 1),
                     copy=False, ndmin=1))
    good_label = np.zeros((count + 1,), bool)
    good_label[1:] = sums > 0
    output_mask = good_label[labels]
    return output_mask, isobel, jsobel, magnitude

