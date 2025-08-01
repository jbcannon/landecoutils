#' Function to clip large LAS object into small pieces
#'
#' This function clips a large LAS object such as that from a
#' UAV-derived point cloud into smaller chunks using the
#' `catalog_map` function.
#' @param las_path character - path to large LAS file
#' @param output_directory character - path to directory for output files.
#' Outputs will be automatically named XLEFT_YBOTTOM.las'
#' @param chunk_size numeric size of the chunks to break the file into
#' @param chunk_alignment coordinates of origin of grid to which files
#' @import lidR
#' @examples
#' library(lidR)
#' library(landecoutils)
#' large_file = 'C:/user/me/large_file.las'
#' output_location = 'C:/user/me/outputs/'
#' tile_large_LAS(large_file, output_location)
#' @export
tile_large_LAS = function(las_path, output_directory, chunk_size=100,
                          chunk_alighment = c(0,0)) {
  ctg = lidR::readLAScatalog(las_path)
  lidR::opt_chunk_size(ctg) = chunk_size
  lidR::opt_chunk_buffer(ctg) = 0
  lidR::opt_chunk_alignment(ctg) = chunk_alignment
  lidR::plot(ctg, chunk_pattern=TRUE)
  lidR::opt_output_files(ctg) = paste0(output_directory, '/{XLEFT}_{YBOTTOM}')
  lidR::catalog_map(ctg, function(las) return(las))
  cat('files output to:\n', output_directory)
  return(NULL)
}

#' Create elevation products from a LAS tile
#'
#' This function takes a LAS object and returns a digital elevation model
#' (dem), a canopy surface model (csm), and a canopy height model (chm)
#' using default algorithms in lidR. Can be used with lidR::catalog_map.
#' Returns a 3-lasyered SpatRast object with names dem, csm, chm.
#' @param las LAS object from lidR package to create model
#' @param res output resolution in m
#' @examples
#' library(terra)
#' library(lidR)
#' las = readLAS('E:/mylas.laz')
#' elev = get_cdem_csm_chm(las)
#' plot(elev$dem)
#' plot(elev$csm)
#' plot(elev$chm)
#'
#' ctg = readLASCatalog('E:/my/las/dir/')
#' elev = catalog_map(ctg, get_dem_csm_chm)
#' plot(elev$dem)
#' plot(elev$csm)
#' plot(elev$chm)
#'
#' @export
get_dem_csm_chm = function(las, res=0.5) {
  las = lidR::classify_ground(las, lidR::csf(class_threshold = 0.5))
  dem = lidR::rasterize_terrain(las, res=res)
  csm = lidR::rasterize_canopy(las, res=res)
  las = lidR::normalize_height(las, lidR::tin())
  chm = lidR::rasterize_canopy(las, res = res)
  out = terra::rast(c(dem=dem,csm=csm,chm=chm))
  return(out)
}

#' Compress a folder of LAS files to LAZ format.
#'
#' This function takes a directory of LAS files and compresses them to LAZ.
#' The original LAS files are deleted assuming you want to save space.
#' @param las_dir path to a directory containing .LAS files to compress
#' @param n_cores number of cores to create doSNOW cluster
#' @examples
#' ## NOT RUN ##
#' # compress_las('E:/my/las/dir/', n_cores=2)
#' @export
compress_las = function(las_dir, n_cores, index=TRUE, delete_old = FALSE) {
  files = list.files(las_dir, '.las', full.names=TRUE)
  cl = parallel::makeCluster(n_cores)
  doSNOW::registerDoSNOW(cl)
  pb = utils::txtProgressBar(max = length(files), style = 3)
  progress = function(n) utils::setTxtProgressBar(pb, n)
  opts = list(progress = progress)
  `%dopar%` = foreach::`%dopar%`
  foreach::foreach(fn=files, .options.snow=opts) %dopar% {
    new_laz_fn = gsub('.las', '.laz', fn)
    if(file.exists(new_laz_fn)) {return(NULL)}
    lidR::writeLAS(lidR::readLAS(fn), new_laz_fn, index=index)
    if(delete_old) {if(file.exists(new_laz_fn)) unlink(fn)}
    return(NULL)
  }
  close(pb)
  parallel::stopCluster(cl)
}

#' Check for and create a lax index from a directory of LAS files
#'
#' Check for LAX index files and create them as necessary.  Indexing greatly
#' increases the processingAn earlier version of this function overwrote the
#' original LAS resulting in possible corruption of the original file if the
#' process was interrupted. This updated version uses temporary files to avoid
#' this problem.
#' @param input can be a LASCatalog from `lidR` package or a path to a directory containing .LAS or .LAZ files to index
#' @param n_cores number of cores to create doSNOW cluster
#' @param write_lax indicates if .lax file should be written (`TRUE`), or only
#' checked for (`FALSE`)
#' @examples
#' check_for_lax('E:/my/las/dir/', n_cores=4)
#' @export
check_for_lax = function(input, n_cores=1, write_lax=TRUE) {
  #check inputs for validity
  if(!write_lax %in% c(TRUE, FALSE)) stop('write_lax must be TRUE/FALSE')
  if(is.character(input))  laz = list.files(dir, pattern='.las$|laz$', full.names = TRUE)
  if(class(input) == 'LASCatalog') laz = input$filename
  if(length(laz) < 0) stop('no .LAS or .LAZ found in `input`')

  #check which files are missing indexes
  lax = list.files(dir, pattern='.lax', full.names = TRUE)
  needs_lax = !gsub('.las$|.laz$', '', laz) %in% gsub('.lax', '', lax)
  cat(sum(needs_lax), 'files need indexing\n')
  if(sum(needs_lax)==0) return(NULL) #exit if none needed.

  #return list of needed indexes if write_lax == false
  if(!write_lax) return(laz[needs_lax])

  # if write_lax == TRUE, add index
  cl = parallel::makeCluster(n_cores)
  doSNOW::registerDoSNOW(cl)
  pb = utils::txtProgressBar(max = length(laz[needs_lax]), style = 3)
  progress = function(n) utils::setTxtProgressBar(pb, n)
  opts = list(progress = progress)
  `%dopar%` = foreach::`%dopar%`
  foreach::foreach(i = laz[needs_lax], .options.snow=opts) %dopar% {
    cat('file', i, '\n...indexing\n')
    rlas::writelax(i)
    return(NULL)
  }
  close(pb)
  parallel::stopCluster(cl)
  cat(length(laz[needs_lax]), ' indexes complete')
  return(NULL)
}

#' Find scan locations from a `LAScatalog`
#'
#' This function does the same as `las_find_centroids` but it works on a
#' `LAScatalog` and returns an `sf` object of point lcoations of scan centroids.
#' Result of this function can be used as an input to `stitch_TLS_dir_to_LAS` and
#' `stitch_TLS_dir_to_LAS_tile` to increase speed.
#'
#' @param ctg a `lidR::LAScatalog` object for which you want to find a centroids of all scan locations
#' @param n_cores a number of cores to use for parallel processing.
#' @param res Resolution to thin the point cloud for efficiency.  See `las_find_centroids`
#' @examples
#' library(sf)
#' ctg = lidR::readLAScatalog('path/to/las/catalog')
#' scan_locations = find_ctg_centroids(ctg, n_cores=4)
#' plot(scan_locations$geom)
#' @export
#'
#'

find_ctg_centroids = function(ctg, n_cores=4, res=1) {
  cl = parallel::makeCluster(n_cores)
  doSNOW::registerDoSNOW(cl)
  pb = utils::txtProgressBar(max = length(ctg$filename), style = 3)
  progress = function(n) utils::setTxtProgressBar(pb, n)
  opts = list(progress = progress)
  `%dopar%` = foreach::`%dopar%`
  centroids = foreach::foreach(fn = ctg$filename, .options.snow=opts) %dopar% {
    gc()
    return(suppressWarnings(landecoutils::find_las_centroid(fn, res=res)))
  }
  close(pb)
  parallel::stopCluster(cl)
  centroids = do.call(rbind, centroids)
  gc()
  return(centroids)
}

#' Get Coordinate Reference System (CRS) from LASCatalog
#'
#' This function takes a LASCatalog as input and returns a `crs` object
#' representing the coordinate reference system of the LASCatalog
#' @param ctg a `LAScatalog` object containing LAS files
#' @examples
#' # Load LAScatalog and capture the crs
#' ctg = lidR::readLAScatalog('path/to/LASfiles/')
#' myCRS = get_crg_crs(ctg)
#' print(myCRS)
#' @export
get_ctg_crs = function(ctg){
    hdr = lidR::readLASheader(ctg$filename[1])
    proj = sf::st_crs(hdr@VLR$`WKT OGC CS`$`WKT OGC COORDINATE SYSTEM`)
    return(proj)
}

#' Stitch TLS scans into single *.las
#'
#' This function takes an input directory of TLS scans that are overlapping
#' and a region of interest. The function loads, stitches, and clips
#' all las data to the region of interest and outputs into a single file.
#' Note that this should only be used on relatively small areas to avoid
#' very large file sizes and memory errors. Plot sizes < 0.2 ha are recommended
#' including buffer unless points are sparse. *Note*: script will
#' keep only data columns which are common among all scans.
#' @param ctg a `LAScatalog` object containing LAS files to be stitched
#' @param out_las path for output *.las file
#' @param roi an `sf` object denoting the region of interest. UTM units preferred
#' @param buffer width of buffer to apply to `roi` using `sf::st_buffer`
#' @param max_scan_distance maximum distance from scan to be included. Evaluated using `find_las_centroid()`
#' @param index boolean. Also write a lax file to index the points in the files. see `lidR::writeLAS`
#' @param scan_locations (optional) an `sf` object of scan centroids index the same as `ctg`.
#' If `scan_locations` is `NULL` then scan locations will be found automatically
#' @examples
#' # Load LAScatalog and clip to a region of interest specified by an sf object
#' ctg = lidR::readLAScatalog('path/to/LASfiles/')
#' roi = sf::st_read('plot_boundary.shp')
#' stitch_TLS_dir_to_LAS(ctg, 'output_path.las', roi)
#' @export
#'

stitch_TLS_dir_to_LAS = function(ctg, out_las, roi, buffer = 10, max_scan_distance=60, index=TRUE, scan_locations=NULL) {

  # Load plot boundaries/buffer and create filter
  hdr = lidR::readLASheader(ctg$filename[1])
  proj = sf::st_crs(hdr@VLR$`WKT OGC CS`$`WKT OGC COORDINATE SYSTEM`)
  suppressMessages(sf::st_crs(roi) <- proj)


  roi_buff = sf::st_buffer(roi, dist=buffer)
  ex = sf::st_bbox(roi_buff)
  #tmp = ex/tile_size
  #ex = c(floor(tmp[1]), floor(tmp[2]), ceiling(tmp[3]), ceiling(tmp[4]))*tile_size
  filt = paste('-keep_xy', ex[1], ex[2], ex[3], ex[4]) #min_x min_y max_x max_y
  lidR::opt_filter(ctg) = filt

  #Display catalog and ROI
  lidR::plot(ctg)
  plot(roi[attr(roi, "sf_column")], add=TRUE, lwd=2, col='white', border='white')
  plot(roi_buff[attr(roi_buff, "sf_column")], add=TRUE, lty=2, border='white')
  Sys.sleep(0.5)

  # Load TLS scans from directory and clip to roi.. rbind, and write to file.
  i = 1
  combined_las = list()
  if(is.null(scan_locations)) scan_locations = suppressWarnings(find_ctg_centroids(ctg))
  for(fn in ctg$filename) {
    cat('reading las', i, 'of', length(ctg$filename), '\n')
    scan_location = sf::st_buffer(scan_locations[i,], dist=max_scan_distance)
    intx = sf::st_intersection(roi_buff, scan_location)
    if(nrow(intx)==0) {cat('...out of bounds\n'); i=i+1; next }
    combined_las[[i]] = lidR::clip_roi(lidR::readLAS(fn, filter=filt), intx)
    i=i+1
  }
  #Get rid of any empty items in the list.
  combined_las = combined_las[unlist(lapply(combined_las, function(x) class(x)!='NULL'))]

  #To avoid rbind errors, Find common columns among scans and keep only those
  common_cols = lapply(combined_las, function(x) colnames(x@data))
  common_cols = Reduce(intersect, common_cols)
  cols_to_keep = c('X', 'Y', 'Z', 'gpstime', 'Amplitude', 'Intensity', 'ReturnNumber', "NumberOfReturns", 'Classification', 'Reflectance', 'Deviation')
  if('Distance' %in% common_cols) cols_to_keep = c(cols_to_keep, 'Distance')
  combined_las = lapply(combined_las, function(x) {
    x@data = x@data[, ..cols_to_keep]
    return(x)})
      #make sure las portion with highest NumberofReturns is listed first so bit count is set correctly
  whichMaxReturns = which.max(sapply(combined_las, function(x) max((x@data$NumberOfReturns))))
  n = c(whichMaxReturns, (1:length(combined_las))[-whichMaxReturns])
  combined_las = do.call(rbind,combined_las[n])
  combined_las@header@VLR = list()
  lidR::st_crs(combined_las) = proj
  lidR::writeLAS(lidR::las_update(combined_las), out_las, index=index)
  cat('combined las written to', out_las, '\n')
  return(combined_las)
}

#' Find centroid of large *.las
#'
#' This function takes a `lidR::LAS` object, and reads it in quickly by dropping
#' most points using `subsample`. The function returns an `sf::POINT` object
#' containing the coordinates of the centroid. *Note* that centroid is only
#' approximate and may be imprecise if subsample is large (i.e., > 1000) or
#' LAS point density is low.
#'
#' @param las a `lidR::LAS` object for which you want to find a centroid
#' @param res resolution to thin the point cloud for efficiency
#' @examples
#' # Load large LAS file and identify the centroid (i.e., scan location)
#' library(lidR)
#' library(sf)
#' #' las = readLAS('large_las_file.las')
#' cent = find_las_centroid(las)
#' plot(cent$geom)
#' @export
find_las_centroid = function(las, res = 1) {
  filt = paste0('-thin_with_grid ', res)
  las_thin = lidR::readLAS(las, filt)
  centroid = apply(las_thin@data[, c('X', 'Y')], 2, mean)
  centroid = as.data.frame(t(centroid))
  centroid$id = 'centroid'
  centroid$fn = las
  centroid = sf::st_as_sf(centroid, coords = c('X', 'Y'))
  sf::st_crs(centroid) = sf::st_crs(las_thin)
  return(centroid)
}

#' Add attribute `Distance` to TLS returns based on distance from scanner
#'
#' This function takes the path of a .las or .laz file and return a `LAS` object
#' containing the new attribute `Distance` represetning the distance to the scanner
#' The scanner location is determined by (1) using the `find_las_centroid()` function,
#' then (2) locates the scanner elevation by adding the distance from the `scanner_ht`
#' parameter to the elevation of the centroid. Elevation is determined using the
#' `lidR::classify_ground` and `lidR::rasterize_terrain` functions.
#' @param las_filename path to a `.las` or `.laz` file
#' @param scanner_ht height of scanner above ground. Defaults to 1.75
#' @param subsample subsamplign factor used in `find_las_centroid()` to reduce lidar resoultion. Defaults to 1e05.
#' @examples
#'
#' las_fn = 'data.laz'
#' las = las_add_scanner_distance(las_fn)
#' lidR::writeLAS(las, 'data_withDistance.laz')
#' @export
las_add_scanner_distance = function(las_filename,
                                    scanner_ht = 1.75,
                                    subsample = 1e5) {
  # Does it already have a Distance column?
  hdr = lidR::readLASheader(las_filename)
  extrabytes = names(hdr@VLR$Extra_Bytes$`Extra Bytes Description`)
  if('Distance' %in% extrabytes) {
    warning('source file already contains Distance column. Returning NULL')
    return(NULL)
  }

  # Find XYZ location of scan from centroid and elevation
  message('Step 1: Finding scanner XY location')
  check_for_lax(las_filename)
  centroid = find_las_centroid(las_filename, subsample=subsample)

  message('Step 2: Mapping elevation to capture scanner Z location')
  filt = paste0('-keep_circle ', paste0(sf::st_coordinates(centroid), collapse = ' '), ' 5')
  las = lidR::readLAS(las_filename, filter=filt)
  las = lidR::classify_ground(las, lidR::csf())
  dem = lidR::rasterize_terrain(las, 1, lidR::tin())
  scanner_elev = terra::extract(dem, centroid)$Z
  scanner_elev = scanner_elev + scanner_ht
  scanner_loc = as.data.frame(sf::st_coordinates(centroid))
  scanner_loc$Z = scanner_elev
  message(paste0('\tScanner location: ', paste0(round(as.numeric(scanner_loc),1), collapse=' ')))

  message('Step 3: Loading full resolution las')
  las = lidR::readLAS(las_filename)

  # Use 3d distance function from lidar returns to scanner location
  message('Step 4: Calculating distance between returns and scanner')
  mat = dplyr::as_tibble((dplyr::select(las@data, X, Y, Z)))
  loc = as.matrix(scanner_loc)
  dx = (mat[,1] - loc[1])^2
  dy = (mat[,2] - loc[2])^2
  dz = (mat[,3] - loc[3])^2
  dist = sqrt(dx + dy + dz)
  dist = dist$X

  # add attribute to las and return
  message('Step 5: Attaching distance measurements to las')
  las = lidR::add_lasattribute(las, dist, 'Distance', 'Distance in meters from scanner')
  return(las)
}

#' Correct Return Reflectance based on Distance using Xu et al. 2017
#'
#' This function takes a `LAS` input containing the fields `Amplitude` in dB and
#' `Distance` in meteres. It applies a distance based correction using the methods
#' of Xu et al. 2017. The formula only applies the correction based on distance, but
#' not incident angle. The study used a REIGL VZ400i, but less sure if this applies
#' to other scanners.
#' @param las a `LAS` object containing attributes for `Amplitude` and `Distance`
#' @examples
#' ## NOT RUN ##
#' #las = readLAS('las1.laz')
#' # corr = correct_Reflectance_Xu2017(las)
#' # par(mfrow = c(1,2))
#' # plot(las$Amplitude ~ las$Distance)
#' # plot(corr, las$Distance)
#' @export
correct_Reflectance_Xu2017 = function(las) {
  stopifnot(class(las) == 'LAS')
  stopifnot('Amplitude' %in% colnames(las@data))
  stopifnot('Distance' %in% colnames(las@data))
  R = las$Distance
  A = las$Amplitude
  # Piecewise distance-base correction alogrith (F1) from Xu et al. 2017 Rem. Sens.
  F11 = 1.623e-3*R^3 - 9.287e-2*R^2 +1.37*R + 25.88
  F12 = 10*log(3.218e5/R^2, base=10)
  Ic = ifelse(R < 20, F11, F12)
  return(A - Ic)
}

#' Combine overlapping TLS scans into tiled LAS scene
#'
#' This function takes an input directory of TLS scans that are overlapping
#' and a region of interest (e.g., large plot). The function loads, stitches, and clips
#' all las data to the region of interest and outputs LAS tiles of a user specified
#' size. *Note*: when binding overlapping scan sections, script will
#' keep only data columns which are common among all scans.
#' @param ctg a `LAScatalog` object containing overlapping TLS scans
#' @param out_dir directory to output LAS tiles (no ending '/' in path).
#' @param bnd an `sf` object denoting the region of interest. Only UTM units tested
#' @param buffer width of buffer to apply to `bnd` using `sf::st_buffer` This
#' allows inclusion of TLS data outside study area to address edge effects.
#' @param n_cores number of cores to use in parallel processing
#' @param max_scan_distance maximum distance from scan to be included. Evaluated using `find_las_centroid()`
#' @param index boolean. Also write a lax file to index the points in the files. see `lidR::writeLAS`
#' @param timeout numeric in seconds. Number of seconds for timeout on a tile
#' before proceeding to the next one. This avoids cores getting hung up on
#' problem tiles, so that that can be trouble shooted separately. This may not
#' be necessary in later versions.
#' @examples
#' # Load LAScatalog, clip, and tile for a large area specified by an sf object
#' ctg = lidR::readLAScatalog('path/to/LASfiles/')
#' bnd = sf::st_read('plot_boundary.shp')
#' stitch_TLS_dir_to_LAS_tiles(ctg, 'output_tiles', bnd, tile_size = 30)
#' @export
stitch_TLS_dir_to_LAS_tiles = function(ctg, out_dir, tile_size, n_cores,
                                       buffer = 10, max_scan_distance=60,
                                       bnd=NULL, index=TRUE, scan_locations=NULL,
                                       timeout = 600) {
  # Check for duplicate scans
  check_for_duplicate_scans(ctg)
  if(!dir.exists(out_dir)) {
    cat(out_dir, 'does not exist, creating')
    dir.create(out_dir, recursive = TRUE)
  }

  # Clean up or create Plot boundaries and place in correct projection.
  hdr = lidR::readLASheader(ctg$filename[1])
  proj = sf::st_crs(hdr@VLR$`WKT OGC CS`$`WKT OGC COORDINATE SYSTEM`)
  # Generate scan locations if they do not exist
  if(is.null(scan_locations)) {
    tmp = landecoutils::find_ctg_centroids(ctg, n_cores = 6, res=1)
    scan_locations = sf::st_tran
    sform(tmp, proj)
  }
  # Generate boundary from scan_locations if it does not exist
  if(is.null(bnd)) {
    bnd = scan_locations
    bnd = sf::st_convex_hull(sf::st_union(bnd))
    sf::st_crs(bnd) = sf::st_crs(scan_locations)
  }

  bnd = sf::st_transform(bnd, proj)
  scan_locations = sf::st_transform(scan_locations, proj)

  # Create ctg filter based on boundary
  bnd_buff = sf::st_buffer(bnd, dist=buffer)
  ex = sf::st_bbox(bnd_buff)
  filt = paste('-keep_xy', ex[1], ex[2], ex[3], ex[4]) #min_x min_y max_x max_y
  tmp = ex/tile_size
  ex = c(floor(tmp[1]), floor(tmp[2]), ceiling(tmp[3]), ceiling(tmp[4]))*tile_size
  lidR::opt_filter(ctg) = filt

  #create fishnet from extent
  ncol = (ex[3] - ex[1])/tile_size
  nrow = (ex[4] - ex[2]) / tile_size
  grid = terra::rast(extent=terra::ext(ex[c(1,3,2,4)]), ncols=ncol, nrows=nrow)
  grid = terra::as.polygons(grid)
  grid = sf::st_as_sf(grid)
  sf::st_crs(grid) = proj
  grid = grid[sf::st_intersects(grid, bnd_buff, sparse=FALSE),]

  #Check to see what all is already complete
  todo_list = c()
  for(t in 1:nrow(grid)){
    tile = grid[t,]
    ex = round(sf::st_bbox(tile),1)
    out_las = paste0(out_dir, '/', ex[1], '_', ex[2], '.laz')
    out_las = gsub('\\/\\/', '\\/', out_las)
    if(!file.exists(out_las))  todo_list = c(todo_list, t)
  }
  if(length(todo_list)<1) {
    cat('all scans already tiled')
    return(NULL)
  }
  todo_list = sample(todo_list) #reorder to avoid pileups

  scan_coverage = sf::st_buffer(scan_locations, dist = max_scan_distance)
  plot(scan_coverage$geometry, col = grDevices::rgb(0, 0, 1, 0.05))
  plot(grid$geometry, add = TRUE, border = "red")
  plot(bnd, lwd = 2, border = "white", add = TRUE)
  Sys.sleep(0.5)

  # run through grid tiles, load proximal TLS scans from directory and clip to bnd. rbind, and write to file.
  cl = parallel::makeCluster(n_cores)
  doSNOW::registerDoSNOW(cl)
  pb = utils::txtProgressBar(max = length(todo_list), style = 3)
  progress = function(n) utils::setTxtProgressBar(pb, n)
  opts = list(progress = progress)
  `%dopar%` = foreach::`%dopar%`
  if(is.null(timeout)) timeout = 1e9
  out = foreach::foreach(t=todo_list, .errorhandling = 'stop', .packages=c('sf', 'lidR'), .options.snow=opts) %dopar% {
    # Load and display tile
    R.utils::withTimeout({
      cat('\nloading tile', t, 'of', nrow(grid), '\n')
      tile = grid[t,]
      ex = sf::st_bbox(tile)
      out_las = paste0(out_dir, '/', ex[1], '_', ex[2], '.laz')
      out_las = gsub('\\/\\/', '\\/', out_las)
      if(file.exists(out_las)) return(NULL)
      scans_to_load = which(sf::st_intersects(tile, scan_coverage, sparse = FALSE))
      scans_to_load = dplyr::slice(scan_locations, scans_to_load)
      if(nrow(scans_to_load)<1) return(NULL)
      extent_filter = paste('-keep_xy', ex[1], ex[2], ex[3], ex[4]) #min_x min_y max_x max_y

      #loop through relevant scans, clip and combine
      combined_las = list()
      for(i in 1:nrow(scans_to_load)) {
        cat('.....appending scan', i, 'of', nrow(scans_to_load), '\n')
        scan_coords = sf::st_coordinates(scans_to_load)[i,c('X','Y')]
        distance_filter = paste('-keep_circle ', round(scan_coords['X'],2), round(scan_coords['Y'],2), max_scan_distance)
        combined_filter = paste(extent_filter, distance_filter)
        combined_las[[i]] = lidR::readLAS(scans_to_load[i,]$fn, filter=combined_filter)
      }

      #To avoid rbind errors, Find common columns among scans and keep only those
      combined_las = combined_las[!sapply(combined_las, is.null)]
      combined_las = combined_las[sapply(combined_las, nrow)>0]
      common_cols = lapply(combined_las, function(x) colnames(x@data))
      common_cols = Reduce(intersect, common_cols)
      cols_to_keep = c('X', 'Y', 'Z', 'gpstime', 'Amplitude', 'Intensity', 'ReturnNumber', "NumberOfReturns", 'Classification', 'Reflectance', 'Deviation')
      if('Distance' %in% common_cols) cols_to_keep = c(cols_to_keep, 'Distance')
      combined_las = lapply(combined_las, function(x) {
        x@data = dplyr::select(x@data, dplyr::all_of(cols_to_keep))
        return(x@data)})
      #make sure las portion with highest NumberofReturns is listed first so bit count is set correctly
      whichMaxReturns = which.max(sapply(combined_las, function(x) max((x$NumberOfReturns))))
      n = c(whichMaxReturns, (1:length(combined_las))[-whichMaxReturns])
      combined_las = do.call(rbind,combined_las[n])
      combined_las = dplyr::filter(combined_las, NumberOfReturns < 8)
      combined_las = lidR::LAS(combined_las, crs = proj, check = TRUE)
      combined_las = lidR::las_update(combined_las)


      #write tile to disk
      cat('.....scans stitched. writing tile to disk')
      lidR::writeLAS(combined_las, out_las, index=TRUE)
      return(NULL)
      }, timeout = timeout)
  }
  close(pb)
  parallel::stopCluster(cl)
  return(out)
}

#' Clip a LAS catalog to a boundary to speed up processing
#'
#' This function takes an input directory of TLS scans that are potentially
#' overlapping and a region of interest (`bnd`). The function loads and clips
#' all las data to `bnd` and outputs LAS tiles with an extent
#' matching `bnd`.
#' @param dir a directory of LAS/LAZ files
#' @param output_dir directory to output LAS tiles
#' @param bnd an `sf` object denoting the region of interest. Only UTM units tested
#' @export
clip_las_catalog = function(dir, bnd, output_dir) {
  #check inputs
  needs_index = length(landecoutils::check_for_lax(dir, write_lax = FALSE)) > 0
  if(needs_index) stop('LAS catalog needs indexing. See function `check_for_lax`')
  if(!'sf' %in% class(bnd)) stop('`bnd` must be an `sf` object')
  laz = list.files(dir, pattern='.las|laz', full.names = TRUE)
  if(length(laz) < 0) stop('no .LAS or .LAZ found in `dir`')
  ctg = lidR::readLAScatalog(dir)
  if(!lidR::st_crs(bnd) == lidR::st_crs(ctg)) {
    warning('crs of `bnd` and `ctg` do not match. Reprojecting')
    bnd = sf::st_transform(bnd, sf::st_crs(ctg))
  }
  # make filter
  filt = paste('-keep_xy', sf::st_bbox(bnd) %>% paste(collapse=' '))
  cat('clipping', length(laz),'files\n')
  #clip laz
  for(i in laz) {
    fn = paste0(output_dir, '/', basename(i))
    if(file.exists(fn)) next
    cat('file', basename(i),'...\n...reading\n')
    x = lidR::readLAS(i, filter = filt)
    cat('...writing\n')
    lidR::writeLAS(x, fn, index=TRUE)
  }
}

#' Check LAS Catalog for duplicate scans, this will reduce errors elsewhere.
#'
#' This function takes a LAS Catalog and throws an error
#' if any duplicate scan positions are detected. These files are sometimes
#' generated if a bad scan was taken in the field and a second scan
#' was generated from the same position
#' @param ctg a LAS Catalog
#' @export
check_for_duplicate_scans = function(ctg) {
  fns = ctg$filename
  scanpos = gsub('ScanPos', '', substr(basename(fns),0,10))
  dups = duplicated(scanpos)
  dups = which(scanpos == scanpos[dups])
  dups = basename(fns)[dups]
  if(length(dups) > 0) {
    cat(paste0('possible duplicate scans detected:\n\t', paste0(dups, collapse='\n\t'), sep=''))
    stop('dup scans detected')
  } else {
    cat('no duplicate scans detected\n')
  }
  return(NULL)
}
