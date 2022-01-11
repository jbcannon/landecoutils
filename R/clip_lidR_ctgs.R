#' Find scan locations from a `LAScatalog`
#' 
#' This function does the same as `las_find_centroids` but it works on a
#' `LAScatalog` and returns an `sf` object of point lcoations of scan centroids. 
#' Result of this function can be used as an input to `stitch_TLS_dir_to_LAS` and
#' `stitch_TLS_dir_to_LAS_tile` to increase speed.
#' 
#' @param ctg a `lidR::LAScatalog` object for which you want to find a centroids of all scan locations
#' @param subsample drop every `nth` return to speed up processing time. See `las_find_centroids`
#' @examples 
#' library(sf)
#' ctg = readTLScatalog('path/to/las/catalog')
#' scan_locations = find_ctg_centroids(ctg)
#' plot(scan_locations$geom)
#' @export
find_ctg_centroids = function(ctg, subsample=1e5) {
  if(subsample>1000) warning('Only ', round(1/subsample*100,3), '% of cloud used, centroid may be imprecise')
  s = as.character(format(subsample, scientific=FALSE))
  filt = paste0('-keep_every_nth ', s)
  centroids = list()
  i=1
  for(fn in ctg$filename) {
    cat('finding centroid for scan', i, 'of', length(ctg$filename), '\n\t', basename(fn),'\n')
    centroids[[fn]] = suppressWarnings(find_las_centroid(fn, subsample=subsample))
    i=i+1
  }
  centroids = do.call(rbind, centroids)
  return(centroids)
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
#' ctg = readLAScatalog('path/to/LASfiles/')
#' roi = sf::st_read('plot_boundary.shp')
#' stitch_TLS_dir_to_LAS(ctg, 'output_path.las', roi)
#' @export
stitch_TLS_dir_to_LAS = function(ctg, out_las, roi, buffer = 10, max_scan_distance=60, index=TRUE, scan_locations=NULL) {
  #require(lidR)  
  #require(sf)
  
  # Load plot boundaries/buffer and create filter
  hdr = lidR::readLASheader(ctg$filename[1])
  proj = sf::st_crs(hdr@VLR$`WKT OGC CS`$`WKT OGC COORDINATE SYSTEM`)
  suppressMessages(sf::st_crs(roi) <- proj)
  roi_buff = sf::st_buffer(roi, dist=buffer)
  ex = sf::st_bbox(roi_buff)
  filt = paste('-keep_xy', ex[1], ex[2], ex[3], ex[4]) #min_x min_y max_x max_y
  lidR::opt_filter(ctg) = filt
  
  #Display catalog and ROI
  lidR::plot(ctg)
  plot(roi$geometry, add=TRUE, lwd=2, col='white', border='white')
  plot(roi_buff$geometry, add=TRUE, lty=2, border='white')
  Sys.sleep(0.5)
  
  # Load TLS scans from directory and clip to roi.. rbind, and write to file.
  i = 1
  combined_las = list()
  if(is.null(scan_locations)) scan_locations = suppessWarnings(find_ctg_centroids(ctg))
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
  combined_las = lapply(combined_las, function(x) {
    x@data = x@data[, c('X', 'Y', 'Z', 'gpstime', 'Intensity', 'ReturnNumber', "NumberOfReturns", 'Classification', 'Reflectance', 'Deviation')]
    return(x)})
  combined_las = do.call(rbind,combined_las)
  combined_las@header@VLR = list()
  lidR::writeLAS(las_update(combined_las), out_las, index=index)
  cat('combined las written to', out_las, '\n')
  return(NULL)
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
#' @param subsample drop every `nth` return to speed up processing time
#' @examples 
#' # Load large LAS file and identify the centroid (i.e., scan location)
#' library(lidR)
#' library(sf)
#' #' las = readLAS('large_las_file.las')
#' cent = find_las_centroid(las)
#' plot(cent$geom)
#' @export
find_las_centroid = function(las, subsample=1e5) {
  s = as.character(format(subsample, scientific=FALSE))
  filt = paste0('-keep_every_nth ', s)
  las_thin = lidR::readLAS(las, filter=filt)
  centroid = apply(las_thin@data[, c('X', 'Y', 'Z')], 2, mean)
  centroid = sf::st_point(centroid[1:2])
  out = data.frame(id='centroid')
  out$geom = sf::st_sfc(centroid)
  centroid = sf::st_as_sf(out)
  centroid$fn = las
  sf::st_crs(centroid) = sf::st_crs(las_thin)
  if(subsample>1000) warning('Only ', round(1/subsample*100,3), '% of cloud used, centroid may be imprecise')
  return(centroid)
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
#' @param max_scan_distance maximum distance from scan to be included. Evaluated using `find_las_centroid()`
#' @param index boolean. Also write a lax file to index the points in the files. see `lidR::writeLAS`
#' @examples 
#' # Load LAScatalog, clip, and tile for a large area specified by an sf object
#' ctg = readLAScatalog('path/to/LASfiles/')
#' bnd = sf::st_read('plot_boundary.shp')
#' stitch_TLS_dir_to_LAS_tiles(ctg, 'output_tiles', bnd, tile_size = 30)
#' @export
stitch_TLS_dir_to_LAS_tiles = function(ctg, out_dir, bnd, tile_size, buffer = 10, max_scan_distance=60, index=TRUE, scan_locations=NULL) {
  #require(lidR)  
  #require(sf)
  
  # Load plot boundaries/buffer and create filter
  hdr = lidR::readLASheader(ctg$filename[1])
  proj = sf::st_crs(hdr@VLR$`WKT OGC CS`$`WKT OGC COORDINATE SYSTEM`)
  suppressMessages(sf::st_crs(bnd) <- proj)
  bnd_buff = sf::st_buffer(bnd, dist=buffer)
  ex = sf::st_bbox(bnd_buff)
  filt = paste('-keep_xy', ex[1], ex[2], ex[3], ex[4]) #min_x min_y max_x max_y
  lidR::opt_filter(ctg) = filt
  
  #create fishnet from extent
  ex = round(terra::ext(bnd_buff), 1)
  ncol = ceiling(diff(ex[1:2]) / tile_size)
  nrow = ceiling(diff(ex[3:4]) / tile_size)
  ex[2] = ex[1]  + ncol * tile_size
  ex[4] = ex[3]  + nrow * tile_size
  grid = terra::rast(extent=ex, ncols=ncol, nrows=nrow)
  grid = terra::as.polygons(grid)
  grid = sf::st_as_sf(grid)
  sf::st_crs(grid) = proj
  ex = sf::st_as_sf(terra::as.polygons(round(terra::ext(bnd_buff), 1)))
  grid = grid[sf::st_intersects(grid, bnd_buff, sparse=FALSE),]
  
  #Display catalog and grid
  lidR::plot(ctg)
  plot(bnd$geometry, add=TRUE, lwd=2, border='white')
  plot(grid$geometry, add=TRUE, border='white')
  Sys.sleep(0.5)
  
  #Get all scan centroids once
  if(is.null(scan_locations)) {
    cat('finding all scan footprints within', max_scan_distance, 'meters\n')
    scan_locations = list()
    i=1
    for(fn in ctg$filename) {
      cat('.....',i, 'of', length(ctg), '\n')
      scan_centroid = suppressMessages(find_las_centroid(fn))
      scan_location = sf::st_buffer(scan_centroid, dist=max_scan_distance)
      scan_location$fn=fn
      scan_locations[[length(scan_locations)+1]] = scan_location
      i=i+1
    }; scan_locations = do.call(rbind , scan_locations)
  }
  
  scan_locations = sf::st_buffer(scan_locations, dist=max_scan_distance)
  plot(grid$geom)
  plot(bnd$geom, lwd=2, border='black', add=TRUE)
  plot(scan_locations$geom, add=TRUE, col=rgb(0,0,1,0.2))
  Sys.sleep(0.5)
  
  # run through grid tiles, load proximal TLS scans from directory and clip to bnd. rbind, and write to file.
  for(t in 1:nrow(grid)) {
    # Load and display tile
    cat('loading tile', t, 'of', nrow(grid), '\n')
    tile = grid[t,]
    scans_to_load = scan_locations[sf::st_intersects(tile, scan_locations, sparse = FALSE),]
    plot(grid$geom)
    plot(bnd$geom, add=TRUE, lwd=2)
    plot(scan_locations$geom, add=TRUE, col=rgb(0,0,1,0.2))
    plot(grid[1:t,]$geom, add=TRUE, col='grey')
    plot(tile,add=TRUE, col='yellow')
    
    #loop through relevant scans, clip and 
    combined_las = list()
    ex = round(sf::st_bbox(tile))
    filt = paste('-keep_xy', ex[1], ex[2], ex[3], ex[4]) #min_x min_y max_x max_y
    for(i in 1:nrow(scans_to_load)) {
      cat('.....appending scan', i, 'of', nrow(scans_to_load), '\n')  
      roi = sf::st_intersection(tile, scans_to_load[i,])
      combined_las[[i]] = lidR::clip_roi(lidR::readTLSLAS(scans_to_load[i,]$fn, filter=filt), roi)
    }
    
    #To avoid rbind errors, Find common columns among scans and keep only those
    common_cols = lapply(combined_las, function(x) colnames(x@data))
    common_cols = Reduce(intersect, common_cols)
    combined_las = lapply(combined_las, function(x) {
      x@data = x@data[, c('X', 'Y', 'Z', 'gpstime', 'Intensity', 'ReturnNumber', "NumberOfReturns", 'Classification', 'Reflectance', 'Deviation')]
      return(x)})
    combined_las = do.call(rbind,combined_las)
    combined_las@header@VLR = list()

    #write tile to disk
    cat('.....scans stitched. writing tile to disk')
    out_las = paste0(out_dir, '/', ex[1], '_', ex[2], '.las')
    lidR::writeLAS(las_update(combined_las), out_las, index=index)
  }
  return(NULL)
}
