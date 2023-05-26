#' Clip a LAS object to a volume of interest
#'
#' This code is a 3-dimensional version of `lidR::clip_roi` which interprets
#' a smaller las object as a volume of interest and clips a larger LAS object
#' to that volume. Some segmentation algorithms result in lower resolution
#' LAS objects than the original files, so this function can be used to clip
#' a low resolution segmentated LAS object and retrieve the volume specified from
#' a higher resolution LAS.
#' @param las a `LAS` usually with a large extent and/or high point density
#' @param voi a `LAS` object often with lower resolution and extent from which
#' to generate a volume to clip `las`
#' @param `res` the resolution at which the `voi` is generated
#' @examples
#' library(landecoutils)
#' library(lidR)
#'
#' # Large full resolution LAS object
#' full_las = readTLSLAS('path/to/highres.las')
#'
#' # smaller lower resolution LAS object
#' voi = readTLSLAS('path/to/lowres.las')
#'
#' # view point density an extent of each las
#' nrow(full_las@data);extent(full_las)
#' nrow(voi@data);extent(voi)
#'
#' # Clip voi and view point density and extent
#' high_res_voi = clip_voi(full_las, voi, res=0.5)
#' nrow(high_res_voi@data); extent(high_res_voi)
#' @export
clip_voi = function(las, voi, res){
  vox = lidR::voxelize_points(voi, res)
  bbox = as.numeric(sf::st_bbox(voi))
  high_res = lidR::filter_poi(las, X >= bbox[1] & X <= bbox[3] & Y >= bbox[2] & Y <= bbox[4])
  grid_temp = terra::rast(terra::ext(vox))
  terra::res(grid_temp) = c(res, res)
  slices = list()
  for(z in seq(min(vox$Z), max(vox$Z) - res, res)) {
    slice_vox = lidR::filter_poi(vox, Z == z)
    hull = lidR::grid_metrics(slice_vox, ~length(Z)>0, res = res, start=terra::ext(grid_temp)[c(1,3)])
    hull = terra::rast(hull)
    hull = terra::as.points(terra::as.polygons(hull))
    hull = sf::st_as_sf(terra::convHull(hull))
    slice_highres = suppressWarnings(lidR::clip_roi(lidR::filter_poi(high_res, Z < (z+res) & Z > (z-res)), hull))
    slices[[length(slices)+1]] = slice_highres
  }
  slices=do.call(rbind, slices)
  return(slices)
}
