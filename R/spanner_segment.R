#' Create a stem map from a TLS scan tile using spanner.
#' 
#' This is code to make it more convenient to pre-process, stem map, and 
#' segment TLS tiles. The function takes a merged TLS scan and a boundary
#' and generates a stem map (.shp) attributed with stem radius, and height. 
#' You may optionally output a LAS with the segmented trees by specifying 
#' a path name (*.las) in the argument `tree_seg_las`. If left NULL, only the 
#' stem map will be written to disk. Note that there are large memory requirements
#' and you may encounter memory issues if you use this on areas
#' larger than 0.2 ha with high scan densities. 
#' @param las a `LAS` object not much larger than 0.2 ha
#' @param bnd an `sf` object representing the *interior* boundary of the study area.
#' Generally `bnd` is smaller than the `bbox` of `las` to adjust edge effects
#' @param stemmap_shp a path name (.shp) with the desired output of the resulting stem map.
#' @param tree_seg_las a path name (.las) for the output of segmented trees
#' @param threads see `lidR::set_lidR_threads`. Positive scalar. Default 0 means use all
#'  CPU available. Values > 1 mean using n cores, values in (0, 1) mean using a 
#'  fraction of the cores e.g. 0.5 = half.
#'  @export
segment_with_spanner = function(las, bnd, stemmap_shp, tree_seg_las=NULL, threads=0.5){
  # Load and pre-process las
  las = lidR::readLAS(in_las, filter='-thin_with_voxel 0.02')
  las = lidR::classify_ground(las, lidR::csf(class_threshold=0.1))
  dem = lidR::grid_terrain(las, res=0.1, lidR::tin())
  las = lidR::normalize_height(las, lidR::tin())
  las = lidR::filter_poi(las, Classification != 2)
  
  #Find tree locations and segment tree points with spanner functions 
  boles = lidR::filter_poi(las, Intensity>40000)
  set_lidr_threads(threads)
  tree_locs = spanner::get_raster_eigen_treelocs(boles, res=0., pt_spacing=0.01,
                                                 dens_threshold = 0, eigen_threshold = 0.9,
                                                 grid_slice_min=1, grid_slice_max=2,
                                                 minimum_polygon_area = 0.01)
  las = spanner::segment_graph(las, tree_locs,
                               k = 50,
                               distance.threshold = 0.5,
                               subsample.graph = 0.2, 
                               return.dense = TRUE)
  # Create stemmap and clip to interior boundary
  las = lidR::add_lasattribute(las, las$treeID, 'TreeID', 'spanner TreeID')
  las = remove_lasattribute(las, 'treeID')
  las = TreeLS::stemPoints(las, method = TreeLS::stm.eigen.voxel())
  inv = TreeLS::tlsInventory(las)
  stemmap = sf::st_as_sf(inv, coords = c('X', 'Y'), crs=sf::st_crs(las))
  st_crs(bnd) = sf::st_crs(las)
  stemmap = subset(stemmap, !is.na(Radius))
  stemmap = sf::st_intersection(stemmap, bnd)
  sf::st_write(stemmap, stemmap_shp, delete_layer=TRUE)
  
  #Cleanup and output las with tree segmentation (optional)
  if(!is.null(tree_seg_las)) {
    clipout_trunks = sf::st_buffer(stemmap, dist=stemmap$Radius*2)
    las = lidR::classify_poi(las, lidR::LASLOWVEGETATION, roi = clipout_trunks, poi = formula(~Z<1.5), inverse_roi=TRUE)
    tree_seg = lidR::filter_poi(las, TreeID %in% stemmap$TreeID & Classification != lidR::LASLOWVEGETATION)
    tree_seg = lidR::las_update(tree_seg)
    lidR::writeLAS(tree_seg, tree_seg_las)
  }
  return(stemmap)
}

