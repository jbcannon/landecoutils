#' Create points for scan locations around a central point
#'
#' This function takes a shapefile containing plot locations,
#' and creates new points in the specified directions at a constant
#' distance. It returns the original points (newid=0) plus the
#' surrounding pionts (newid > 0)
#' @param input_shp path to a .SHP containing points
#' @param output_shp path to a .SHP to output points
#' @param dist distance (in meters) of desired distance to each point
#' @param dirs directions of desired new points (in degrees)
#' @examples
#' in_file = 'R:/landscape_ecology/test.shp'
#' out_file = 'R:/landscape_ecology/output.shp'
#' dirs = c(0,90,180,270)
#' dist = 20
#' make_scan_points(in_file, out_file, dist, dirs)
#'
#' @export
make_scan_points = function(input_shp, output_shp, dist, dirs=c(0,90,180,270)) {
  if(!file.exists(input_shp)) stop('input file not found')
  if(!dir.exists(dirname(output_shp))) stop('output directory does not exist')
  if(is.null(dist)) stop('dist must be provided')
  if(is.null(dirs)) stop('dirs must be provided')

  shp = sf::read_sf(input_shp) # read in points
  orig_trans = sf::st_crs(shp) # get original coordinate system
  shp = sf::st_transform(shp, 32616) #transform to UTMs
  coords = sf::st_coordinates(shp) #etract coordinates

  # Loop through all points and get new point in all directions

  # function to output new coordinate given original
  my_function = function(x, y, dir, dist) {
    rads = dir * pi/180 #convert to radians
    x_dist = sin(rads)*dist
    y_dist = cos(rads)*dist
    new_x = x+x_dist
    new_y = y+y_dist
    return(data.frame(X=new_x, Y=new_y))
  }

  all_out = list()
  for(i in 1:nrow(coords)) { # loop through all input coords, i
    pt_out = list()
    # apply function for each direction
    for(d in dir) {  # loop through all directions, d
      new_pt = my_function(coords[i,'X'], coords[i,'Y'], d, dist)
      pt_out[[length(pt_out)+1]] = new_pt
    }
    pt_out = do.call(rbind, pt_out) #combine outputs
    #attach old column info to new pts
    colnames = colnames(st_drop_geometry(shp))
    for(c in colnames) {
      pt_out[,c] = NA
      pt_out[,c] = sf::st_drop_geometry(shp[i,c])
    }
    # Attach new points to old point
    pt_out = sf::st_as_sf(pt_out, coords=c('X', 'Y'), crs=32616)
    pt_out = sf::st_transform(pt_out, orig_trans)
    pt_out = rbind(sf::st_transform(shp[i,],orig_trans), pt_out)
    #add new id number
    pt_out$newid = 0:(nrow(pt_out)-1)
    all_out[[length(all_out)+1]] = pt_out
  }
  all_out = do.call(rbind, all_out)
  sf::write_sf(all_out, output_shp)
  cat('output file written to:\n', output_shp)
  return(all_out)
}


