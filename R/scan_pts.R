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
    for(d in dirs) {  # loop through all directions, d
      new_pt = my_function(coords[i,'X'], coords[i,'Y'], d, dist)
      pt_out[[length(pt_out)+1]] = new_pt
    }
    pt_out = do.call(rbind, pt_out) #combine outputs
    #attach old column info to new pts
    colnames = colnames(sf::st_drop_geometry(shp))
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

#' Covert .kml or .shp to gps coordinates
#'
#' This function takes a kml or other shapefile, transforms them to
#' WGS84, and export GPS coordinates
#' @param in_file path to a .KML or .SHP containing points
#' @param csv_out path to a .CSV or .TXT to output points. If not provided, an
#' interactive file chooser will be used
#' @param proj epsg projection for output. Defaults to 32316 (UTM 16N)
#' @examples
#' in_file = 'R:/landscape_ecology/test.shp'
#' out_file = 'R:/landscape_ecology/output.csv'
#' spatial_to_coords(in_file, out_file)
#'
#' # Use an interactive file chooser
#' in_file = file.choose()
#' out_file = file.choose()
#' spatial_to_coords(in_file, out_file)
#'
#' @export
spatial_to_coords = function(in_file, csv_out = NULL, proj=32616) {
  # load and transform kml
  if(!file.exists(in_file)) stop('input file does not exist')
  if(!is.null(csv_out)) {
    if(file.exists(csv_out)) stop('output file already exists')
    if(tools::file_ext(csv_out) %in% c('txt', 'csv')) stop('output file must be *.txt or *.csv')
  }
  kml = sf::read_sf(in_file)
  kml = sf::st_transform(kml, sf::st_crs(proj))
  coords = sf::st_coordinates(kml)[,c('X', 'Y')]
  coords = as.data.frame(coords)
  coords$proj =   suppressWarnings(sf::st_crs(proj)$input)
  colnames(coords) = c('X', 'Y', 'proj')
  kml = sf::st_drop_geometry(kml)
  kml = cbind(kml, coords)
  print(knitr::kable(kml))
  if(is.null(csv_out)) {
    readr::write_csv(kml, file.choose(new=TRUE))
  } else {
    readr::write_csv(kml, csv_out)
  }
  return(kml)
}
