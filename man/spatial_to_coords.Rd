% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/scan_pts.R
\name{spatial_to_coords}
\alias{spatial_to_coords}
\title{Covert .kml or .shp to gps coordinates}
\usage{
spatial_to_coords(in_file, csv_out = NULL, proj = 32616)
}
\arguments{
\item{in_file}{path to a .KML or .SHP containing points}

\item{csv_out}{path to a .CSV or .TXT to output points. If not provided, an
interactive file chooser will be used}

\item{proj}{epsg projection for output. Defaults to 32316 (UTM 16N)}
}
\description{
This function takes a kml or other shapefile, transforms them to
WGS84, and export GPS coordinates
}
\examples{
in_file = 'R:/landscape_ecology/test.shp'
out_file = 'R:/landscape_ecology/output.csv'
spatial_to_coords(in_file, out_file)

# Use an interactive file chooser
in_file = file.choose()
out_file = file.choose()
spatial_to_coords(in_file, out_file)

}
