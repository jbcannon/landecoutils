# landecoutils
This R package contains functions to assist in common tasks in the Jones Center Landscape Ecology lab. Current functionality includes
1. tiling terrestrial lidar scans (`stitch_TLS_dir_to_LAS()` and `stitch_TLS_dir_to_LAS_tiles()`
2. stem mapping terrestrial lidar scans (`segment_with_spanner()`)
3. clipping las to 3d convext hull (`clip_roi()`)

These functions are primarily built for my personal use and thus have undergone minimal testing and error-proofing. The only guarantee is that you'll likely run into some issues if you try to use them :)

# Install `landecoutils`
If you haven't installed packages from Github before or used the devtools package, you'll need to first ensure that Rtools is properly installed. It is installed separately from typical packages.

Install the correct version of Rtools using this link: [https://cran.r-project.org/bin/windows/Rtools/rtools40.html](https://cran.r-project.org/bin/windows/Rtools/rtools40.html)
Be sure to follow the directions at the bottom regarding Putting Rtools on the PATH. This is critical for a proper installation. You'll need to create a file in your Documents folder named '.Renviron' that contains the following line: `PATH="${RTOOLS40_HOME}\usr\bin;${PATH}"` and restart R.
Next, you should be able to get the latest released version of hurrecon from github

```
install.packages('devtools')
devtools::install_github('jbcannon/landecoutils')
```
