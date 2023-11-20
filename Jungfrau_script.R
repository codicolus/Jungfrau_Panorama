# --------------------------------------------------------------------------------------------------------------------- #
# --------------------------------------------------------------------------------------------------------------------- #
#                                               ALETSCH-ARENA MOVIE                                                     #
#                                                                                                                       #
#                    Script for rayrendering sequence of images used for compiling a short                              #
#                                              360° panorama movie                                                      #
#                                                                                                                       #
#                   Required data is provided by (c) Swissstopo under CC-BY                                             #
#                   1) SwissALTI3D:  https://www.swisstopo.admin.ch/de/geodata/height/alti3d.html                       #
#                   2) SwissIMAGE: https://www.swisstopo.admin.ch/de/geodata/images/ortho/swissimage10.html             #
#                                                                                                                       #
#                                          Created by: Christoph von Matt                                               #
#                                                   15-03-2022                                                          #
#                                          Licence (Script): GNU GPL v3.0                                               #
#                                                                                                                       #
#                 Some steps follow or are adapted from Tutorials or rayshader-function descriptions                    #
#                                         provided by Tyler Morgan-Wall                                                 #
#                                                                                                                       #
#                              Tiles were manually merged & cropped using Quantum GIS.                                  #
#                                                                                                                       #
# --------------------------------------------------------------------------------------------------------------------- #
# --------------------------------------------------------------------------------------------------------------------- #

# required libraries
library(raster)
library(rayshader)
library(ggplot2)
library(sf)
library(tidyverse)

# load DEM + SwissIMAGE
dem <- raster("temp/merged_large.tif")
image_cropped <- stack("temp/merged_image_large.tif")

#plot(dem)
#plotRGB(image)

# get rgb channels and rescale for contrast
image_cropped <- image_cropped[[1:3]]
names(image_cropped) <- c("r", "g", "b")
image_array <- as.array(image_cropped) / 255
image_array_contrast <- scales::rescale(image_array, to = c(0,1))


# Custom assembled control points (using rgl window + freehand-sketch)
# -------------------------------------------------------------------
# --> by now probably more convenient by using "render_path"
# see: https://rdrr.io/cran/rayshader/man/render_path.html

# data frame with control points
df <- as.data.frame(do.call(rbind, list(
  c(90, 0, 1.19),
  c(90, 0, 0.8),
  c(90, 0, 0.359),
  c(90, 0, 0.236),
  c(45, 0, 0.359),
  c(30.00272, -40, 0.3972), # evtl 25? - first turn
  c(30.00272, -45, 0.24095), # prolong a bit! (in final df? x-steps)
  c(30.00272, -45, 0.3972),
  c(30.00272, -70, 0.3972), # Helper points to make turns slower (30 steps)
  c(30.00272, -100, 0.3972), # Helper points to make turns slower (30 steps)
  c(30.00272, -130, 0.3972), # Helper points to make turns slower (30 steps)
  c(30.00272, -135, 0.3972), # second fast turn (too fast!)
  c(20, -135, 0.1895), # prolong a bit (e.g. final df? x-steps)
  c(10, -135, 0.1895), # phi-change point
  c(20, -135, 0.1895),
  c(20, -135, 0.33),
  c(20, -160, 0.33), # Helper point for slower turing
  c(20, -190, 0.33), # Helper point for slower turing
  c(25, -220, 0.33), # 0.2884
  # c(25, -200, 0.2884), # original point - must eventually be restored?
  c(25, -235, 0.2884),
  c(20, -235, 0.2315), # before: zoom = 0.25078 but now more + prolong a bit (final df, x-steps)
  c(20, -265, 0.2315), # Helper point for slower turning
  c(20, -295, 0.2315), # phi evtl. 30, zoom evtl 0.2269
  c(20, -310, 0.2315), # helper point for slower turning
  c(25, -325, 0.2315), # phi evtl. 30, zoom evtl. 0.2269
  c(45, -345, 0.8),
  c(90, -360, 1.19)
)))
colnames(df) <- c("phi", "theta", "zoom")

# get camera vectors
phi <- df$phi
theta <- df$theta
zoom <- df$zoom
nums <- seq_along(zoom)
df <- df %>% 
  as_tibble() %>% 
  mutate(n = row_number()) %>% 
  select(n, everything()) %>% 
  pivot_longer(cols = 2:4)

# quick look at camera vectors
ggplot(df) +
  geom_point(aes(n, value)) +
  geom_smooth(aes(n, value), method = "gam") +
  facet_wrap(~name, ncol = 1, scales = "free_y")

# -----------------------------------------------------------------------------
# desired number of images / timesteps
timesteps = 1800
# -----------------------------------------------------------------------------

# Smoothing camera paths following manual control points
# ------------------------------------------------------
# check theta-fit
fit <- lm(theta~poly(nums, 15))
plot(nums, theta, pch = 15)
theta_smooth <- spline(nums, theta, n = timesteps, method = "hyman") # natural, fmm, periodic, monoH.FC, hyman
lines(theta_smooth$x, theta_smooth$y, col = "red", lwd = 2)
lines(thetavec$x, thetavec$y, col = "blue", lwd = 2)

# check zoom fit
fit <- lm(zoom~poly(nums, 15))
plot(nums, zoom, pch = 15)
zoom_smooth <- spline(nums, zoom, n = timesteps, method = "periodic") # natural, fmm, periodic, monoH.FC, hyman
lines(zoom_smooth$x, zoom_smooth$y, col = "red", lwd = 2)
lines(zoomvec$x, zoomvec$y, col = "blue", lwd = 2)

# check phi fit
fit <- lm(phi~poly(nums, 15))
plot(nums, phi, pch = 15)
phi_smooth <- spline(nums, phi, n = timesteps, method = "periodic") # natural, fmm, periodic, monoH.FC, hyman
lines(phi_smooth$x, phi_smooth$y, col = "red", lwd = 2)
lines(phivec$x, phivec$y, col = "blue", lwd = 2)


# Function to create a combined version of both smoothed fits for smoother
# transitions (red + blue lines)
flatten_fits <- function(orig_x, orig_vals, fit_x, fit_y){
  out_x <- c()
  out_y <- c()
  for(i in seq_along(orig_x)){
    # if last one skip
    if(i == length(orig_x)){
      out_x <- c(out_x, fit_x[length(fit_x)])
      out_y <- c(out_y, fit_y[length(fit_y)])
    }else{
      # get current
      temp_x1 <- orig_x[i]
      temp_y1 <- orig_vals[i]
      # get next values
      temp_x2 <- orig_x[i+1]
      temp_y2 <- orig_vals[i+1]
      
      # get contorol and follow up points
      # 1) if fit and orig equal - take orig x-position
      # 2) if not equal take fit x-position
      indizes <- which(fit_x >= temp_x1 & fit_x < temp_x2)
      values_x <- fit_x[indizes]
      values_y <- fit_y[indizes]
      
      # case if y2 == y1 then all values the same as x2 for orig_vals
      if(temp_y2 == temp_y1){
        values_y <- rep(temp_y2, length(indizes))
      }
      out_x <- c(out_x, values_x)
      out_y <- c(out_y, values_y)
    }
  }
  return(data.frame(x = out_x, y = out_y, stringsAsFactors = FALSE))
}

# combined fit + orig camera vectors
phivec <- flatten_fits(nums, phi, phi_smooth$x, phi_smooth$y)
thetavec <- flatten_fits(nums, theta, theta_smooth$x, theta_smooth$y)
zoomvec <- flatten_fits(nums, zoom, zoom_smooth$x, zoom_smooth$y)
phivec_custom <- phivec$y
theta_custom <- thetavec$y
custom_zoom <- zoomvec$y

# smoothed camera vector dataset
df <- data.frame(
  n = seq_along(theta_custom),
  phi = phivec_custom,
  theta = theta_custom,
  zoom = custom_zoom,
  stringsAsFactors = FALSE
) %>% 
  pivot_longer(cols = 2:4)

# quick look at camera vectors
ggplot(df) + 
  geom_line(aes(n, value)) +
  facet_wrap(~name, ncol = 1, scales = "free_y")

# Prolong at stopping points for x steps
# -------------------------------------
# TODO: evtl adjust later if fps = 60
howmany_prolong = 30 # 30 = 1s
df <- df %>% 
  pivot_wider(names_from = name, values_from = value) %>% 
  as.data.frame()
# Point 1 at zoom = 0.24
temp_df <- df[df$n > 400 & df$n < 500,]
temp_nums <- temp_df$n
temp_row <- temp_nums[which(temp_df$zoom == min(temp_df$zoom))]
rm(temp_nums)
prolong_df <- do.call(rbind, lapply(seq(howmany_prolong), function(x, temp_row, data){return(data[temp_row,])},
                                    temp_row = temp_row, data = df))
df <- rbind(df[1:temp_row,], prolong_df, df[(temp_row + 1):dim(df)[1],])
df <- df %>% mutate(n = row_number())


# get vectors again
phivec_custom <- df$phi
theta_custom <- df$theta
custom_zoom <- df$zoom

render_movie(filename = "movie/custom_testlong_adjusted", type = "custom",
             frames = length(custom_zoom) , fps = 30, 
             phi = phivec_custom, zoom = custom_zoom, theta = theta_custom)
