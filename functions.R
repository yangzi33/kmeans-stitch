### STA314 Assignment 1: Cross-Stitch ###
# Author: Ziyue Yang                    #
# Student Number: 1004804759            #
# Contact: ziyue.yang@mail.utoronto.ca  #
#########################################

# Install the required libraries
# Uncomment the following if you don't have them
# install.packages("imager")
# install.packages("tidyverse")
# install.packages("tidymodels")
# install.packages("sp")
# install.packages("scales")
# install.packages("cowplot") 
# devtools::install_github("sharlagelfand/dmc")

# Loading libraries
library(imager)
library(tidyverse)
library(tidymodels)
library(sp)
library(scales)
library(cowplot) 
library(dmc)

# From change_solution.R
change_resolution <- function(image_df, x_size)
{
  ## change_resolution(image_df, x_size) subsamples an image to produce
  ## a lower resolution image. Any non-coordinate columns in the data
  ## frame are summarized with their most common value in the larger
  ## grid cell.
  ##
  ## Input:
  ## - image_df: A data frame in wide format. The x-coordinate column MUST
  ##             be named 'x' and the y-coordinate column MUST be named 'y'.
  ##             Further columns have no naming restrictions.
  ## - x_size:   The number of cells in the x-direction. The number of cells
  ##             in the vertical direction will be computed to maintain the 
  ##             perspective. There is no guarantee that the exact number
  ##             of cells in the x-direction is x_size
  ##
  ## Output:
  ## - A data frame with the same column names as image_df, but with fewer 
  ##   entries that corresponds to the reduced resolution image.
  ##
  ## Example:
  ##   library(imager)
  ##   library(dplyr)
  ##   fpath <- system.file('extdata/Leonardo_Birds.jpg',package='imager') 
  ##   im <- load.image(fpath)
  ##   im_dat<- as.data.frame(im,wide = "c") %>% rename(R = c.1, G = c.2, B = c.3) %>%
  ##            select(x,y,R,G,B)
  ##   agg_image <- change_resolution(im_dat, 50)
  
  if(!require(sp)) {
    stop("The sp packages must be installed. Run install.packages(\"sp\") and then try again.")
  }
  if(!require(dplyr)) {
    stop("The dplyr packages must be installed. Run install.packages(\"dplyr\") and then try again.")
  }
  
  sp_dat <- image_df 
  gridded(sp_dat) = ~x+y
  
  persp = (gridparameters(sp_dat)$cells.dim[2]/gridparameters(sp_dat)$cells.dim[1])
  y_size = floor(x_size*persp)
  orig_x_size = gridparameters(sp_dat)$cells.dim[1]
  orig_y_size = gridparameters(sp_dat)$cells.dim[2]
  
  x_res = ceiling(orig_x_size/x_size)
  y_res = ceiling(orig_y_size/y_size)
  
  gt = GridTopology(c(0.5,0.5), c(x_res, y_res),
                    c(floor(orig_x_size/x_res), floor(orig_y_size/y_res)))
  SG = SpatialGrid(gt)
  agg = aggregate(sp_dat, SG, function(x) names(which.max(table(x)))[1] )
  agg@grid@cellsize <- c(1,1)
  df <- agg %>% as.data.frame %>% rename(x = s1, y = s2)  %>% select(colnames(image_df))
  
  return(df)
}

# Function a
process_image <- function(image_file_name, k_list) {
  ###########################################################################
  ## process_image gets k-cluster info of input image, where k is the
  ## numbers of cluster centers in input list k_list
  ## Input:
  ##        - image_file_name: path of the image to cluster as a string
  ## 
  ##        - k_list: list of cluster centers. For instance, if we want
  ##                  to cluster the image with 1, 2, 3 centers respec-
  ##                  tively, then k_list will be the list c(1:3). 
  ## 
  ## Output: A tibble `cluster_info` containing cluster information for 
  ##         each k in k_list. For each k, cluster_info[[k]] contains the 
  ##         following info:
  ##
  ##        - k_clust: the original output of the kmeans() call on image with k.
  ## 
  ##        - glanced: the glance of cluster_info[[k]], which will be used to 
  ##                   produce the scree plot.
  ## 
  ##        - tidied: tidied clusters info, with their associated RGB value and
  ##          the DMC value closest to the RGB color.
  ## 
  ###########################################################################
  
  # Loading image and store color info as data frame image_dat
  img <- imager::load.image(image_file_name)
  image_dat <- as.data.frame(img, wide = "c") %>% 
    rename(R = c.1, G = c.2, B = c.3)

  kclusts <- 
    tibble(k = k_list) %>%
    mutate(
      kclust = map(k, ~kmeans(x = select(image_dat, c(-x, -y)), centers = .x, nstart = 5)),
      glanced = map(kclust, glance),
      tidied = map(k, ~(
        tidy(kclust[[.x]]) %>% mutate(RGB = rgb(R, G, B),
                                  # Initialize the DMC column with null entries
                                      DMC = vector("character", length = .x)))
        )
    )

  # Adding DMC
  for (k in k_list) {
    for (cen in 1:length(kclusts$tidied[[k]]$RGB)) {
      kclusts$tidied[[k]]$DMC[cen] <- dmc(kclusts$tidied[[k]]$RGB[cen])$hex
    }
    kclusts$kclust[[k]] <- augment(kclusts$kclust[[k]], image_dat)
  }
  
  return(kclusts)
}

# Function b
scree_plot <- function(cluster_info) {
  ###########################################################################
  ## Produces the scree plot of cluster_info - the output of `process_image`. 
  ## Input:
  ##        - cluster_info: the output of `process_image`. Details can be found
  ##                        in documentation of function A: process_image().
  ## 
  ## Output: 
  ##        -A scree plot based on `glanced` for each k in cluster_info.
  ## 
  ###########################################################################
  ggplot(cluster_info %>% unnest(cols = c(glanced)),
         aes(x = k, y = tot.withinss)) + geom_line() + geom_point()
}


# Function C
colour_strips <- function(cluster_info) {
  ###########################################################################
  ## For each k in cluster_info$k, produces a DMC colour strip.
  ## Input:
  ##        - cluster_info: the output of `process_image`. Details can be found
  ##                        in documentation of function A: process_image().
  ## 
  ## Output: 
  ##        - k colour strips. 
  ##          The i^th colour strip ( where 1 <= i <= length(k_list)) contains
  ##          the closest DMC colour of centers
  ## 
  ###########################################################################
  for (k in cluster_info$k) {
    n_col = length(cluster_info$tidied[[k]]$DMC)
    rect_dat <- tibble(x1 = c(0: (n_col - 1)), x2 = c(1: n_col), y1 = rep(0, n_col),
                       y2 = rep(1, n_col), colour = cluster_info$tidied[[k]]$DMC)
    k_strip <- ggplot(rect_dat) + coord_fixed() +
      geom_rect(aes(xmin = x1, xmax = x2, ymin = y1, ymax = y2, fill = colour)) + 
      geom_text(aes(x = x1 + (x2 - x1) / 2, y = y1 + (y2 - y1) / 2, label = colour),
                size = 24 / n_col, colour = "white") + 
      scale_fill_manual(values = rect_dat$colour) + theme_void() + 
      theme(legend.position = "none")
  }
  plot_grid(k_strip)
}


# Function D
make_pattern <- function(cluster_info, k, x_size,
                         black_white = FALSE, background_colour = NULL) {
  ###########################################################################
  ## Make a cross-stitch pattern based on image inputted to `process_image`
  ## Inputs:
  ##        - cluster_info: Output of `process_image()`
  ##
  ##        - k: The chosen number of clusters
  ##
  ##        - x_size: The approximate total number of possible stitches in the
  ##                  horizontal direction.
  ##
  ##        - black_white: Logical (or boolean) value that indicates whether 
  ##                       the cross-stitch produced is colourful (FALSE) or
  ##                       black-and-white (TRUE). Default is set to FALSE.
  ##
  ##        - background_colour: Background colour of the cross-stitch plot. 
  ##                             Default is set to NULL, which produces a 
  ##                             transparent background. Input takes colour
  ##                             names such as "blue", or colour codes such 
  ##                             as "#0000FF".
  ##
  ## Output: The cross-stitch pattern of image input with k centers. The pattern
  ##         is "stitched" on surface with colour of input `background_colour`.
  ## 
  ###########################################################################
  
  # Resize the image
  resize_image_info <- cluster_info$kclust[[k]] %>% select(x, y, R, G, B, .cluster)
  resized_image <- change_resolution(resize_image_info, x_size)
  
  dmc_names <- c()
  for (i in 1:k) {
    dmc_names <- append(dmc_names, dmc(cluster_info$tidied[[k]]$RGB[i])$name)
  }
  
  if (!black_white) {  # Colourful
    resized_image %>% ggplot(aes(x = x, y = y)) + 
      geom_point(aes(shape = .cluster,
                     col = .cluster)) + 
      scale_y_reverse() + theme_void() +
      scale_colour_manual(values = cluster_info$tidied[[k]] %>% select(DMC) %>% deframe,
                          label = cluster_info$kclust[[k]] %>% select(.cluster) %>% deframe,) +
      theme(panel.background = element_rect(fill = background_colour, colour = background_colour)) +
      labs(shape = "Cluster", col = "DMC Colour")
  } else { # Black and white
    resized_image %>% ggplot(aes(x = x, y = y)) + 
      geom_point(aes(shape = .cluster)) + 
      scale_y_reverse() + theme_void() +
      theme(panel.background = element_rect(fill = background_colour, colour = background_colour)) + 
      labs(shape = "Cluster")
  }
}
