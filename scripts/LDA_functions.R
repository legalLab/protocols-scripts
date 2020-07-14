#!/usr/bin/env Rscript

# functions to process LDA for visualization
# Tomas Hrbek July 2019

# load required libraries
library(dplyr)
library(jamba)

# function - generate list of species occuring within each cell
spCellPresence <- function(rstr, path) {
  # generate list of all cells (including NA), converto to tibble
  arc_cells <- as_tibble(paste('cell', seq(1, length(values(rstr))), sep='_', collapse=NULL))
  names(arc_cells)[1] <- 'cell'
  
  # set up a loop
  for (i in 1:length(path)) {
    # read species polygon
    sp_p <- readOGR(path[i])
    # add projection only if missing or different from base raster
    #proj4string(sp_p) <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
    # mask species polygon, absence=0 (updatevalue=0)
    r <- mask(x=rstr, mask=sp_p)
    # add species to list and change column names to species
    arc_cells <- add_column(arc_cells, values(r) %>% replace(!is.na(.), sp_names[i]))
    names(arc_cells)[i+1] <- sp_names[i]
    }
  
  # remove all NA columns - species that do not occur within an area
  arc_cells <- arc_cells %>% select_if(function(x){!all(is.na(x))})
  
  return(arc_cells)
}

# function - generate cell colors based on relative contributions of ecoregions, add to tibble
cellColors <- function(cell_ecoregions, ecoregions_colors) {
  cell_colors <- cell_ecoregions[,2:ncol(cell_ecoregions)] %>% 
    as.matrix() %*% t(col2rgb(ecoregions_colors)) %>% 
    rgb2col()
  
  # add color of NA cells
  cell_colors <- c("#FFFFFF00", cell_colors)
  
  return(cell_colors)
}

# function - generate cell shade based on relative occurance of ecoregion
cellShades <- function(cell_ecoregions, ecoregions_colors, cell_colors, i) {
  # remove first color "#FFFFFF00 representing NA
  cell_colors <- cell_colors[-1]
  
  z <- ecoregions_colors[i] %>% 
    col2rgb(alpha = FALSE) %>% 
    t()
  z2 <- matrix(z, nrow = length(cell_colors), ncol = 3, byrow = TRUE)
  z3 <- trunc(255 * cell_ecoregions[i+1])
  cell_shades <- cbind(z2, z3) %>% 
    as.matrix() %>% 
    rgb2col()
  
  # add color of NA cells
  cell_shades <- c("#FFFFFF00", cell_shades)
  
  return(cell_shades)
}


# function - harmonic mean
harmonicMean <- function(logLikelihoods) {
  1/mean(1/logLikelihoods)
}

# function geometric mean
geometricMean <- function(logLikelihoods) {
  prod(logLikelihoods)^(1/length(logLikelihoods))
}

# examples
