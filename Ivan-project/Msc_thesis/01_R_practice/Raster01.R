rm(list=ls(all=TRUE))
####################################################################
##THE PURPOSE OF THIS SCRIPT IS TO PRACTICE CREATING RANDOM RASTER## 
###################################################################
library(raster)

###################
#Create a random spatial grid
#Run library(raster)

#Output file names
fout="randomgrid_practice_01.csv"

#Define extent

xmin = 0
ymin = 0
xmax = 25
ymax = 25
aoi_extent1 = extent(xmin, xmax, ymin, ymax)

#Create grid and get xy coordinates from which to extract these random points

rast1 = raster(ext=aoi_extent1, resolution=1)
values(rast1) = 1:ncell(rast1)

#Extract cell centre coordinates

x_centres=xFromCol(rast1)
y_centres=yFromRow(rast1)

#Select some random points
random_point_count = 0
random_point_sample_number = 100 #the number of points you want

while(random_point_count < random_point_sample_number) {
  
  easting_random = sample(x_centres, 1)
  northing_random = sample(y_centres, 1)
  
  if(exists("random_points")==TRUE){
    
    #random_points dataframe already exists 
      random_points = rbind(random_points,c("EASTING"=easting_random,
                                            "NORTHING"=northing_random))
  }else{
    #random_points dataframe doesn't exists yet...
    random_points = data.frame("EASTING"=easting_random, "NORTHING" = northing_random)
  }
  
  random_point_count=random_point_count + 1
}

write.csv(random_points, fout, row.names = FALSE)

#Plotting#
library(ggplot2)
library(reshape2)

#Output file name
plot_out = "random_points_CELL_CENTRES_1.png"

#Convert raster to matrix

mat1 = as.matrix(rast1)
dat1 = melt(mat1)
names(dat1) = c("easting","northing","value")

#Shift over the easting and northing values by half pixel (ggplot uses coordinates as midpoints, whereas we want it to plot from the mid-point minus half a pixel lenght or width)

dat1$easting = dat1$easting - 0.5
dat1$northing = dat1$northing - 0.5

p = ggplot(dat1, aes(x = easting, y = northing)) +
          geom_tile(aes(fill=value), colour="grey20") +
          scale_fill_gradientn(colours = terrain.colors(10)) +
          geom_point(data=random_points, mapping = aes(x=EASTING, y=NORTHING), colour="black") +
          labs(x = "Easting", y = "Northing") #+ 
          #theme(legend.position = "none")
#p

ggsave(plot_out, plot = p, height = 6, width = 7, dpi = 200)



####INCREASE SIZE####

#Output file names
fout="randomgrid_practice_02.csv"

#Define extent

xmin = 0
ymin = 0
xmax = 100
ymax = 100
aoi_extent2 = extent(xmin, xmax, ymin, ymax)

#Create grid and get xy coordinates from which to extract these random points

rast2 = raster(ext=aoi_extent2, resolution=1)
values(rast2) = 1:ncell(rast2) #rnorm to do a random grit too!!!

#Extract cell centre coordinates

x_centres=xFromCol(rast2)
y_centres=yFromRow(rast2)

#Select some random points
random_point_count = 0
random_point_sample_number = 100 #the number of points you want

while(random_point_count < random_point_sample_number) {
  
  easting_random = sample(x_centres, 1)
  northing_random = sample(y_centres, 1)
  
  if(exists("random_points")==TRUE){
    
    #random_points dataframe already exists 
    random_points = rbind(random_points,c("EASTING"=easting_random,
                                          "NORTHING"=northing_random))
  }else{
    #random_points dataframe doesn't exists yet...
    random_points = data.frame("EASTING"=easting_random, "NORTHING" = northing_random)
  }
  
  random_point_count=random_point_count + 1
}

write.csv(random_points, fout, row.names = FALSE)

#Plotting#
library(ggplot2)
library(reshape2)

#Output file name
plot_out = "random_points_CELL_CENTRES_2.png"

#Convert raster to matrix

mat2 = as.matrix(rast2)
dat2 = melt(mat2)
names(dat2) = c("easting","northing","value")

#Shift over the easting and northing values by half pixel (ggplot uses coordinates as midpoints, whereas we want it to plot from the mid-point minus half a pixel lenght or width)

dat2$easting = dat2$easting - 0.5
dat2$northing = dat2$northing - 0.5

p = ggplot(dat2, aes(x = easting, y = northing)) +
  geom_tile(aes(fill=value), colour="grey20") +
  scale_fill_gradientn(colours = terrain.colors(10)) +
  geom_point(data=random_points, mapping = aes(x=EASTING, y=NORTHING), colour="black") +
  labs(x = "Easting", y = "Northing") #+ 
#theme(legend.position = "none")
#p

ggsave(plot_out, plot = p, height = 6, width = 7, dpi = 500)













#For Data functions# 
#This 
library(Matrix)
library(spdep) #install.packages("spdep") 
#
#
r = matrix(sample(1:9,100, replace = TRUE),10,10)
#
Aggregate = function(x,aggLevel, binaryAgg = FALSE){

  if(aggLevel <= 1) return(x)
  
  matSize = length(x)^0.5
  m = (ceiling(matSize/aggLevel)-1)
  Xnew = rep(0, ceiling(matSize/aggLevel)^2)
  
  for (j in 0:m){                        # Row Position (aggregated matrix)
    for (i in 1:(m+1)){                # Column Position (aggregated matrix)
      for(dj in 0:(aggLevel-1)){     # Row Position (original matrix)
        for(di in 1:aggLevel){  # Column Position (original matrix)
          Xnew[(j*(m+1)) + i] = Xnew[(j*(m+1)) + i] + X[((j*aggLevel)+dj)*matSize + ((aggLevel*(i-1))+di)]
        }
      }         
      Xnew[(j*(m+1)) + i] = Xnew[(j*(m+1)) + i] / (aggLevel^2)
      if(binaryAgg){if(Xnew[(j*(m+1)) + i] > 0.5){Xnew[(j*(m+1)) + i] = 1} else {Xnew[(j*(m+1)) + i] = 0}}
    }
  }
  
  return (Xnew)
  
}

#Review Spatial Panel Model#---- install.packages("spml")
#Spatial Panel Model (by Maximum Likelihood)#















































#

#Simple examples of grid
#1
#r =raster(ncol = 30, nrow = 20) #Raster layer to test
#r[] = 1:(30*20) #giving values to cells 
#plot(r)
#r[runif(30*20) >= 0.30] = NA #Randomly *unselect* 70% of data
#plot(r)

#2
#r = raster(ncol=5,nrow =5)
#r[] = 1:(5*5)
#plot(r)