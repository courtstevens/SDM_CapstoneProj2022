#############################################
###### CStevens Maxent for SDM Script #######
# Adapted from Drew Kerkhoff's SDM tutorial #
#https://cmerow.github.io/RDataScience/3_6_Teaching_Ecoinformatics.html

#for installing any of the following packages, follow the format:
#install.packages("dismo")
#then run library(dismo) as shown below:

library(dismo)#species distribution modelling
#dismo also loads raster and sp
library(dplyr) #general data manipulation
library(terra) #for reading tif files & raster manipulation
library(maps) #draw geographical maps
library(mapdata) #extra maps
library(maptools) #handling spatial objects
library(rJava) #to use MaxEnt
library(sf) #reading in garden polygon
library(RColorBrewer) #for nice colours in maps
library(corrplot) #for looking at correlation
library(ggplot2) #for final histograms


##Data Loading & Processing:
setwd("C:/Users/csteve02.stu/Desktop/CStevens_FCOR599_AcerSDM_20220408")

#species data for A. griseum and A. pentaphyllum ... & any other .csv data
tree_data = read.csv(file = "Acer_griseum_both.csv")
tree_data = read.csv(file = "Acer_pentaphyllum_both.csv")

tree_data = dplyr::select(tree_data, #selecting location data only
                   "long" = decimalLongitude,
                   "lat" = decimalLatitude)
#removing duplicate points from dataset
tree_data_dups = duplicated(tree_data[, c("long", "lat")]) #isolate duplicates
tree_data <- tree_data[!tree_data_dups,] #remove duplicates from data



## using GBIF to get A. macrophyllum and A. circinatum data
# can use this for any other species not considered scarce/endangered
tree_data_raw <- gbif("Acer macrophyllum") #ran it with AM first
tree_data_raw <- gbif("Acer circinatum") #repeat with AC

#Cleaning GBIF data
tree_data <- subset(tree_data_raw, (!is.na(lat)) & (!is.na(lon))) #remove NA values
tree_data_dups = duplicated(tree_data[, c("lon", "lat")]) #isolate duplicates
tree_data <- tree_data[!tree_data_dups,] #remove duplicates from data
#removing trees considered 'preserved specimen' - only care about living records
tree_data <- subset(tree_data, 
                    basisOfRecord != "PRESERVED_SPECIMEN")
tree_data <- dplyr::select(tree_data, #we only care about location
                    "long" = lon,
                    "lat" = lat)
plot(tree_data)
#only want species that are on western coast (for my project)
tree_data <- tree_data[which(tree_data$long > -130 & tree_data$long < -90), ]


#elevation data ... just took this from WorldClim SRTM
elevation <- rast("2.5m-elevation/wc2.1_2.5m_elev.tif")
names(elevation) <- c("elev")

#present environmental data for extraction
present_climate = getData(name = "worldclim", #downloading climate data 'standard' b/w 1970-2000
                          var = "bio", #all 19 bioclimatic variables; or tmax, precip...
                          res = "2.5", #resolution at 2.5 minutes
                          path = "climatedata") #where data is stored

#future environmental data for prediction
#downloaded each from WorldClim's future projections
#https://www.worldclim.org/data/cmip6/cmip6climate.html
#reading SSP245 for 2040-2060
r245_50 <- rast("climatedata/ssp245_2041-2060/spatial03/worldclim/cmip6/7_fut/_/2.5m/MIROC-ES2L/ssp245/wc2.1_2.5m_bioc_MIROC-ES2L_ssp245_2041-2060.tif") #rasterize tif file using terra::rast

#reading SSP245 for 2080-2100
r245_90 <- rast("climatedata/ssp245_2081-2100/spatial03/worldclim/cmip6/7_fut/_/2.5m/MIROC-ES2L/ssp245/wc2.1_2.5m_bioc_MIROC-ES2L_ssp245_2081-2100.tif")

#reading SSP370 for 2040-2060
r370_50 <- rast("climatedata/ssp370_2041-2060/spatial03/worldclim/cmip6/7_fut/_/2.5m/MIROC-ES2L/ssp370/wc2.1_2.5m_bioc_MIROC-ES2L_ssp370_2041-2060.tif")

#reading SSP370 for 2080-2100
r370_90 <- rast("climatedata/ssp370_2081-2100/spatial03/worldclim/cmip6/7_fut/_/2.5m/MIROC-ES2L/ssp370/wc2.1_2.5m_bioc_MIROC-ES2L_ssp370_2081-2100.tif")

#setting list with bioclimatic variable names for cleaning future environment layers
bio_list <- c("elev", "bio1","bio2","bio3","bio4","bio5","bio6","bio7","bio8","bio9","bio10",
              "bio11","bio12","bio13","bio14","bio15","bio16","bio17","bio18","bio19")

#Creating a function to combine elevation raster with climate rasters & make them accessible
elev_fun <- function (elev, myrast){
  elev <- stack(elev) #stack elevation
  myrast <- stack(myrast) #stack raster
  elev <- resample(elev, myrast) #resample elevation so it's the same paramaters as rasters
  fullrast <- stack(elev, myrast) #stack resampled elevation with raster
  names(fullrast) <- bio_list #clean the names
  return(fullrast) #return raster
}

present_climate <- elev_fun(elevation, present_climate)
climate_ssp245_2050 <- elev_fun(elevation, r245_50) #SSP245 for 2040-2060
climate_ssp245_2090 <- elev_fun(elevation, r245_90) #SSP245 for 2080-2100
climate_ssp370_2050 <- elev_fun(elevation, r370_50) #SSP370 for 2040-2060
climate_ssp370_2090 <- elev_fun(elevation, r370_90) #SSP370 for 2080-2100

#examining the data:
data("wrld_simpl")
plot(wrld_simpl, 
     xlim=c(min(tree_data$long)-3,
            max(tree_data$long)+3), 
     ylim=c(min(tree_data$lat)-3,
            max(tree_data$lat)+3), 
     axes=TRUE)
points(tree_data$long, tree_data$lat, pch = 2, col = "darkgreen")

#establish extent of study area
#Want a buffer of the min & max values of points so background points can be generated properly
map_extent <- extent(min(tree_data$long) - round((max(tree_data$long) - min(tree_data$long))/2),
                     max(tree_data$long) + round((max(tree_data$long) - min(tree_data$long))/2),
                     min(tree_data$lat) - round((max(tree_data$lat) - min(tree_data$lat))/2), 
                     max(tree_data$lat) + round((max(tree_data$lat) - min(tree_data$lat))/2))
present_climate_crop <- crop(present_climate, extent(map_extent)) #crop climate vars to extent

#plotting to check it worked
plot(present_climate_crop$bio11, main = "A. pentaphyllum distribution [bio11]")
map('worldHires',
    xlim=c(min(tree_data$long)-30, 
           max(tree_data$long)+30), 
    ylim=c(min(tree_data$lat)-30,
           max(tree_data$lat)+30),
    fill = FALSE,
    add = TRUE)
points(tree_data$long, tree_data$lat, pch = 2, col = "darkgreen")


##Species Distribution Modelling Time

#establish empty dataframe for storing information
treevar <- data.frame(row.names = bio_list) %>% #same number of rows as rasters' layers
  tibble::rownames_to_column(var = "bv") %>% #make biovariable a column for sorting
  arrange(bv) #sort alphabetically so the following loop can be integrated seamlessly

fin_auc <- c() #empty list for final AUC values

#for loop to run Maxent bootstrap method - DO NOT USE FOR GBIF DATA
for (i in 1:nrow(tree_data)) {
  if (i == nrow(tree_data)){ #break the loop after 9 iterations
    treevar <- tibble::column_to_rownames(treevar, "bv...1") #return biovariables to row names
    treevar <- select(treevar, -contains("bv")) #remove biovariable columns
    treevar <- treevar %>% 
      mutate(avgcont = rowMeans(treevar)) %>% #average row values for final assessment
      arrange(avgcont) #sort by least to most significant variables
  } else{
  treetest <- tree_data[i,] #isolate just one row from data
  treetrain <- tree_data[-i,] #isolate all other values for training data
  tree.me <- maxent(present_climate_crop, treetrain) #MaxEnt model using training data
  me.results <- data.frame("varc" = plot(tree.me)) %>% #create dataframe with maxent results
    tibble::rownames_to_column(var = "bv") %>% #transfer row names to column; bv = 'biovariable'
    arrange(bv) #sort names of dataframe so binding can happen appropriately
  treevar <- bind_cols(treevar, me.results) #bind variables to table
  #evaluate the model with the test
  bg <- randomPoints(present_climate_crop, #generate background points for evaluating MaxEnt 
                     n = nrow(tree_data)) #background points are same number as presence points
  model.eval <- evaluate(tree.me, #evaluate model
                         p = treetest, #using reserved testing presence points
                         a = bg, #and generated background 'absence' points
                         x = present_climate) #present climate global
  fin_auc <- bind_cols(fin_auc, model.eval@auc) #create table with AUC values
  i <- i + 1
  }
}

#save the variables! :) 
write.csv(treevar, "C:/Users/csteve02.stu/Desktop/FCOR 599/maxent/ag_maxent_treevariable.csv")
write.csv(fin_auc, "C:/Users/csteve02.stu/Desktop/FCOR 599/maxent/ag_auc.csv")

avg_auc <- rowMeans(fin_auc) #average AUC values from model iterations
avg_auc #average AUC to report
plot(model.eval, "ROC") #plotting the last AUC curve

#dropping the variables that had 0 variable contribution to the model
#variables that had 0% contribution for 4 or more runs were excluded for AP / AG
#variables that contributed more than 5% to the model were included for AM
#variables contributed more than 3% on average and did not collineate for AC

droplist <- c("bio1","bio3","bio4","bio5","bio7","bio8","bio10",
              "bio12","bio13","bio15", "bio16","bio17","bio18","bio19") 
#create list of var to drop

present_climate_croptrim <- dropLayer(present_climate_crop, droplist) #drop the layers
present_climatetrim <- dropLayer(present_climate, droplist)

tree.mefinal <- maxent(present_climate_croptrim, treetrain) #final model w dropped layers

#looking at correlations... see if any variables correlate.
#If variables correlate, the model may be overfitting
tree_corr <- extract(present_climate_croptrim, tree_data)
tree_corr <- cor(tree_corr)
corrplot.mixed(tree_corr)

#Predict current distribution - TAKES AWHILE, GRAB TEA
tree.predict <- dismo::predict(tree.mefinal, #predict with final model 
                               present_climatetrim, #give it present climate data
                               progress = "text") #progress bar in console

x11() #open new window to view plot
plot(tree.predict, main = "Present Climate") #plotting present climate global
plot(tree.predict, ext = extent(-124,-122,48,50), main = "UBCBG Extent") #present climate to UBCBG
writeRaster(tree.predict, "predrasters/ag_current.tif") #save as tif just for reference


#Predicting the future... 
#cleaning up future climate data sets

climate_ssp245_2050trim <- dropLayer(climate_ssp245_2050, droplist)
climate_ssp245_2090trim <- dropLayer(climate_ssp245_2090, droplist)

climate_ssp370_2050trim <- dropLayer(climate_ssp370_2050, droplist)
climate_ssp370_2090trim <- dropLayer(climate_ssp370_2090, droplist)


#Predicting the future... GRAB SOME TEA THIS TAKES AWHILE
tree.ssp245.50 <- dismo::predict(tree.mefinal, #using dismo::predict function, give it the model
                                 climate_ssp245_2050trim, #give it the raster file to apply predictions to 
                                 progress = "text") #progress bar so we know how long it takes
writeRaster(tree.ssp245.50, "predrasters/ag_ssp245_50.tif") #save as raster ... MIND THE FILE NAME

tree.ssp245.90 <- dismo::predict(tree.mefinal, climate_ssp245_2090trim, progress = "text")
writeRaster(tree.ssp245.90, "predrasters/ag_ssp245_90.tif")

tree.ssp370.50 <- dismo::predict(tree.mefinal, climate_ssp370_2050trim, progress = "text")
writeRaster(tree.ssp370.50, "predrasters/ag_ssp370_50.tif")

tree.ssp370.90 <- dismo::predict(tree.mefinal, climate_ssp370_2090trim, progress = "text")
writeRaster(tree.ssp370.90, "predrasters/ag_ssp370_90.tif")

#In case the environment has been erased and you don't want to re-run everything:
#climate_ssp245_2050 <- rast("C:/Users/csteve02.stu/Desktop/FCOR 599/maxent/predrasters/ag_ssp245_50.tif")
#tree.ssp245.50 <- stack(climate_ssp245_2050)

#plotting to see if our future predictions worked
plot(tree.ssp245.50, ext = extent(-124.5,-123,49,49.5))
map('worldHires',
    xlim=c(min(-124), 
           max(-122)), 
    ylim=c(min(48),
           max(50)),
    fill = FALSE,
    add = TRUE)

points(tree_data$long, tree_data$lat, pch = 2, col = "darkgreen")


#Looking at the differences between present and future
pres50_245diff <- tree.ssp245.50 - tree.predict #future raster minus cells from 'present' raster
plot(pres50_245diff,
     main = "2050/245 to Present Difference",
     ext = extent(-124,-122,48,50)) #limiting to UBCBG extent

pres50_370diff <- tree.ssp370.50 - tree.predict
plot(pres50_370diff,
     main = "2050/370 to Present Difference",
     ext = extent(-124,-122,48,50))

pres90_245diff <- tree.ssp245.90 - tree.predict
plot(pres90_245diff,
     main = "2090/245 to Present Difference",
     ext = extent(-124,-122,48,50))

pres90_370diff <- tree.ssp370.90 - tree.predict
plot(pres90_370diff,
     main = "2090/370 to Present Difference",
     ext = extent(-124,-122,48,50))

#Reading in this shapefile that UBCBG gave me... 
ubcbg <- st_read("shapefiles/UBCBotanicalGardenProperty_M.shp")
ubcbg <- ubcbg$geometry #limit to just geography of shapefile

plot(tree.ssp245.50, #plotting the future
     ext = extent(-124,-122,49,49.5),
     col = brewer.pal(n=9, name="OrRd"))
plot(ubcbg, add = TRUE, border = "black") #adding UBCBG shape to the plot

plot(tree.ssp245.90,
     ext = extent(-124,-122,49,49.5),
     col = brewer.pal(n=9, name="OrRd"))
plot(ubcbg, add = TRUE, border = "black")
  
plot(tree.ssp370.50,
     ext = extent(-124,-122,49,49.5),
     col = brewer.pal(n=9, name="OrRd"))
plot(ubcbg, add = TRUE, border = "black")

plot(tree.ssp370.90,
     ext = extent(-124,-122,49,49.5),
     col = brewer.pal(n=9, name="OrRd"))
plot(ubcbg, add = TRUE, border = "black")


#creating plots looking at UBCBG for different variables
treevar <- data.frame(extract(present_climate, tree_data)) #extract variables from present climate
ubcbgcoord <- data.frame(st_coordinates(ubcbg)) %>% #extract coordinates from ubcbg polygon
  summarise(across(everything(), mean)) %>% #finding approximate coordinate of ubcbg
  select("long" = X,
         "lat" = Y)
ubcbgvar <- data.frame(extract(present_climate, ubcbgcoord)) #extract variables for ubcbg from present climate

#extract variables from the future climates too ... 
ubcbgvar24550 <- data.frame(extract(climate_ssp245_2050, ubcbgcoord))
ubcbgvar37050 <- data.frame(extract(climate_ssp370_2050, ubcbgcoord))
ubcbgvar24590 <- data.frame(extract(climate_ssp245_2090, ubcbgcoord))
ubcbgvar37090 <- data.frame(extract(climate_ssp370_2090, ubcbgcoord))

#creating plot of the treevariable curve with UBCBG markers - visualizing limit of climate in Garden
ggplot(treevar, aes(x=bio1)) + 
  geom_density() + 
  geom_vline(xintercept = ubcbgvar$bio1/10, col = "darkgreen") +
  geom_vline(xintercept = ubcbgvar24590$bio1, col = "blue") +
  geom_vline(xintercept = ubcbgvar37090$bio1, col = "red") +
  labs(x = "Average Temperature of the Driest Quarter (Celsius)", y = "Frequency")


##Isolating points from UBCBG tree inventory data
trees <- read.csv("C:/Users/csteve02.stu/Desktop/FCOR 599/GIS/Acer_ItemData_Oct2021.csv")

#select the columns that we care about
trees <- select(trees,
                "lat" = LocationCoordLatDD,
                "long" = LocationCoordLongDD,
                "species" = Species)
#filter the rows for species that we care about
treesfin <- filter(trees, species == "griseum" | species == "pentaphyllum" |
                   species == "macrophyllum" | species == "circinatum")

#removing duplicates
dups = duplicated(treesfin[, c("lat", "long")])
treesfin <- treesfin[!dups,]

#save as CSV
write.csv(treesfin, "C:/Users/csteve02.stu/Desktop/FCOR 599/maxent/ubcbgpts.csv")