////////////////////////////////////////////////////////////////
// Hi there! This is the official README for my project! ///////
// I'm Courtney Stevens, author and manufacturer 	 ///////
//    This project is licensed under a CC-BY Creative    ///////
//  Commons Attribution 4.0 International License. Users ///////
//    may distribute, remix, adapt, and build upon the 	 ///////
//    material in any medium so long as attribution is   ///////
// 	given to the creator (which is me, Courtney). 	 ///////
////////////////////////////////////////////////////////////////

// My project, "Modelling Species Distribution for 4 Species of Maple (Acer): Climate adaptation scenarios for a resilient future at UBC Botanical Gardens" was done entirely in R and final maps were created using QGIS. 

// I started this project end of August of 2021 and finished April of 2022 as part of my Master's degree at UBC: Masters in 
Geomatics for Environmental Management. 

// When opening and using the CStevens_Maxent.R script, having the working directory set to this file should enable you to run the script flawlessly. 

// Questions and grievances can be sent to :
	ms.courtneystevens@gmail.com

	//////////////////////////////////////////////////
	////////          ///////////////   ((    ,///////
	/////      /#((      /////////     #%###      ////
	///. .   /###((%#      /////       #%###       ///
	///  (##(. #%#(( ##%##  ///         #%##        //
	// ,#####%(#%##(##%%#   ///   (###(( . (##%###( //
	///     ####%%%%%%     ////  (#%#####%##%##**   //
	////  /((###  ##%%%## //////.    ./  (         ///
	//////              //////////               /////
	//////////      ///////////////////     ./////////
	//////////////////////////////////////////////////
	///////((  (((/    ,///////////       ,      /////
	//// . ######## (((   ///////       *          ///
	///  /(###########(    /////.   &/%(/    ** (&.*//
	//, .((##%#########((   ///    (*** ,/ ((/(  /( ,/
	//*  ((#%#%#######(     /// /#   ,*, // (   .(/  /
	///       , (((###     /////     *,  (*     (/  //
	////*                 ///////   *    (,     (  ///
	///////            /////////////            //////
	//////////////////////////////////////////////////

>>>^ This was meant to be 4 leaves of the 4 species I modelled

List of Files:

2.5m-elevation [folder]
	This is just the downloaded elevation file from WorldClim. https://www.worldclim.org/data/worldclim21.html - This elevation was used in creating WorldClim 2.1, derived from SRTM elevation data. Theres a simple tif file in this folder at 2.5 minute resolution. 

climatedata [folder]
	This contains ALL climate data, present and future, used in the study. Within this folder there are more folders that, guess what, contain more folders. These are all essentially downloaded from WorldClim at 2.5 minute resolution: https://www.worldclim.org/data/cmip6/cmip6_clim2.5m.html

maps_final [folder]
	Contains 6 final maps (as .png - not shapefiles) made with QGIS: four are the results from SDMs; insitu_China is a simple map showing the points I was given in China; UBCBG is the map used to showcase my study area. 

prediction_rasters_Acer [folder]
	This is the creme-de-la-creme of my project. There are approximately 20 tif files of rasters from predicting probable suitability of where these 4 species of Acer will exist in the future. Whats that? All the SSP_ files look identical? Yeah, I'm aware. 

shapefiles [folder]
	A conglomerate of relevant shapefiles used (and not used) in this project. UBCBG shapefile is in here. A 'Natural_Earth_quick_start' has a handful of different shapefiles for global parameters - country boundaries, prominent cities, roadways, etc. Suitable to get started in map-making within QGIS. 

Acer_griseum_both.csv
	UBCBG-provided points of one of the endangered species. I believe this could have also been downloaded from the IUCN Red List: https://www.iucnredlist.org/species/193593/2244567

Acer_pentaphyllum_both.csv
	UBCBG-provided points of the other endangered species. Also could have been downloaded from the IUCN Red List: https://www.iucnredlist.org/species/193850/2285958

all_model_results.xlsx
	Excel file containing results from running Maxent appropriately - each sheet is for a different species, with the first sheet representing the average bioclimatic variable contribution and final AUC. 

CStevens_Maxent.R
	My final script, my pride, my glory. Heavily annotated - like every step is annotated - you're welcome. You should have full ability to run this script so long as you set the working directory to be this folder and have all the applicable libraries downloaded. 

raster_metadata.xml
	The related metadata for all of the rasters located in 'prediction_rasters_Acer' generated with ArcGIS Pro. 

README
	What you're reading now

shps_metadata.xml
	The related metadata for all point data including Acer_griseum_both, Acer_pentaphyllum_both, and acircpts and amacpts located in the 'shapefiles' folder; generated with ArcGIS Pro

ubcbg_treepts.csv
	Cleaned dataset of the tree inventory that UBCBG sent me, containing georeferenced locations of the trees I was interested in mapping. 

UBCBotanicalGardenProperty.xml
	The related metadata for the polygon parameter of UBC Botanical Garden. Generated with ArcGIS Pro. 