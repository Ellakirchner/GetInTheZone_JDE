// Script to extract spatial data from Google Earth Engine. 
// Instructions: 
// Section 1 is extracting the agri-environmental variables at the crop cut locations (hereafter referred to as "oaf points"). For data privacy reasons, the dataset with the locations is not published as part of this replication package, but any other dataset with up to ca. 20k locations that is saved in your personal GEE Account under the name "oaf_points" can be used to run the script.  
// Section 2 can be replicated without the original dataset as it generates 50k random locations on maize crop land (section 2.0) and extracts the corresponding agri-environmental inputs (sections 2.1-2.5).
// To replicate section 2, please follow these steps:
// 1) Import the shapefiles for admin 1 level data in Kenya as a GEE asset to your GEE account.
// 2) Replace the words "your_username" with your actual username. 
// 3) Run scripts one by one. The beginning of a new script is marked by two lines of "%" that enclose the script title which are numbered in the format "1.1".
// Attention has to be paid to the indicated timeframes in each section. Some sections need to be run several times on different years to extract the datasets for the entire relevant period (e.g., for CHIRPS data only two years of data can be extracted at a time).

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// 1.1 - Extract CHIRPS for oaf data points
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

// Purpose: Extract and export precipitation data for crop cut locations in study region in Kenya 
// Input files: CHIRPS (pentad), crop cut locations (GEE asset) 
// Export files: CSV file of precipitation data for four years at a time (from 2000-2020)

// 1. Import necessary files:
// - geopoints of cropcuts: import the asset, 
// - chirps data
// - polygon of region of interest
var roi_reg = ee.Geometry.Polygon(
        [[[33.57858523739146, 1.627217257564576],
          [33.57858523739146, -1.5803713344511332],
          [38.21481570614146, -1.5803713344511332],
          [38.21481570614146, 1.627217257564576]]], null, false);   
var points = ee.FeatureCollection("users/your_username/oaf_points");

//Import CHIRPS and select precipitation band: 
var chirpdata_pentad = ee.ImageCollection('UCSB-CHG/CHIRPS/PENTAD').select('precipitation');

// 2. Visualize imports: 
Map.addLayer(points,{},'Crop cut location' ); 
Map.addLayer(roi_reg,{}, 'ROI'); 
Map.centerObject(points.geometry(), 8);
Map.setOptions("HYBRID");

// 3. Define start and endyear: 
// - for this dataset, the timeframe needs to be limited to 2 years (otherwise the export dataset is too big)
// --> Decide on start and end year, uncomment export lines at end of the script, run task, repeat with new year of interest
var startyear = 2016;
var endyear = 2020;

// 4. Filter for the timeframe of interest:
var chirpday= chirpdata_pentad.filterDate(startyear + '-05-01', endyear + '-09-01'); 
var chirpdaily = chirpday.map(function(img){
  return img.copyProperties(img,['system:time_start','system:time_end']);
});

// 5. Stack all the images together:
var chirpstack = ee.Image(chirpdaily.iterate(function(imgIn, imgOut){
  var dateString = ee.Date(imgIn.get('system:time_start')).format('yyyy-MM-dd');
  return ee.Image(imgOut).addBands(imgIn.select(['precipitation']).clip(roi_reg)
          .rename(dateString));
},ee.Image()));
print('CHIRPS collection', chirpstack);

// 6. Extract the Precipitation data:
var chirps_done = chirpstack.reduceRegions({
  collection: points,
  reducer: ee.Reducer.first(),
  scale: 250,});
print("This works:", chirps_done.limit(1));

// 7. Export the data to csv.
Export.table.toDrive({
  collection: chirps_done,
  description: 'chirpstack_'+ startyear + '_' + endyear+ '_oaf',
  folder: 'gee_oaf',
  fileFormat: 'CSV'
});

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// 1.2 - Extract EVI for oaf data points
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// Purpose: Extract and export EVI data for crop cut locations in study region in Kenya 
// Input files: MODIS Aqua, crop cut locations (GEE asset) 
// Export files: Two CSV files of EVI data, each for 10 years (2000-2010, 2010-2020)
// Comment: The script is organized in two sections (marked by %%%-lines)
// The first section exports the data for the years 2000-2010. 
// The second section does the exact same operation for the years 2011-2020.

// 1. Import necessary datasets 
// - geopoints of cropcuts: import the asset, 
// - polygon of region of interest
var roi_reg = ee.Geometry.Polygon(
        [[[33.57858523739146, 1.627217257564576],
          [33.57858523739146, -1.5803713344511332],
          [38.21481570614146, -1.5803713344511332],
          [38.21481570614146, 1.627217257564576]]], null, false);
Map.centerObject(roi_reg, 8);
var points = ee.FeatureCollection("users/your_username/oaf_points");

// 2. Visualize inputs on Map: 
// Choose one image:
var displayMODIS = ee.ImageCollection('MODIS/061/MYD13Q1').filterDate('2010-08-01', '2010-08-30').select('EVI').first();
// MODIS EVI is scaled by 0.0001, so we need to rescale it
var eviRescaled = displayMODIS.multiply(0.0001).clip(roi_reg);  // use your ROI
// Visualization parameters (white and yellow - bare soil or sparse vegetation; green to darkgreen - higher vegetation density)
var visParams = {
  min: 0,
  max: 1,
  palette: ['white', 'yellow', 'green', 'darkgreen']
};
// Display on the map
Map.addLayer(eviRescaled, visParams, 'EVI MODIS AQUA (Aug 2010)');
Map.addLayer(roi_reg, {}, 'ROI');
Map.addLayer(points, {}, 'Crop cut locations');
Map.centerObject(points.geometry(), 8);
Map.setOptions("HYBRID");

// 3. Select EVI from general MODIS collection:   
var modis_aqua = ee.ImageCollection('MODIS/061/MYD13Q1').select('EVI');

// 4. Filter data in two steps (2000-2010 first and in the later section 2011-2020) to place in two csvs 
var startyear = 2000;
var endyear = 2010;
var MODISeviday= modis_aqua.filterDate(startyear + '-04-01', endyear + '-09-25');

// 5. Extract the imagery date from each image
var modEVIday = MODISeviday.map(function(img){
  return img.copyProperties(img,['system:time_start','system:time_end']);
});

// 6. Stack all the images together
var evistack = ee.Image(modEVIday.iterate(function(imgIn, imgOut){
  var dateString = ee.Date(imgIn.get('system:time_start')).format('yyyy-MM-dd');
  return ee.Image(imgOut).addBands(imgIn.select(['EVI']).clip(roi_reg)
          .rename(dateString));
},ee.Image()));

print('EVI stack', evistack);

// 7. Extract the EVI points for every crop cut location
var EVI = evistack.reduceRegions({
  collection: points,
  reducer: ee.Reducer.first(),
  scale: 250,});
  
print("This works:", EVI.limit(1));

// Uncomment this block if you want to export the data: 
// 8. Export the data to csv.
Export.table.toDrive({
  collection: EVI,
  description: 'evistack_oaf_'+ startyear + '_' + endyear,
  folder: 'gee_oaf',
  fileFormat: 'CSV'
});

// Repeat for different years: 

// 1. Filter data in two steps (now 2011-2020) to place in two csvs (one is too large)
var startyear = 2011;
var endyear = 2020;
var MODISeviday= modis_aqua.filterDate(startyear + '-04-01', endyear + '-09-25');

// 2. Extract the imagery date from each image
var modEVIday = MODISeviday.map(function(img){
  return img.copyProperties(img,['system:time_start','system:time_end']);
});

// 3. Stack all the images together
var evistack = ee.Image(modEVIday.iterate(function(imgIn, imgOut){
  var dateString = ee.Date(imgIn.get('system:time_start')).format('yyyy-MM-dd');
  return ee.Image(imgOut).addBands(imgIn.select(['EVI']).clip(roi_reg)
          .rename(dateString));
},ee.Image()));
print('EVI stack (2nd period)', evistack);

// 4. Extract the EVI for every crop cut location
var EVI = evistack.reduceRegions({
  collection: points,
  reducer: ee.Reducer.first(),
  scale: 250,});
print("This works:", EVI.limit(1));

// Uncomment this block for export. 
// 5. Export the data to csv.
Export.table.toDrive({
  collection: EVI,
  description: 'evistack_oaf_'+ startyear + '_' + endyear,
  folder: 'gee_oaf',
  fileFormat: 'CSV'
});

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// 1.3 - Extract NDVI for oaf data points
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// Purpose: Extract and export NDVI data for crop cut locations in study region in Kenya 
// Input files: MODIS Aqua, crop cut locations (GEE asset) 
// Export files: Two CSV files of NDVI data, each for 10 years (2000-2010, 2010-2020)
// Comment: The script is organized in two sections (marked by %%%-lines)
// The first section exports the data for the years 2000-2010. 
// The second section does the exact same operation for the years 2011-2020.

// 1. Import necessary datasets 
// - geopoints of cropcuts: import the asset, 
// - polygon of region of interest
var roi_reg = ee.Geometry.Polygon(
        [[[33.57858523739146, 1.627217257564576],
          [33.57858523739146, -1.5803713344511332],
          [38.21481570614146, -1.5803713344511332],
          [38.21481570614146, 1.627217257564576]]], null, false);
Map.centerObject(roi_reg, 8);
var points = ee.FeatureCollection("users/your_username/oaf_points");

// 2. Visualize inputs on Map: 
// Choose one image:
var displayMODIS = ee.ImageCollection("MODIS/061/MYD13Q1").filterDate('2010-08-01', '2010-08-30').select('NDVI').first();
// MODIS NDVI is scaled by 0.0001, so we need to rescale it
var ndviRescaled = displayMODIS.multiply(0.0001).clip(roi_reg);  // use your ROI
// Visualization parameters (white and yellow - bare soil or sparse vegetation; green to darkgreen - higher vegetation density)
var visParams = {
  min: 0,
  max: 1,
  palette: ['white', 'yellow', 'green', 'darkgreen']
};
// Display on the map
Map.addLayer(ndviRescaled, visParams, 'NDVI MODIS AQUA (Aug 2010)');
Map.addLayer(roi_reg, {}, 'ROI');
Map.addLayer(points, {}, 'Crop cut locations');
Map.centerObject(points.geometry(), 8);
Map.setOptions("HYBRID");

// 3. Select NDVI from general MODIS collection:   
var modis_aqua = ee.ImageCollection('MODIS/061/MYD13Q1').select('NDVI');

// 4. Filter data in two steps (2000-2010 first and in the later section 2011-2020) to place in two csvs 
var startyear = 2000;
var endyear = 2010;
var MODISndviday= modis_aqua.filterDate(startyear + '-04-01', endyear + '-09-25');

// 5. Extract the imagery date from each image
var modNDVIday = MODISndviday.map(function(img){
  return img.copyProperties(img,['system:time_start','system:time_end']);
});

// 6. Stack all the images together
var ndvistack = ee.Image(modNDVIday.iterate(function(imgIn, imgOut){
  var dateString = ee.Date(imgIn.get('system:time_start')).format('yyyy-MM-dd');
  return ee.Image(imgOut).addBands(imgIn.select(['NDVI']).clip(roi_reg)
          .rename(dateString));
},ee.Image()));
print('NDVI stack', ndvistack);

// 7. Extract the NDVI points for every crop cut location
var NDVI = ndvistack.reduceRegions({
  collection: points,
  reducer: ee.Reducer.first(),
  scale: 250,});
print("This works:", NDVI.limit(1));

//Uncomment this block if you want to export the data: 
// 8. Export the data to csv.
Export.table.toDrive({
  collection: NDVI,
  description: 'ndvistack_oaf_' + startyear + '_' + endyear,
  folder: 'gee_oaf',
  fileFormat: 'CSV'
});

// Repeat for different years: 

// 1. Filter data in two steps (now 2011-2020) to place in two csvs (one is too large)
var startyear = 2011;
var endyear = 2020;
var MODISndviday= modis_aqua.filterDate(startyear + '-04-01', endyear + '-09-25');

// 2. Extract the imagery date from each image
var modNDVIday = MODISndviday.map(function(img){
  return img.copyProperties(img,['system:time_start','system:time_end']);
});

// 3. Stack all the images together
var ndvistack = ee.Image(modNDVIday.iterate(function(imgIn, imgOut){
  var dateString = ee.Date(imgIn.get('system:time_start')).format('yyyy-MM-dd');
  return ee.Image(imgOut).addBands(imgIn.select(['NDVI']).clip(roi_reg)
          .rename(dateString));
},ee.Image()));
print('NDVI stack (2nd period)', ndvistack);

// 4. Extract the NDVI for every crop cut location
var NDVI = ndvistack.reduceRegions({
  collection: points,
  reducer: ee.Reducer.first(),
  scale: 250,});
print("This works:", NDVI.limit(1));

// Uncomment this block for export. 
// 5. Export the data to csv.
Export.table.toDrive({
  collection: NDVI,
  description: 'ndvistack_oaf_' + startyear + '_' + endyear,
  folder: 'gee_oaf',
  fileFormat: 'CSV'
});

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// 1.4 - Extract Temp for oaf data points
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// Purpose: Extract and export temperature data for crop cut locations in study region in Kenya 
// Input files: ERA5 Daily, crop cut locations (GEE asset) 
// Export files: CSV file of temperature data from one season at a time (2000-2020)

// 1. Define start and endyear
// - for this dataset, the timeframe needs to be limited to one year (otherwise the export dataset is too big)
// --> Decide on year, uncomment export lines at end of the script, run task, repeat with new year of interest
var year=2000;

// 2. Import necessary files:
// - geopoints of crop cuts: import the asset, 
// - temp data: from ERA5 
// - polygon of region of interest
var roi_reg = ee.Geometry.Polygon(
        [[[33.57858523739146, 1.627217257564576],
          [33.57858523739146, -1.5803713344511332],
          [38.21481570614146, -1.5803713344511332],
          [38.21481570614146, 1.627217257564576]]], null, false); 
var points = ee.FeatureCollection("users/your_username/oaf_points");
var temp=ee.ImageCollection("ECMWF/ERA5_LAND/DAILY_AGGR").filterDate(year+'-05-01', year+'-09-01');
print(temp);

// 3. Compute GDD: 
function gdd(temp){
  var mintemp = temp.select('temperature_2m_min').subtract(273.15);
  var maxtemp = temp.select('temperature_2m_max').subtract(273.15);
  var meantemp= temp.select('temperature_2m').subtract(273.15);
  var doubletemp = mintemp.add(maxtemp);
  var meanDouble = doubletemp.divide(2);
  var meanDoubleAdj = meanDouble.where(meanDouble.gt(29), 29);
  var GDD = meanDoubleAdj.subtract(10); 
  var GDDNonNeg = GDD.where(GDD.lt(0), 0).rename('GDD');
  return GDDNonNeg.addBands(meantemp); 
  //return GDDNonNeg.addBands(mintemp).addBands(maxtemp); 
}
// Map the function over the collection
var gddtemp = temp.map(gdd);
print('gdd', gddtemp);

// 5. Extract the imagery date from each image
var tempday = gddtemp.map(function(img){
  return img.copyProperties(img,['system:time_start','system:time_end']);
});
print('tempday', tempday);

// 6. Stack all the images together
var tempstack = tempday.toBands().clip(roi_reg);

// 7. Print to console:
print('tempstack', tempstack);

// 8. Extract the temp points for every point included in the points dataset
var temp = tempstack.reduceRegions({
  collection: points,
 reducer: ee.Reducer.first(),
  scale: 250,});
print("This works:", temp.limit(5));

// 9. Export the data to csv.
Export.table.toDrive({
  collection: temp,
  description: 'gdd_'+year+'_oaf',
  folder: 'gee_oaf',
  fileFormat: 'CSV'
});

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// 1.5 - Extract Soil for oaf data points
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// Purpose: Extract and export soil data for crop cut locations in study region in Kenya 
// Input files: Different soil datasets saved as GEE assets, crop cut locations (also GEE asset) 
// Export files: CSV file of soil data, time-invariant data, hence only one output file


// 1. Import all the necessary datasets: 
// A. Soil data first: 
    var awcp30 = ee.Image ("users/kirchnere/AWC_1km"); 
        //based on Jin's file which was named: projects/jin-digitalag-lab/AfricaSoilGrids/af_agg_30cm_TAWCpF23mm__M_1km, but fresh download on 23-10-11 from https://files.isric.org/public/gyga/GYGA_Results/
    var erzd = ee.Image ("users/kirchnere/ERZD_1km");
        //based on Jin's file which was named: projects/jin-digitalag-lab/AfricaSoilGrids/af_ERZD_M_1km, but fresh download on 23-10-11 from https://files.isric.org/public/gyga/GYGA_Results/
    var al = ee.Image ("projects/jin-digitalag-lab/AfricaSoilGrids/af250m_nutrient_al_m_agg30cm");
    var soc = ee.Image ("projects/jin-digitalag-lab/AfricaSoilGrids/af_ORCDRC_T__M_sd1_250m");
    var totN = ee.Image ("projects/jin-digitalag-lab/AfricaSoilGrids/af250m_nutrient_n_m_agg30cm");
// B. Geopoints of crop cuts: 
    var points = ee.FeatureCollection("users/your_username/oaf_points");
  
// 2. Combine all the bands: 
var soil = ee.Image([])
  .addBands(erzd.rename('erzd'))
  .addBands(awcp30.rename('awcp30'))
  .addBands(al.rename('al'))
  .addBands(soc.rename('soc'))
  .addBands(totN.rename('totN'));
print(soil);

// 3. Reduce the datasets to the values of the geopoints: 
var soilpoints = soil.reduceRegions({
  collection: points,
 reducer: ee.Reducer.first(),
  scale: 250,});


//4. Export soil data:  
Export.table.toDrive({
  collection: soilpoints,
  description: 'soil_oaf',
  folder: 'gee_oaf',
  fileFormat: 'CSV'
});

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// 2.0 - Create 50k random points on maize cropland
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// Purpose of script: Generate 50k points on maize cropland in study region 
// Input files: GADM Admin1 level data for Kenya (saved as GEE asset)
// Export files: Feature Collection of the generated points exported as GEE asset

// Define ROI based on 14 counties in Kenya using GADM v4.1
var counties = [
  "Bungoma", "Busia", "Homa Bay", "Kakamega", "Kericho", "Kisumu", "Kisii",
  "Migori", "Nandi", "Nyamira", "Siaya", "Trans Nzoia", "Uasin Gishu", "Vihiga"
];
var gadm = ee.FeatureCollection("users/your_username/gadm41_KEN_1")
              .filter(ee.Filter.inList("NAME_1", counties));
var roi = gadm.geometry();

// Load the WorldCereal dataset
var dataset = ee.ImageCollection("ESA/WorldCereal/2021/MODELS/v100");
function mask_other(img) {
  return img.updateMask(img.neq(0))
}

// Apply the mask_other function to the collection
dataset = dataset.map(mask_other);

// Get the maize crop mask:
var maizeCollection = dataset
  .filter('product == "maize"');
  //.filter('season == "tc-maize-main"')

// Mosaic the filtered images to create a single image
var maizeImage = maizeCollection.mosaic().clip(roi);
print(maizeImage, 'Maize Image');

// Visualization specifics
var visualization_class = {
  bands: ["classification"],
  max: 100,
  palette: ["ff0000"]
};
Map.addLayer(maizeImage, visualization_class, 'Maize cropland');

// Filter image
var cropland = maizeImage.select('classification').eq(100).selfMask();

// Select point locations
var maizePoints = maizeImage.select('classification').stratifiedSample({
  numPoints: 50000,
  classBand: 'classification',
  classValues: [100],
  classPoints: [50000],
  region: cropland.geometry(),
  scale: 30,
  seed: 42,
  geometries: true
});

// Count and visualize:
var count = maizePoints.size();
Map.addLayer(maizePoints, {color: 'yellow'}, 'Cropland Points');
print(count, 'Number of Maize Points');

// Export Point feature: 
Export.table.toAsset({
  collection: maizePoints,               
  description: 'maize_points_asset_task', 
  assetId: 'users/your_username/maizePoints' 
});

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// 2.1 - Extract CHIRPS for 50k points
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// Purpose: Extract and export precipitation data for crop cut locations in study region in Kenya 
// Input files: CHIRPS (pentad), 50k geopoints randomly placed on maize cropland (GEE asset) 
// Export files: CSV file of precipitation data for two years at a time (2000-2020)

// 1. Import necessary files:
// - 50k random geopoints: import the asset, 
// - chirps data
// - polygon of region of interest
var roi_reg = ee.Geometry.Polygon(
        [[[33.57858523739146, 1.627217257564576],
          [33.57858523739146, -1.5803713344511332],
          [38.21481570614146, -1.5803713344511332],
          [38.21481570614146, 1.627217257564576]]], null, false);    
var points = ee.FeatureCollection("users/your_username/maizePoints");
//Import CHIRPS and select precipitation band: 
var chirpdata_pentad = ee.ImageCollection('UCSB-CHG/CHIRPS/PENTAD').select('precipitation');

// 2. Visualize imports: 
Map.addLayer(points,{},'Crop cut location' ); 
Map.addLayer(roi_reg,{}, 'ROI'); 
Map.centerObject(points.geometry(), 8);
Map.setOptions("HYBRID");

// 3. Define start and endyear: 
// - for this dataset, the timeframe needs to be limited to 2 years (otherwise the export dataset is too big)
// --> Decide on start and end year, uncomment export lines at end of the script, run task, repeat with new year of interest
var startyear = 2000;
var endyear = 2002;

// 4. Filter for the timeframe of interest:
var chirpday= chirpdata_pentad.filterDate(startyear + '-05-01', endyear + '-09-01'); 
var chirpdaily = chirpday.map(function(img){
  return img.copyProperties(img,['system:time_start','system:time_end']);
});

// 5. Stack all the images together:
var chirpstack = ee.Image(chirpdaily.iterate(function(imgIn, imgOut){
  var dateString = ee.Date(imgIn.get('system:time_start')).format('yyyy-MM-dd');
  return ee.Image(imgOut).addBands(imgIn.select(['precipitation']).clip(roi_reg)
          .rename(dateString));
},ee.Image()));
print('CHIRPS collection', chirpstack);

// 6. Extract the Precipitation data:
var chirps_done = chirpstack.reduceRegions({
  collection: points,
  reducer: ee.Reducer.first(),
  scale: 250,});
 print("This works:", chirps_done.limit(1));

// 7. Export the data to csv.
Export.table.toDrive({
  collection: chirps_done,
  description: 'chirpstack_'+ startyear + '_' + endyear+ '_50k',
  folder: 'gee_oaf',
  fileFormat: 'CSV'
});

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// 2.2 - Extract EVI for 50k points
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// Purpose: Extract and export EVI data for crop cut locations in study region in Kenya 
// Input files: MODIS Aqua, 50k geopoints randomly placed on maize cropland (GEE asset) 
// Export files: Four CSV files of EVI data, each for 5 (or 6) years (2000-2005, 2006-2010, 2011-2015, 2016-2020)

// Comment: The script is organized in four sections (marked by lines with five %-signs) to export the four datasets.

// 1. Import necessary datasets 
// - 50k random geopoints on maize cropland: import the asset, 
// - polygon of region of interest
var roi_reg = ee.Geometry.Polygon(
        [[[33.57858523739146, 1.627217257564576],
          [33.57858523739146, -1.5803713344511332],
          [38.21481570614146, -1.5803713344511332],
          [38.21481570614146, 1.627217257564576]]], null, false);
Map.centerObject(roi_reg, 8);
var points = ee.FeatureCollection("users/your_username/maizePoints");

// 2. Visualize inputs on Map: 
// Choose one image:
var displayMODIS = ee.ImageCollection("MODIS/061/MYD13Q1").filterDate('2010-08-01', '2010-08-30').select('EVI').first();
// MODIS EVI is scaled by 0.0001, so we need to rescale it
var eviRescaled = displayMODIS.multiply(0.0001).clip(roi_reg);  // use your ROI
// Visualization parameters (white and yellow - bare soil or sparse vegetation; green to darkgreen - higher vegetation density)
var visParams = {
  min: 0,
  max: 1,
  palette: ['white', 'yellow', 'green', 'darkgreen']
};
// Display on the map
Map.addLayer(eviRescaled, visParams, 'EVI MODIS AQUA (Aug 2010)');
Map.addLayer(roi_reg, {}, 'ROI');
Map.addLayer(points, {}, 'Crop cut locations');
Map.centerObject(points.geometry(), 8);
Map.setOptions("HYBRID");

// 3. Select EVI from general MODIS collection:   
var modis_aqua = ee.ImageCollection('MODIS/061/MYD13Q1').select('EVI');

// 4. Filter data in two steps (2000-2010 first and in the later section 2011-2020) to place in two csvs 
var startyear = 2000;
var endyear = 2005;
var MODISeviday= modis_aqua.filterDate(startyear + '-04-01', endyear + '-09-25');

// 5. Extract the imagery date from each image
var modEVIday = MODISeviday.map(function(img){
  return img.copyProperties(img,['system:time_start','system:time_end']);
});

// 6. Stack all the images together
var evistack = ee.Image(modEVIday.iterate(function(imgIn, imgOut){
  var dateString = ee.Date(imgIn.get('system:time_start')).format('yyyy-MM-dd');
  return ee.Image(imgOut).addBands(imgIn.select(['EVI']).clip(roi_reg)
          .rename(dateString));
},ee.Image()));
print('EVI stack', evistack);

// 7. Extract the EVI points for every crop cut location
var EVI = evistack.reduceRegions({
  collection: points,
  reducer: ee.Reducer.first(),
  scale: 250,});
print("This works:", EVI.limit(1));

// Uncomment this block if you want to export the data: 
// 8. Export the data to csv.
Export.table.toDrive({
  collection: EVI,
  description: 'evistack_50k_'+ startyear + '_' + endyear,
  folder: 'gee_oaf',
  fileFormat: 'CSV'
});

// %%%%%
// Repeat for different years: 
// 1. Filter data in two steps (now 2011-2020) to place in two csvs (one is too large)
var startyear = 2006;
var endyear = 2010;
var MODISeviday= modis_aqua.filterDate(startyear + '-04-01', endyear + '-09-25');



// 2. Extract the imagery date from each image
var modEVIday = MODISeviday.map(function(img){
  return img.copyProperties(img,['system:time_start','system:time_end']);
});



// 3. Stack all the images together
var evistack = ee.Image(modEVIday.iterate(function(imgIn, imgOut){
  var dateString = ee.Date(imgIn.get('system:time_start')).format('yyyy-MM-dd');
  return ee.Image(imgOut).addBands(imgIn.select(['EVI']).clip(roi_reg)
          .rename(dateString));
},ee.Image()));

print('EVI stack (2nd period)', evistack);



// 4. Extract the EVI for every crop cut location
var EVI = evistack.reduceRegions({
  collection: points,
  reducer: ee.Reducer.first(),
  scale: 250,});
  
print("This works:", EVI.limit(1));

// Uncomment this block for export. 
// 5. Export the data to csv.
Export.table.toDrive({
  collection: EVI,
  description: 'evistack_50k_'+ startyear + '_' + endyear,
  folder: 'gee_oaf',
  fileFormat: 'CSV'
});

// %%%%%
// Repeat for different years: 
// 1. Filter data in two steps (now 2011-2020) to place in two csvs (one is too large)
var startyear = 2011;
var endyear = 2015;
var MODISeviday= modis_aqua.filterDate(startyear + '-04-01', endyear + '-09-25');

// 2. Extract the imagery date from each image
var modEVIday = MODISeviday.map(function(img){
  return img.copyProperties(img,['system:time_start','system:time_end']);
});

// 3. Stack all the images together
var evistack = ee.Image(modEVIday.iterate(function(imgIn, imgOut){
  var dateString = ee.Date(imgIn.get('system:time_start')).format('yyyy-MM-dd');
  return ee.Image(imgOut).addBands(imgIn.select(['EVI']).clip(roi_reg)
          .rename(dateString));
},ee.Image()));
print('EVI stack (3rd period)', evistack);

// 4. Extract the EVI for every crop cut location
var EVI = evistack.reduceRegions({
  collection: points,
  reducer: ee.Reducer.first(),
  scale: 250,});
print("This works:", EVI.limit(1));

// Uncomment this block for export. 
// 5. Export the data to csv.
Export.table.toDrive({
  collection: EVI,
  description: 'evistack_50k_'+ startyear + '_' + endyear,
  folder: 'gee_oaf',
  fileFormat: 'CSV'
});

// %%%%%
// Repeat for different years: 
// 1. Filter data in two steps (now 2011-2020) to place in two csvs (one is too large)
var startyear = 2016;
var endyear = 2020;
var MODISeviday= modis_aqua.filterDate(startyear + '-04-01', endyear + '-09-25');

// 2. Extract the imagery date from each image
var modEVIday = MODISeviday.map(function(img){
  return img.copyProperties(img,['system:time_start','system:time_end']);
});

// 3. Stack all the images together
var evistack = ee.Image(modEVIday.iterate(function(imgIn, imgOut){
  var dateString = ee.Date(imgIn.get('system:time_start')).format('yyyy-MM-dd');
  return ee.Image(imgOut).addBands(imgIn.select(['EVI']).clip(roi_reg)
          .rename(dateString));
},ee.Image()));
print('EVI stack (4th period)', evistack);

// 4. Extract the EVI for every crop cut location
var EVI = evistack.reduceRegions({
  collection: points,
  reducer: ee.Reducer.first(),
  scale: 250,});
print("This works:", EVI.limit(1));

// Uncomment this block for export. 
// 5. Export the data to csv.
Export.table.toDrive({
  collection: EVI,
  description: 'evistack_50k_'+ startyear + '_' + endyear,
  folder: 'gee_oaf',
  fileFormat: 'CSV'
});

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// 2.3 - Extract NDVI for 50k points
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// Purpose: Extract and export NDVI data for crop cut locations in study region in Kenya 
// Input files: MODIS Aqua, 50k geopoints randomly placed on maize cropland (GEE asset) 
// Export files: Four CSV files of NDVI data, each for 10 years (2000-2005, 2006-2010,2011-2015, 2016-2020)

// Comment: The script is organized in four sections (marked by lines with five %-signs) to export the four datasets.

// 1. Import necessary datasets 
// - 50k geopoints: import the asset, 
// - polygon of region of interest
var roi_reg = ee.Geometry.Polygon(
        [[[33.57858523739146, 1.627217257564576],
          [33.57858523739146, -1.5803713344511332],
          [38.21481570614146, -1.5803713344511332],
          [38.21481570614146, 1.627217257564576]]], null, false);
Map.centerObject(roi_reg, 8);
var points = ee.FeatureCollection("users/your_username/maizePoints");

// 2. Visualize inputs on Map: 
// Choose one image:
var displayMODIS = ee.ImageCollection("MODIS/061/MYD13Q1").filterDate('2010-08-01', '2010-08-30').select('NDVI').first();
// MODIS NDVI is scaled by 0.0001, so we need to rescale it
var ndviRescaled = displayMODIS.multiply(0.0001).clip(roi_reg);  // use your ROI
// Visualization parameters (white and yellow - bare soil or sparse vegetation; green to darkgreen - higher vegetation density)
var visParams = {
  min: 0,
  max: 1,
  palette: ['white', 'yellow', 'green', 'darkgreen']
};
// Display on the map
Map.addLayer(ndviRescaled, visParams, 'NDVI MODIS AQUA (Aug 2010)');
Map.addLayer(roi_reg, {}, 'ROI');
Map.addLayer(points, {}, 'Crop cut locations');
Map.centerObject(points.geometry(), 8);
Map.setOptions("HYBRID");

// 3. Select NDVI from general MODIS collection:   
var modis_aqua = ee.ImageCollection('MODIS/061/MYD13Q1').select('NDVI');

// 4. Filter data in two steps (2000-2010 first and in the later section 2011-2020) to place in two csvs 
var startyear = 2000;
var endyear = 2005;
var MODISndviday= modis_aqua.filterDate(startyear + '-04-01', endyear + '-09-25');

// 5. Extract the imagery date from each image
var modNDVIday = MODISndviday.map(function(img){
  return img.copyProperties(img,['system:time_start','system:time_end']);
});

// 6. Stack all the images together
var ndvistack = ee.Image(modNDVIday.iterate(function(imgIn, imgOut){
  var dateString = ee.Date(imgIn.get('system:time_start')).format('yyyy-MM-dd');
  return ee.Image(imgOut).addBands(imgIn.select(['NDVI']).clip(roi_reg)
          .rename(dateString));
},ee.Image()));
print('NDVI stack', ndvistack);

// 7. Extract the NDVI points for every crop cut location
var NDVI = ndvistack.reduceRegions({
  collection: points,
  reducer: ee.Reducer.first(),
  scale: 250,});
print("This works:", NDVI.limit(1));

//Uncomment this block if you want to export the data: 
// 8. Export the data to csv.
Export.table.toDrive({
  collection: NDVI,
  description: 'ndvistack_50k_' + startyear + '_' + endyear,
  folder: 'gee_oaf',
  fileFormat: 'CSV'
});

// %%%%%
// Repeat for different years: 

// 1. Filter data in two steps (now 2011-2020) to place in two csvs (one is too large)
var startyear = 2006;
var endyear = 2010;
var MODISndviday= modis_aqua.filterDate(startyear + '-04-01', endyear + '-09-25');

// 2. Extract the imagery date from each image
var modNDVIday = MODISndviday.map(function(img){
  return img.copyProperties(img,['system:time_start','system:time_end']);
});

// 3. Stack all the images together
var ndvistack = ee.Image(modNDVIday.iterate(function(imgIn, imgOut){
  var dateString = ee.Date(imgIn.get('system:time_start')).format('yyyy-MM-dd');
  return ee.Image(imgOut).addBands(imgIn.select(['NDVI']).clip(roi_reg)
          .rename(dateString));
},ee.Image()));
print('NDVI stack (2nd period)', ndvistack);

// 4. Extract the NDVI for every crop cut location
var NDVI = ndvistack.reduceRegions({
  collection: points,
  reducer: ee.Reducer.first(),
  scale: 250,});
print("This works:", NDVI.limit(1));

// Uncomment this block for export. 
// 5. Export the data to csv.
Export.table.toDrive({
  collection: NDVI,
  description: 'ndvistack_50k_' + startyear + '_' + endyear,
  folder: 'gee_oaf',
  fileFormat: 'CSV'
});

// %%%%%
// Repeat for different years: 

// 1. Filter data in two steps (now 2011-2020) to place in two csvs (one is too large)
var startyear = 2011;
var endyear = 2015;
var MODISndviday= modis_aqua.filterDate(startyear + '-04-01', endyear + '-09-25');

// 2. Extract the imagery date from each image
var modNDVIday = MODISndviday.map(function(img){
  return img.copyProperties(img,['system:time_start','system:time_end']);
});

// 3. Stack all the images together
var ndvistack = ee.Image(modNDVIday.iterate(function(imgIn, imgOut){
  var dateString = ee.Date(imgIn.get('system:time_start')).format('yyyy-MM-dd');
  return ee.Image(imgOut).addBands(imgIn.select(['NDVI']).clip(roi_reg)
          .rename(dateString));
},ee.Image()));
print('NDVI stack (3rd period)', ndvistack);

// 4. Extract the NDVI for every crop cut location
var NDVI = ndvistack.reduceRegions({
  collection: points,
  reducer: ee.Reducer.first(),
  scale: 250,});
print("This works:", NDVI.limit(1));

// Uncomment this block for export. 
// 5. Export the data to csv.
Export.table.toDrive({
  collection: NDVI,
  description: 'ndvistack_50k_' + startyear + '_' + endyear,
  folder: 'gee_oaf',
  fileFormat: 'CSV'
});

// %%%%%
// Repeat for different years: 
// 1. Filter data in two steps (now 2011-2020) to place in two csvs (one is too large)
var startyear = 2016;
var endyear = 2020;
var MODISndviday= modis_aqua.filterDate(startyear + '-04-01', endyear + '-09-25');

// 2. Extract the imagery date from each image
var modNDVIday = MODISndviday.map(function(img){
  return img.copyProperties(img,['system:time_start','system:time_end']);
});

// 3. Stack all the images together
var ndvistack = ee.Image(modNDVIday.iterate(function(imgIn, imgOut){
  var dateString = ee.Date(imgIn.get('system:time_start')).format('yyyy-MM-dd');
  return ee.Image(imgOut).addBands(imgIn.select(['NDVI']).clip(roi_reg)
          .rename(dateString));
},ee.Image()));
print('NDVI stack (4th period)', ndvistack);

// 4. Extract the NDVI for every crop cut location
var NDVI = ndvistack.reduceRegions({
  collection: points,
  reducer: ee.Reducer.first(),
  scale: 250,});
print("This works:", NDVI.limit(1));

// Uncomment this block for export. 
// 5. Export the data to csv.
Export.table.toDrive({
  collection: NDVI,
  description: 'ndvistack_50k_' + startyear + '_' + endyear,
  folder: 'gee_oaf',
  fileFormat: 'CSV'
});

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// 2.4 - Extract temp for 50k points
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// Purpose: Extract and export temperature data for crop cut locations in study region in Kenya 
// Input files: ERA5 Daily, 50k geopoints randomly placed on maize cropland (GEE asset) 
// Export files: CSV file of temperature data from one season at a time (2000-2020)

// 1. Define start and endyear
// - for this dataset, the timeframe needs to be limited to one year (otherwise the export dataset is too big)
// --> Decide on year, uncomment export lines at end of the script, run task, repeat with new year of interest
var year=2000;

// 2. Import necessary files:
// - 50k geopoints: import the asset, 
// - temp data: from ERA5 
// - polygon of region of interest
var roi_reg = ee.Geometry.Polygon(
        [[[33.57858523739146, 1.627217257564576],
          [33.57858523739146, -1.5803713344511332],
          [38.21481570614146, -1.5803713344511332],
          [38.21481570614146, 1.627217257564576]]], null, false);
var points = ee.FeatureCollection("users/your_username/maizePoints");
var temp=ee.ImageCollection("ECMWF/ERA5_LAND/DAILY_AGGR").filterDate(year+'-05-01', year+'-09-01');

// 3. Visualize imported datasets 
var meantemp = temp.select('temperature_2m').filterDate(year + '-05-01',year+ '-09-01');
var meanTempSeason = meantemp.mean().subtract(273.15).clip(roi_reg);  // Average over time, convert to °C
Map.addLayer(meanTempSeason, 
  {min: 10, max: 35, palette: ['blue', 'green', 'yellow', 'orange', 'red']},
  'Mean Temp (°C, May–July)');
Map.addLayer(points,{},'Crop cut location' ); 
Map.addLayer(roi_reg,{}, 'ROI'); 
Map.centerObject(points.geometry(), 8);
//Map.setOptions("HYBRID");
print('temp', meantemp);

// 4. Compute GDD: 
function gdd(temp){
  var mintemp = temp.select('temperature_2m_min').subtract(273.15);
  var maxtemp = temp.select('temperature_2m_max').subtract(273.15);
  var meantemp= temp.select('temperature_2m').subtract(273.15);
  var doubletemp = mintemp.add(maxtemp);
  var meanDouble = doubletemp.divide(2);
  var meanDoubleAdj = meanDouble.where(meanDouble.gt(29), 29);
  var GDD = meanDoubleAdj.subtract(10); 
  var GDDNonNeg = GDD.where(GDD.lt(0), 0).rename('GDD');
  return GDDNonNeg.addBands(meantemp); 
  //return GDDNonNeg.addBands(mintemp).addBands(maxtemp); 
}
// Map the function over the collection
var gddtemp = temp.map(gdd);
print('gdd', gddtemp);

// 5. Separate out GDD and meantemp to reduce number of bands for reducer
// And extract the imagery date from each image:
var gdd_only = gddtemp.map(function(img) {
  return img.select('GDD').copyProperties(img, ['system:time_start', 'system:time_end']);
});
var meantemp_only = gddtemp.map(function(img) {
  return img.select('temperature_2m').copyProperties(img, ['system:time_start', 'system:time_end']);
});
print('GDD stack', gdd_only);
print('Meantemp stack', meantemp_only);

// 6. Stack the images together (GDD and meantemp separately)
var gdd_stack = gdd_only.toBands().clip(roi_reg);
var meantemp_stack = meantemp_only.toBands().clip(roi_reg);

// 7. Extract the values for all points in the dataset
var gdd_points = gdd_stack.reduceRegions({
  collection: points,
  reducer: ee.Reducer.first(),
  scale: 250
});
var meantemp_points = meantemp_stack.reduceRegions({
  collection: points,
  reducer: ee.Reducer.first(),
  scale: 250
});

// 8. Export both tables to Drive
Export.table.toDrive({
  collection: gdd_points,
  description: 'gdd_'+ year +'_50k',
  folder: 'gee_oaf',
  fileFormat: 'CSV'
});
Export.table.toDrive({
  collection: meantemp_points,
  description: 'meantemp_'+ year +'_50k',
  folder: 'gee_oaf',
  fileFormat: 'CSV'
});

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// 2.5 - Extract soil for 50k points
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// Purpose: Extract and export soil data for crop cut locations in study region in Kenya 
// Input files: Different soil datasets saved as GEE assets, 50k geopoints randomly placed on maize cropland (GEE asset) 
// Export files: CSV file of soil data, time-invariant data, hence only one output file

// 1. Import all the necessary datasets: 
// A. Soil data first: 
    var awcp30 = ee.Image ("users/kirchnere/AWC_1km"); 
        //based on Jin's file which was named: projects/jin-digitalag-lab/AfricaSoilGrids/af_agg_30cm_TAWCpF23mm__M_1km, but fresh download on 23-10-11 from https://files.isric.org/public/gyga/GYGA_Results/
    var erzd = ee.Image ("users/kirchnere/ERZD_1km");
        //based on Jin's file which was named: projects/jin-digitalag-lab/AfricaSoilGrids/af_ERZD_M_1km, but fresh download on 23-10-11 from https://files.isric.org/public/gyga/GYGA_Results/
    var al = ee.Image ("projects/jin-digitalag-lab/AfricaSoilGrids/af250m_nutrient_al_m_agg30cm");
    var soc = ee.Image ("projects/jin-digitalag-lab/AfricaSoilGrids/af_ORCDRC_T__M_sd1_250m");
    var totN = ee.Image ("projects/jin-digitalag-lab/AfricaSoilGrids/af250m_nutrient_n_m_agg30cm");
// B. Geopoints of crop cuts: 
     var points = ee.FeatureCollection("users/your_username/maizePoints");

// 2. Combine all the bands: 
var soil = ee.Image([])
  .addBands(erzd.rename('erzd'))
  .addBands(awcp30.rename('awcp30'))
  .addBands(al.rename('al'))
  .addBands(soc.rename('soc'))
  .addBands(totN.rename('totN'));
print(soil);

// 3. Reduce the datasets to the values of the geopoints: 
var soilpoints = soil.reduceRegions({
  collection: points,
 reducer: ee.Reducer.first(),
  scale: 250,});

//4. Export soil data:  
Export.table.toDrive({
  collection: soilpoints,
  description: 'soil_50k',
  folder: 'gee_oaf',
  fileFormat: 'CSV'
});
