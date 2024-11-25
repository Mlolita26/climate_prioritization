if(!require("pacman", character.only = TRUE)){install.packages("pacman",dependencies = T)}

pacman::p_load(data.table,terra,sf,ecmwfr,pbapply,install = T)

# Rainfall CV ####
rain_cv<-terra::rast("raw_data/chirps_cv/cv_2010-2023_annual.tif")

# Classify
# Create a classification matrix for terra::classify
# Columns: from, to, becomes
classification_matrix <- matrix(c(
  0, 20, 1,   # Low Variability
  20, 50, 2,  # Moderate Variability
  50, 80, 3,  # High Variability
  80, 200, 4  # Extreme Variability
), ncol = 3, byrow = TRUE)

# Use classify with the classification matrix
rain_cv_class <- classify(rain_cv, classification_matrix)

# Assign levels (optional, to add labels)
cv_levels <-data.frame(value =1:4, label =  c("Low Variability", "Moderate Variability", "High Variability", "Extreme Variability"))
levels(rain_cv_class)<-cv_levels

# Aqueduct ####
aqueduct_url<-"https://files.wri.org/aqueduct/aqueduct-4-0-water-risk-data.zip"
bws_file<-"raw_data/aqueduct/bws_baseline.tif"
rfr_file<-"raw_data/aqueduct/rfr_baseline.tif"
bwd_file<-"raw_data/aqueduct/bwd_baseline.tif"
cfr_file<-"raw_data/aqueduct/cfr_baseline.tif"
sev_file<-"raw_data/aqueduct/sev_baseline.tif"

if(!file.exists(bws_file)){
aqueduct_gdb<-"raw_data/aqueduct/Aqueduct40_waterrisk_download_Y2023M07D05/GDB/Aq40_Y2023D07M05.gdb"
  
sf::st_layers(aqueduct_gdb)
aq_base_annual<-sf::st_read(aqueduct_gdb,layer="baseline_annual")
aq_base_annual<-terra::vect(aq_base_annual)
# https://github.com/wri/Aqueduct40/blob/master/data_dictionary_country-rankings.md#country-baseline

bws<-terra::rasterize(aq_base_annual,rain_cv,field="bws_cat")
rast_levels<-unique(data.table(data.frame(aq_base_annual))[,.(bws_cat,bws_label)])
levels(bws) <- data.frame(value = rast_levels$bws_cat, label = rast_levels$bws_label)
terra::writeRaster(bws,filename=bws_file,overwrite=T)

bwd<-terra::rasterize(aq_base_annual,rain_cv,field="bwd_cat")
rast_levels<-unique(data.table(data.frame(aq_base_annual))[,.(bwd_cat,bwd_label)])
levels(bwd) <- data.frame(value = rast_levels$bwd_cat, label = rast_levels$bwd_label)
terra::writeRaster(bwd,filename=bwd_file,overwrite=T)

rfr<-terra::rasterize(aq_base_annual,rain_cv,field="rfr_cat")
rast_levels<-unique(data.table(data.frame(aq_base_annual))[,.(rfr_cat,rfr_label)])[order(rfr_cat)]
levels(rfr) <- data.frame(value = rast_levels$rfr_cat, label = rast_levels$rfr_label)
terra::writeRaster(rfr,filename=rfr_file,overwrite=T)

cfr<-terra::rasterize(aq_base_annual,rain_cv,field="cfr_cat")
rast_levels<-unique(data.table(data.frame(aq_base_annual))[,.(cfr_cat,cfr_label)])
levels(cfr) <- data.frame(value = rast_levels$cfr_cat, label = rast_levels$cfr_label)
terra::writeRaster(cfr,filename=cfr_file,overwrite=T)

sev<-terra::rasterize(aq_base_annual,rain_cv,field="sev_cat")
rast_levels<-unique(data.table(data.frame(aq_base_annual))[,.(sev_cat,sev_label)])
levels(sev) <- data.frame(value = rast_levels$sev_cat, label = rast_levels$sev_label)
terra::writeRaster(sev,filename=sev_file,overwrite=T)


}else{
  bws<-rast(bws_file)
  bwd<-rast(bwd_file)
  rfr<-rast(rfr_file)
  cfr<-rast(cfr_file)
  sev<-rast(sev_file)
}

# JRC flood hazard map ####
# values are flood depth in m
# 100 yr return period
# https://data.jrc.ec.europa.eu/dataset/jrc-floods-floodmapgl_rp100y-tif
# https://jeodpp.jrc.ec.europa.eu/ftp/jrc-opendata/CEMS-GLOFAS/flood_hazard/

jrc_floodrisk_100y<-terra::rast("raw_data/jrc_floods/floodMapGL_rp100y.tif")
jrc_floodrisk_100y<-terra::resample(jrc_floodrisk_100y,rain_cv)

# 20 year return period
# https://data.jrc.ec.europa.eu/dataset/jrc-floods-floodmapgl_rp20y-tif
jrc_floodrisk_20y<-terra::rast("raw_data/jrc_floods/floodMapGL_rp20y.tif")
jrc_floodrisk_20y<-terra::resample(jrc_floodrisk_20y,rain_cv)

# Global drought observatory ####
  ## Read in drought #####
files<-list.files("raw_data/drought_observatory","tif$",full.names = T,recursive = T)
layer_names<-gsub("_t.tif","",gsub("rdria_m_wld_","",basename(files)))
layer_names<-data.table(year=substr(layer_names,1,4),
                        month=substr(layer_names,5,6),
                        day=substr(layer_names,7,8),
                        file=files)
years<-unique(substr(layer_names,1,4))

layer_names<-layer_names[order(year,month,day)]
drought<-terra::rast(layer_names$file)
drought<-drought+0

  ## Read in crop calendars  ####
# Function to convert Julian day to month
convert_to_month <- function(julian_day) {
  julian_day<-as.numeric(julian_day)
  date <- as.Date(julian_day - 1, origin = paste0(year, "-01-01"))
  as.numeric(format(date, "%m"))  # Return numeric month
}

    ### Rainfed  #####
cc_maize_rf<-terra::rast("raw_data/gccmi/mai_rf_ggcmi_crop_calendar_phase3_v1.01.nc4")
# Resample to drought
cc_maize_rf<-terra::resample(cc_maize_rf,drought,method="near")

# Convert days to month
cc_maize_rf$planting_month<- app(cc_maize_rf$planting_day, convert_to_month)
cc_maize_rf$maturity_month<- app(cc_maize_rf$maturity_day, convert_to_month)
# Optional: Add names for clarity
names(cc_maize_rf$planting_month) <- "month"
names(cc_maize_rf$maturity_month) <- "month"

    ### Irrigated ######
cc_maize_ir<-terra::rast("raw_data/gccmi/mai_ir_ggcmi_crop_calendar_phase3_v1.01.nc4")
# Resample to drought
cc_maize_ir<-terra::resample(cc_maize_ir,drought,method="near")

# Convert days to month
cc_maize_ir$planting_month<- app(cc_maize_ir$planting_day, convert_to_month)
cc_maize_ir$maturity_month<- app(cc_maize_ir$maturity_day, convert_to_month)
# Optional: Add names for clarity
names(cc_maize_ir$planting_month) <- "month"
names(cc_maize_ir$maturity_month) <- "month"

  ## Summarize drought for cropping seasons ####
    ### Rainfed ######
plant<-(cc_maize_rf$planting_month*3)-2
harvest<-cc_maize_rf$maturity_month*3
# Where plant>harvest (e.g. plant = 11 harvest = 3) add 12 to harvest (e.g. plant = 11 harvest = 15)
harvest[plant[]>harvest[]]<-harvest[plant[]>harvest[]]+36

plant_min<-min(plant[],na.rm=T)
harvest_max<-max(harvest[],na.rm=T)

years<-layer_names[,unique(year)]

drought_yrs_rf<-terra::rast(pblapply(1:(length(years)-1),FUN=function(m){

  # Subset drought to increase efficiency
  x = drought[[(plant_min+36*(m-1)):(harvest_max+36*(m-1))]]
  
  drought1<-terra::rapp(x,
                         first=plant,
                         last=harvest,
                         fun="mean",
                         na.rm=T)
  
  names(drought1)<-years[m]
  drought1
}))

# Calculate proportion of seasons where mean drought risk is >1 (1 = low)
d_fun<-function(x,threshold){
  x<-x[!is.na(x)]
  x<-sum(x>threshold)/length(x)
  return(x)
}

drought_rf<-app(drought_yrs_rf,d_fun,threshold=1)

    ### Irrigated ######
plant<-(cc_maize_ir$planting_month*3)-2
harvest<-cc_maize_ir$maturity_month*3
# Where plant>harvest (e.g. plant = 11 harvest = 3) add 12 to harvest (e.g. plant = 11 harvest = 15)
harvest[plant[]>harvest[]]<-harvest[plant[]>harvest[]]+36

plant_min<-min(plant[],na.rm=T)
harvest_max<-max(harvest[],na.rm=T)

years<-layer_names[,unique(year)]

drought_yrs_ir<-terra::rast(pblapply(1:(length(years)-1),FUN=function(m){
  
  # Subset drought to increase efficiency
  x = drought[[(plant_min+36*(m-1)):(harvest_max+36*(m-1))]]
  
  drought1<-terra::rapp(x,
                        first=plant,
                        last=harvest,
                        fun="mean",
                        na.rm=T)
  
  names(drought1)<-years[m]
  drought1
}))
drought_ir<-app(drought_yrs_ir,d_fun,threshold=1)

# Agroclim ####
  # https://cds.climate.copernicus.eu/datasets/sis-agroclimatic-indicators?tab=overview
  # WSDI - warm spell duration index
  wsdi<-terra::rast("raw_data/agroclim/CDD_C3S-glob-agric_WFDEI_hist_season_19810101-20101231_v1.1.nc")
  # data is 3-monthly
  time(wsdi)
  
# Combine Hazards ####

drought_rf

## Rain cv ####
var<-rain_cv_class
var[var[]==1]<-NA
var[var[]>1]<-1
levels(var)<-data.frame(value=1,label="v")
names(var)<-"v = rainfall_cv>20"

## Flood ####
flood<-jrc_floodrisk_20y
flood[flood[]>0]<-10
levels(flood)<-data.frame(value=10,label="f")
names(flood)<-"f = flood_risk_20y"

# Drought ####
dry_ir<-drought_ir
dry_ir<-resample(dry_ir,rain_cv)
dry_ir[dry_ir[]<0.25]<-NA
dry_ir[dry_ir[]>=0.25]<-100
levels(dry_ir)<-data.frame(value=100,label="d")
names(dry_ir)<-"d = drought risk>0.25"


# Combine hazards
haz_comb<-sum(c(var,flood,dry_ir),na.rm=T)

# Add labels
haz_levels<-data.frame(
  data.table(expand.grid(var=c(0,1),flood=c(0,10),dry=c(0,100)))[,.(value=var+dry+flood)],
  data.table(expand.grid(var=c(NA,"v"),flood=c(NA,"f"),dry=c(NA,"d")))[,.(label=paste(na.omit(c(var,flood,dry)),collapse = "+")),by=.(var,flood,dry)][,"label"]
)

levels(haz_comb)<-haz_levels
names(haz_comb)<-"combined hazards"

plot(c(haz_comb,var,flood,dry_ir))
haz_comb_f<-data.table(freq(haz_comb))[order(count,decreasing = T)]
haz_comb_f[,perc:=round(100*(count/sum(count)),1)]
