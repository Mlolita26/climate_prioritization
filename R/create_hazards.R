if(!require("pacman", character.only = TRUE)){install.packages("pacman",dependencies = T)}

pacman::p_load(data.table,terra,sf,ecmwfr,pbapply,viridis,MetBrewer,wesanderson,install = T)

Viridis<-data.table(Source="viridis",Palette=c("magma","inferno","plasma","viridis","cividis","rocket","mako","turbo"))
Met<-data.table(Source="MetBrewer",Palette=names(MetBrewer::MetPalettes))
Wes<-data.table(Source="wesanderson",Palette=names(wesanderson::wes_palettes))

Palettes<-rbind(Viridis,Met,Wes)

PalFun<-function(PalName,N,Names) {
  Viridis<-data.table(Source="viridis",Palette=c("magma","inferno","plasma","viridis","cividis","rocket","mako","turbo"))
  Met<-data.table(Source="MetBrewer",Palette=names(MetBrewer::MetPalettes))
  Palettes<-rbind(Viridis,Met)
  
  if(Palettes[Palette==PalName,Source]=="viridis"){
    PAL<-viridis::viridis(N,option=PalName)
  }
  
  if(Palettes[Palette==PalName,Source]=="MetBrewer"){
    PAL<-MetBrewer::met.brewer(name=PalName, n=N, type="continuous")
  }
  
  if(Palettes[Palette==PalName,Source]=="Wes"){
    PAL<-wesanderson::wes_palette(name=PalName, n=N, type="continuous")
  }
  names(PAL)<-Names
  
  return(PAL)
}

# GCCMI crop calendars  ####
# Function to convert Julian day to month
convert_to_month <- function(julian_day) {
  julian_day<-as.numeric(julian_day)
  date <- as.Date(julian_day - 1, origin = paste0(year, "-01-01"))
  as.numeric(format(date, "%m"))  # Return numeric month
}

  ## Rainfed  #####
  cc_maize_rf<-terra::rast("raw_data/gccmi/mai_rf_ggcmi_crop_calendar_phase3_v1.01.nc4")
  # Resample to drought
  cc_maize_rf<-terra::resample(cc_maize_rf,drought,method="near")
  
  # Convert days to month
  cc_maize_rf$planting_month<- app(cc_maize_rf$planting_day, convert_to_month)
  cc_maize_rf$maturity_month<- app(cc_maize_rf$maturity_day, convert_to_month)
  # Optional: Add names for clarity
  names(cc_maize_rf$planting_month) <- "month"
  names(cc_maize_rf$maturity_month) <- "month"
  
  ## Irrigated ######
  cc_maize_ir<-terra::rast("raw_data/gccmi/mai_ir_ggcmi_crop_calendar_phase3_v1.01.nc4")
  # Resample to drought
  cc_maize_ir<-terra::resample(cc_maize_ir,drought,method="near")
  
  # Convert days to month
  cc_maize_ir$planting_month<- app(cc_maize_ir$planting_day, convert_to_month)
  cc_maize_ir$maturity_month<- app(cc_maize_ir$maturity_day, convert_to_month)
  # Optional: Add names for clarity
  names(cc_maize_ir$planting_month) <- "month"
  names(cc_maize_ir$maturity_month) <- "month"
  
# Rainfall CV ####
# 2010 - 2023
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

# Agroclim - WSDI ####
  # https://cds.climate.copernicus.eu/datasets/sis-agroclimatic-indicators?tab=overview
  # The Warm Spell Duration Index (WSDI) measures the number of consecutive days with maximum temperatures above the 90th percentile of a reference period. For crop production, stress thresholds depend on the specific crop, its growth stage, and local climatic conditions. However, general guidelines can be provided:
  #  General Thresholds for WSDI and Crop Stress
    #1.	Baseline Stress	WSDI > 6 days in a season is often indicative of heat stress, especially for heat-sensitive crops.
    #2.	Moderate Stress:WSDI > 15 days: Many crops begin to show moderate reductions in yield due to prolonged exposure to high temperatures.
    #3.	Severe Stress: WSDI > 30 days: Prolonged heat waves can cause severe damage to most crops, leading to significant yield losses or even crop failure.
  wsdi<-terra::rast("raw_data/agroclim/CDD_C3S-glob-agric_WFDEI_hist_season_19810101-20101231_v1.1.nc")
  # data is 3-monthly
  time(wsdi)
  # First Q of 1981 appears to be missing, remove first Qs
  wsdi<-wsdi[[-(1:3)]]
  years<-as.numeric(unique(format(time(wsdi),"%Y")))
  
  # Repeat each wdsi quarter 3 times
  wsdi_m<-wsdi[[rep(1:nlyr(wsdi),each=3)]]

  ## Irrigated #####
  plant<- cc_maize_ir$planting_month
  harvest <- cc_maize_ir$maturity_month
  # Resample to wsdi layer
  plant<-terra::resample(plant,wsdi,method="near")
  harvest<-terra::resample(harvest,wsdi,method="near")
  
  # Where plant>harvest (e.g. plant = 11 harvest = 3) add 12 to harvest (e.g. plant = 11 harvest = 15)
  harvest[plant[]>harvest[]]<-harvest[plant[]>harvest[]]+12
    
  plant_min<-min(plant[],na.rm=T)
  harvest_max<-max(harvest[],na.rm=T)
  
  wsdi_yrs_ir<-terra::rast(pblapply(1:(length(years)-1),FUN=function(m){
    
    # Subset drought to increase efficiency
    x = wsdi_m[[(plant_min+12*(m-1)):(harvest_max+12*(m-1))]]
    
    y<-terra::rapp(x,
                  first=plant,
                  last=harvest,
                  fun="mean",
                  na.rm=T)
    
    names(y)<-years[m]
    y
  }))
  
  # Classify each year
  wsdi_yrs_ir_class<-classify(wsdi_yrs_ir,data.frame(from=c(0,6,15,30),to=c(6,15,30,999),becomes=c(0,1,2,3)))
  wsdi_levels<-data.frame(value=0:3,hazard=c("No Stress","Baseline Stress","Moderate Stress","Severe Stress"))
  
  for(i in 1:nlyr(wsdi_yrs_ir_class)){
  levels(wsdi_yrs_ir_class[[i]])<-wsdi_levels
  }
  names(wsdi_yrs_ir_class)<-names(wsdi_yrs_ir)
  
  # Get frequency of moderate hazard or higher
  wsdi_ir_freq<-wsdi_yrs_ir_class
  wsdi_ir_freq<-classify(wsdi_ir_freq,data.frame(from=c(0:3),to=c(0,0,1,1)))
  wsdi_ir_freq<-sum(wsdi_ir_freq,na.rm=T)/nlyr(wsdi_ir_freq)
  names(wsdi_ir_freq)<-"wsdi_mod-haz-freq_ir-cc"
  plot(wsdi_ir_freq)
  
  ## Rainfed #####

  plant<- cc_maize_rf$planting_month
  harvest <- cc_maize_rf$maturity_month
  # Resample to wsdi layer
  plant<-terra::resample(plant,wsdi,method="near")
  harvest<-terra::resample(harvest,wsdi,method="near")
  
  # Where plant>harvest (e.g. plant = 11 harvest = 3) add 12 to harvest (e.g. plant = 11 harvest = 15)
  harvest[plant[]>harvest[]]<-harvest[plant[]>harvest[]]+12
  
  plant_min<-min(plant[],na.rm=T)
  harvest_max<-max(harvest[],na.rm=T)
  
  wsdi_yrs_rf<-terra::rast(pblapply(1:(length(years)-1),FUN=function(m){
    
    # Subset drought to increase efficiency
    x = wsdi_m[[(plant_min+12*(m-1)):(harvest_max+12*(m-1))]]
    
    y<-terra::rapp(x,
                   first=plant,
                   last=harvest,
                   fun="mean",
                   na.rm=T)
    
    names(y)<-years[m]
    y
  }))
  
  # Classify each year
  wsdi_yrs_rf_class<-classify(wsdi_yrs_rf,data.frame(from=c(0,6,15,30),to=c(6,15,30,999),becomes=c(0,1,2,3)))
  wsdi_levels<-data.frame(value=0:3,hazard=c("No Stress","Baseline Stress","Moderate Stress","Severe Stress"))
  
  for(i in 1:nlyr(wsdi_yrs_rf_class)){
    levels(wsdi_yrs_rf_class[[i]])<-wsdi_levels
  }
  names(wsdi_yrs_rf_class)<-names(wsdi_yrs_rf)
  
  # Get frequency of moderate hazard or higher
  wsdi_rf_freq<-wsdi_yrs_rf_class
  wsdi_rf_freq<-classify(wsdi_rf_freq,data.frame(from=c(0:3),to=c(0,0,1,1)))
  wsdi_rf_freq<-sum(wsdi_rf_freq,na.rm=T)/nlyr(wsdi_rf_freq)
  plot(wsdi_rf_freq)
  names(wsdi_rf_freq)<-"wsdi_mod-haz-freq_rf-cc"
  
# AgERA5 - NTx35 ####
  ntx35<-terra::rast("raw_data/agera5/NTx35.tif")
  years<-unique(format(as.Date(names(ntx35)),"%Y"))
  ## Irrigated #####
  plant<- cc_maize_ir$planting_month
  harvest <- cc_maize_ir$maturity_month
  # Resample to wsdi layer
  plant<-terra::resample(plant,ntx35,method="near")
  harvest<-terra::resample(harvest,ntx35,method="near")
  
  # Where plant>harvest (e.g. plant = 11 harvest = 3) add 12 to harvest (e.g. plant = 11 harvest = 15)
  harvest[plant[]>harvest[]]<-harvest[plant[]>harvest[]]+12
    
  plant_min<-min(plant[],na.rm=T)
  harvest_max<-max(harvest[],na.rm=T)
  
  ntx35_yrs_ir<-terra::rast(pblapply(1:(length(years)-1),FUN=function(m){
    
    # Subset drought to increase efficiency
    x = ntx35[[(plant_min+12*(m-1)):(harvest_max+12*(m-1))]]
    
    y<-terra::rapp(x,
                  first=plant,
                  last=harvest,
                  fun="mean",
                  na.rm=T)
    
    names(y)<-years[m]
    y
  }))

  # Classify each year
  ntx35_yrs_ir_class<-classify(ntx35_yrs_ir,data.frame(from=c(0,10,20,25),to=c(10,20,25,999),becomes=c(0,1,2,3)))
  ntx35_levels<-data.frame(value=0:3,hazard=c("No Stress","Moderate Stress","Severe Stress","Extreme Stress"))
  
  for(i in 1:nlyr(ntx35_yrs_ir_class)){
  levels(ntx35_yrs_ir_class[[i]])<-ntx35_levels
  }
  names(ntx35_yrs_ir_class)<-names(ntx35_yrs_ir)
  
  # Get frequency of moderate hazard or higher
  ntx35_ir_freq<-ntx35_yrs_ir_class
  ntx35_ir_freq<-classify(ntx35_ir_freq,data.frame(from=c(0:3),to=c(0,1,1,1)))
  ntx35_ir_freq<-sum(ntx35_ir_freq,na.rm=T)/nlyr(ntx35_ir_freq)
  names(ntx35_ir_freq)<-"ntx35_mod-haz-freq_ir-cc"
  plot(ntx35_ir_freq)
  
  # Resample
  ntx35_ir_freq<-terra::resample(ntx35_ir_freq,rain_cv)
  
  ## Rainfed #####
  plant<- cc_maize_rf$planting_month
  harvest <- cc_maize_rf$maturity_month
  # Resample to wsdi layer
  plant<-terra::resample(plant,ntx35,method="near")
  harvest<-terra::resample(harvest,ntx35,method="near")
  
  # Where plant>harvest (e.g. plant = 11 harvest = 3) add 12 to harvest (e.g. plant = 11 harvest = 15)
  harvest[plant[]>harvest[]]<-harvest[plant[]>harvest[]]+12
  
  plant_min<-min(plant[],na.rm=T)
  harvest_max<-max(harvest[],na.rm=T)
  
  ntx35_yrs_rf<-terra::rast(pblapply(1:(length(years)-1),FUN=function(m){
    
    # Subset drought to increase efficiency
    x = ntx35[[(plant_min+12*(m-1)):(harvest_max+12*(m-1))]]
    
    y<-terra::rapp(x,
                   first=plant,
                   last=harvest,
                   fun="mean",
                   na.rm=T)
    
    names(y)<-years[m]
    y
  }))
  
  # Classify each year
  ntx35_yrs_rf_class<-classify(ntx35_yrs_rf,data.frame(from=c(0,10,20,25),to=c(10,20,25,999),becomes=c(0,1,2,3)))
  ntx35_levels<-data.frame(value=0:3,hazard=c("No Stress","Moderate Stress","Severe Stress","Extreme Stress"))
  
  for(i in 1:nlyr(ntx35_yrs_rf_class)){
    levels(ntx35_yrs_rf_class[[i]])<-ntx35_levels
  }
  names(ntx35_yrs_rf_class)<-names(ntx35_yrs_rf)
  
  # Get frequency of moderate hazard or higher
  ntx35_rf_freq<-ntx35_yrs_rf_class
  ntx35_rf_freq<-classify(ntx35_rf_freq,data.frame(from=c(0:3),to=c(0,1,1,1)))
  ntx35_rf_freq<-sum(ntx35_rf_freq,na.rm=T)/nlyr(ntx35_rf_freq)
  names(ntx35_rf_freq)<-"ntx35_mod-haz-freq_rf-cc"
  plot(ntx35_rf_freq)
  
# GAEZ - LGP #
  file<-list.files("raw_data/gaez/len-longest-component-LGP","rcp2p6_2050s",full.names = T)
  lgp_fut<-terra::rast(file)
  file<-list.files("raw_data/gaez/len-longest-component-LGP","Hist",full.names = T)
  lgp_base<-terra::rast(file)
  lgp_change<-(lgp_fut-lgp_base)
  lgp_change[lgp_change>=0]<-NA  
  lgp_change<-100*(-lgp_change/lgp_base)
# Combine Hazards ####
  ## Rain cv #####
  var<-rain_cv_class
  var[var[]==1]<-0
  var[var[]>1]<-1
  levels(var)<-data.frame(value=1,label="v")
  names(var)<-"v = rainfall_cv>20"
  
  ## Flood #####
  flood<-jrc_floodrisk_20y
  flood[flood[]>0]<-10
  levels(flood)<-data.frame(value=10,label="f")
  names(flood)<-"f = flood_risk_20y"
  
  ## Drought #####
  cutoff<-0.2
  dry<-drought_rf
  dry<-resample(dry,rain_cv)
  dry[dry[]<cutoff]<-0
  dry[dry[]>=cutoff]<-100
  levels(dry)<-data.frame(value=100,label="d")
  names(dry)<-paste0("d = drought risk>",cutoff)
  plot(dry)
  
  ## Heat #####
  cutoff<-0.2
  heat_rf<-ntx35_rf_freq
  heat_rf<-resample(heat_rf,rain_cv)
  heat_rf[heat_rf[]<cutoff]<-0
  heat_rf[heat_rf[]>=cutoff]<-1000
  levels(heat_rf)<-data.frame(value=1000,label="h")
  names(heat_rf)<-paste0("d = heat risk>",cutoff)
  plot(heat_rf)
  
  cutoff<-0.2
  heat_ir<-ntx35_ir_freq
  heat_ir<-resample(heat_ir,rain_cv)
  heat_ir[heat_ir[]<cutoff]<-0
  heat_ir[heat_ir[]>=cutoff]<-1000
  levels(heat_ir)<-data.frame(value=1000,label="h")
  names(heat_ir)<-paste0("d = heat risk>",cutoff)
  plot(heat_ir)
  
  ## Water Stress #####
  water_stress<-bwd
  water_stress<-resample(water_stress,rain_cv,method="near")
  levels(water_stress)
  water_stress<-terra::classify(water_stress,data.table(from=c(-9999,-1,0,1,2,3,4),to=c(0,0,0,0,100,100,100)))
  levels(water_stress)<-data.frame(value=100,label="ws")
  plot(water_stress)
  ## Growing Season Reduction ####
  gsr<-lgp_change
  gsr<-resample(gsr,rain_cv)
  gsr[gsr<10]<-0
  gsr[gsr>=10]<-10000
  levels(gsr)<-data.frame(value=10000,label="r")
  plot(gsr)
  
  ## Combine hazards #####
    ### Rainfed ######
     haz_comb<-sum(c(var,flood,dry,heat_rf,gsr),na.rm=T)
    # Add labels
    haz_levels<-data.frame(
      data.table(expand.grid(var=c(0,1),flood=c(0,10),dry=c(0,100),heat=c(0,1000),gsr=c(0,10000)))[,.(value=var+dry+flood+heat+gsr)],
      data.table(expand.grid(var=c(NA,"v"),
                             flood=c(NA,"f"),
                             dry=c(NA,"d"),
                             heat=c(NA,"h"),
                             gsr=c(NA,"r")))[,.(label=paste(na.omit(c(var,flood,dry,heat,gsr)),collapse = "+")),by=.(var,flood,dry,heat,gsr)][,"label"]
    )
    
    haz_levels$label[haz_levels$value==0]<-"none"
    
    levels(haz_comb)<-haz_levels
    names(haz_comb)<-"combined hazards"
    
    plot(c(haz_comb,var,flood,heat_rf,dry_rf,gsr))
  
    # Improve palette
    color_map<-PalFun(PalName="turbo",N=nrow(haz_levels),Names = haz_levels$value)
    color_map<-data.frame(value=names(color_map),color=color_map)
    color_map$color[1]<-NA
    
    # Assign the colormap to the raster
    coltab(haz_comb) <- as.matrix(color_map)
    
    plot(haz_comb)  
    
    # Calc freq
    haz_comb_f<-data.table(freq(haz_comb))[order(count,decreasing = T)]
    haz_comb_f[,perc:=round(100*(count/sum(count)),1)]
    haz_comb_f
    
    # Aggregrate small percentages
    agg_names<-haz_comb_f[perc<2,value]
    agg_vals<-data.table(levels(haz_comb)[[1]])[label %in% agg_names,value]      
    haz_comb_agg<-classify(haz_comb,data.frame(from=agg_vals,to=999999))
    hca_levels<-data.table(levels(haz_comb)[[1]])
    hca_levels<-rbind(hca_levels,data.table(value=999999,label="other"))
    levels(haz_comb_agg)<-hca_levels
    names(haz_comb_agg)<-"combined hazards"
    
    # Create legend tables
    used_vals<-freq(haz_comb_agg)$value
    HazTab<-data.table(terra::levels(haz_comb_agg)[[1]])[label %in% used_vals]         
    levels(haz_comb_agg)<-HazTab
    
    # Improve palette
    color_map<-PalFun(PalName="turbo",N=nrow(hca_levels[label %in% used_vals]),Names = hca_levels[label %in% used_vals,value])
    color_map<-data.frame(value=names(color_map),color=color_map)
    color_map$color[1]<-NA
    
    # Assign the colormap to the raster
    coltab(haz_comb_agg) <- as.matrix(color_map)
    plot(haz_comb_agg)
    
    colnames(HazTab)<-c("Code","ShortName")
    HazTab<-HazTab[,Hazard:=ShortName
    ][ShortName=="v",Hazard:="Rain variability (v)"
    ][ShortName=="f",Hazard:="Flood (f)"
    ][ShortName=="d",Hazard:="Drought (d)"
    ][ShortName=="h",Hazard:="Heat (h)"
    ][ShortName=="r",Hazard:="Reduced season (r)"]
    
    used_vals<-freq(haz_comb)$value
    HazTab2<-data.table(terra::levels(haz_comb)[[1]])[label %in% used_vals]         
    levels(haz_comb)<-HazTab2
    colnames(HazTab2)<-c("Code","ShortName")
    HazTab2<-HazTab2[,Hazard:=ShortName
    ][ShortName=="v",Hazard:="Rain variability (v)"
    ][ShortName=="f",Hazard:="Flood (f)"
    ][ShortName=="d",Hazard:="Drought (d)"
    ][ShortName=="h",Hazard:="Heat (h)"
    ][ShortName=="r",Hazard:="Reduced season (r)"]

    
      #### Save results ####
    terra::writeRaster(haz_comb,filename = "raw_data/haz_comb/haz_full_rf.tif",overwrite=T)
    fwrite(HazTab2,"raw_data/haz_comb/haz_full_rf.csv")
    
    terra::writeRaster(haz_comb_agg,filename = "raw_data/haz_comb/haz_agg_rf.tif",overwrite=T)
    fwrite(HazTab,"raw_data/haz_comb/haz_agg_rf.csv")
    ### Irrigated ######
    haz_comb<-sum(c(flood,water_stress,heat_ir,gsr),na.rm=T)
    # Add labels
    haz_levels<-data.frame(
      data.table(expand.grid(flood=c(0,10),
                             waterstress=c(0,100),
                             heat=c(0,1000),
                             gsr=c(0,10000)))[,.(value=waterstress+flood+heat+gsr)],
      data.table(expand.grid(flood=c(NA,"f"),
                             waterstress=c(NA,"ws"),
                             heat=c(NA,"h"),
                             gsr=c(NA,"r")))[,.(label=paste(na.omit(c(flood,waterstress,heat,gsr)),collapse = "+")),by=.(flood,waterstress,heat,gsr)][,"label"]
    )
    
    haz_levels$label[haz_levels$value==0]<-"none"
    unique(values(haz_comb))
    
    levels(haz_comb)<-haz_levels
    names(haz_comb)<-"combined hazards"
    
    plot(c(haz_comb,flood,water_stress,heat_ir,gsr))
    
    # Improve palette
    color_map<-PalFun(PalName="turbo",N=nrow(haz_levels),Names = haz_levels$value)
    color_map<-data.frame(value=names(color_map),color=color_map)
    color_map$color[1]<-NA
    
    # Assign the colormap to the raster
    coltab(haz_comb) <- as.matrix(color_map)
    
    plot(haz_comb)  
    
    # Calc freq
    haz_comb_f<-data.table(freq(haz_comb))[order(count,decreasing = T)]
    haz_comb_f[,perc:=round(100*(count/sum(count)),1)]
    haz_comb_f
    
    # Aggregrate small percentages
    agg_names<-haz_comb_f[perc<2,value]
    agg_vals<-data.table(levels(haz_comb)[[1]])[label %in% agg_names,value]      
    haz_comb_agg<-classify(haz_comb,data.frame(from=agg_vals,to=999999))
    hca_levels<-data.table(levels(haz_comb)[[1]])
    hca_levels<-rbind(hca_levels,data.table(value=999999,label="other"))
    levels(haz_comb_agg)<-hca_levels
    names(haz_comb_agg)<-"combined hazards"
    
    # Create legend tables
    used_vals<-freq(haz_comb_agg)$value
    HazTab<-data.table(terra::levels(haz_comb_agg)[[1]])[label %in% used_vals]         
    levels(haz_comb_agg)<-HazTab
    
    # Improve palette
    color_map<-PalFun(PalName="turbo",N=nrow(hca_levels[label %in% used_vals]),Names = hca_levels[label %in% used_vals,value])
    color_map<-data.frame(value=names(color_map),color=color_map)
    color_map$color[1]<-NA
    
    # Assign the colormap to the raster
    coltab(haz_comb_agg) <- as.matrix(color_map)
    plot(haz_comb_agg)
    
    colnames(HazTab)<-c("Code","ShortName")
    HazTab<-HazTab[,Hazard:=ShortName
    ][ShortName=="r",Hazard:="Reduced season (r)"
    ][ShortName=="f",Hazard:="Flood (f)"
    ][ShortName=="ws",Hazard:="Water stress (ws)"
    ][ShortName=="h",Hazard:="Heat (h)"]
    
    used_vals<-freq(haz_comb)$value
    HazTab2<-data.table(terra::levels(haz_comb)[[1]])[label %in% used_vals]         
    levels(haz_comb)<-HazTab2
    colnames(HazTab2)<-c("Code","ShortName")
    HazTab2<-HazTab2[,Hazard:=ShortName
    ][ShortName=="r",Hazard:="Reduced season (r)"
    ][ShortName=="f",Hazard:="Flood (f)"
    ][ShortName=="ws",Hazard:="Water stress (ws)"
    ][ShortName=="h",Hazard:="Heat (h)"]
    
    # Save results ####
    terra::writeRaster(haz_comb,filename = "raw_data/haz_comb/haz_full_ir.tif",overwrite=T)
    fwrite(HazTab2,"raw_data/haz_comb/haz_full_ir.csv")
    
    terra::writeRaster(haz_comb_agg,filename = "raw_data/haz_comb/haz_agg_ir.tif",overwrite=T)
    fwrite(HazTab,"raw_data/haz_comb/haz_agg_ir.csv")
    
    