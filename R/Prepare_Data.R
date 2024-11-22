if(!require("pacman", character.only = TRUE)){install.packages("pacman",dependencies = T)}

p_load(data.table,terra,wbstats,install = T)

# Use the isciences version and not the CRAN version of exactextractr
if(!require("exactextractr")){
  remotes::install_github("isciences/exactextractr")
}

# 1) Prepare Boundaries ####
cgiar_countries <- terra::vect('raw_data/CGIAR_countries_simplified.geojson')
region_file<-"Data/regions.geojson"
country_file<-"Data/countries.geojson"

# Remove high income countries from regions
# Retrieve country information
country_info <- data.table(wb_countries())[,.(iso3c,income_level,region)][,LMIC:=T
][income_level=="High income",LMIC:=F
][income_level %in% c("Not Classified","Aggregates"),LMIC:=NA]

# Merge income level with country vector
cgiar_countries<-merge(cgiar_countries,country_info,by.x="ADM0_A3",by.y="iso3c")
cgiar_regions<-terra::aggregate(cgiar_countries[cgiar_countries$LMIC==T],by="CG_REG")
cgiar_regions[,c("region","income_level","agg_n","LMIC")]<-NULL

# Create larger regions and merge back with cgiar_regions

middle_east_countries<-c("IRN","IRQ","JOR","KWT","LBN","OMN","SAU","SYR","ARE","YEM","PSE","BHR","QAT","ISR","MLT")
cgiar_regions2<-cgiar_countries
cgiar_regions2$region2<-NA
cgiar_regions2$region2[grepl("Asia",cgiar_regions2$region)]<-"Asia"
cgiar_regions2$region2[grepl("Africa",cgiar_regions2$region) & !cgiar_regions2$ADM0_A3 %in% middle_east_countries]<-"Africa"
cgiar_regions2$region2[grepl("Latin",cgiar_regions2$region)]<-"Latin America"
cgiar_regions2<-terra::aggregate(cgiar_regions2[!is.na(cgiar_regions2$region2) & cgiar_regions2$LMIC==T],by="region2")
cgiar_regions2$CG_REG<-cgiar_regions2$region2
cgiar_regions2[,c("region","region2","income_level","agg_n","LMIC")]<-NULL
cgiar_regions<-rbind(cgiar_regions,cgiar_regions2)

cgiar_regions2$region<-cgiar_regions2$CG_REG
cgiar_regions2[,c("CG_REG","ADM0_A3","ADMIN")]<-NULL
cgiar_countries$region<-NULL

region_info <- terra::extract(cgiar_regions2, cgiar_countries)
cgiar_countries$region <- region_info$region

# Save processed vector data
terra::writeVector(cgiar_regions,region_file,overwrite=T)
terra::writeVector(cgiar_countries,country_file,overwrite=T)

# Load mapspam data
spam_africa<-terra::rast("raw_data/SPAM/ssa_crop_vop15_intd15.tif")
spam_global<-terra::rast("raw_data/SPAM/global_crop_vop15_int15.tif")

# 2) PrepareCGIAR regions ####

# 3) Extract VoP
  ## Regions
    ### NOT SSA #####
  Regions<-cgiar_regions[!cgiar_regions$CG_REG %in% c("ESA","WCA","Africa"),]
  SPAMbyRegion<-data.table(exactextractr::exact_extract(spam_global,sf::st_as_sf(Regions),fun="sum",append_cols=c("CG_REG")))
  setnames(SPAMbyRegion,"CG_REG","Region")
  SPAMbyRegion<-data.table::melt(SPAMbyRegion,id.vars = "Region",value.name = "VoP",variable.name = "Crop")
  SPAMbyRegion[,Crop:=gsub("sum.","",Crop)]
  # Add rank
  SPAMbyRegion<-SPAMbyRegion[order(Region,VoP,decreasing = T)][,Rank:=1:.N,by=Region][VoP==0,Rank:=NA]
  # Add labels
  SPAMbyRegion[,Label:=paste0(Region,"-","","-",Crop,"-",Rank)]
  
  data<-list()
  data$global_regions<-SPAMbyRegion
  
    ### SSA #####
  Regions<-cgiar_regions[cgiar_regions$CG_REG %in% c("ESA","WCA","Africa"),]
  SPAMbyRegion<-data.table(exactextractr::exact_extract(spam_africa,sf::st_as_sf(Regions),fun="sum",append_cols=c("CG_REG")))
  setnames(SPAMbyRegion,"CG_REG","Region")
  SPAMbyRegion<-data.table::melt(SPAMbyRegion,id.vars = "Region",value.name = "VoP",variable.name = "Crop")
  SPAMbyRegion[,Crop:=gsub("sum.","",Crop)]
  # Add rank
  SPAMbyRegion<-SPAMbyRegion[order(Region,VoP,decreasing = T)][,Rank:=1:.N,by=Region][VoP==0,Rank:=NA]
  # Add labels
  SPAMbyRegion[,Label:=paste0(Region,"-","","-",Crop,"-",Rank)]
  
  data$ssa_regions<-SPAMbyRegion
  
  
  ## CGIAR countries ####
    ### NOT SSA #####
  Regions<-cgiar_countries[!cgiar_countries$CG_REG %in% c("ESA","WCA"),]
  SPAMbyRegion<-data.table(exactextractr::exact_extract(spam_global,sf::st_as_sf(Regions),fun="sum",append_cols=c("CG_REG","ADMIN","ADM0_A3")))
  setnames(SPAMbyRegion,c("CG_REG","ADMIN","ADM0_A3"),c("Region","Country","iso3"))
  SPAMbyRegion<-data.table::melt(SPAMbyRegion,id.vars = c("Region","Country","iso3"),value.name = "VoP",variable.name = "Crop")
  SPAMbyRegion[,Crop:=gsub("sum.","",Crop)]
  # Add rank
  SPAMbyRegion<-SPAMbyRegion[order(Region,Country,VoP,decreasing = T)][,Rank:=1:.N,by=.(Region,Country)][VoP==0,Rank:=NA]
  # Add labels
  SPAMbyRegion[,Label:=paste0(Region,"-",Country,"-",Crop,"-",Rank)]
  
  data$global_countries<-SPAMbyRegion
  
    ### SSA #####
  Regions<-cgiar_countries[cgiar_countries$CG_REG %in% c("ESA","WCA"),]
  SPAMbyRegion<-data.table(exactextractr::exact_extract(spam_africa,sf::st_as_sf(Regions),fun="sum",append_cols=c("CG_REG","ADMIN","ADM0_A3")))
  setnames(SPAMbyRegion,c("CG_REG","ADMIN","ADM0_A3"),c("Region","Country","iso3"))
  SPAMbyRegion<-data.table::melt(SPAMbyRegion,id.vars = c("Region","Country","iso3"),value.name = "VoP",variable.name = "Crop")
  SPAMbyRegion[,Crop:=gsub("sum.","",Crop)]
  # Add rank
  SPAMbyRegion<-SPAMbyRegion[order(Region,Country,VoP,decreasing = T)][,Rank:=1:.N,by=.(Region,Country)][VoP==0,Rank:=NA]
  # Add labels
  SPAMbyRegion[,Label:=paste0(Region,"-",Country,"-",Crop,"-",Rank)]
  
  data$ssa_countries<-SPAMbyRegion
  
  data_merge<-rbindlist(data,fill=T)[,version:=2020]
  
  fwrite(data_merge,file=file.path("Data","SPAMextracted.csv"))
  
  
  
  
  
# 4) Prepare Hazards Data
  Haz_Ext_File<-"Data/Haz_Ext.parquet"
  ## Load hazard layer #####
  hazard_file<-'Data/combi_haz_table_cg.asc/combi_haz_table_cg.asc'
  Hazard<- terra::rast(hazard_file)
  Hazard[Hazard$combi_haz_table_cg==0] <-NA
  names(Hazard) <-'hazard'
  
  HazTab<-data.table(Code=1:9,
                     Hazard=c('Drought (D)','Flood (F)','Climate variability (V)','D + V', 'Growing season reduction (R)', 'High growing season temperature (T)','F + T','V + T','Other combination'),
                     ShortName=c('D','F','V','D+V', 'R', 'T','F+T','V+T','Other'))
  
  HazardCrop<-terra::rast(hazard_file)
  
  ## load global crop land extent #####
  #   Cropland extent data comes from https://glad.umd.edu/dataset/croplands
    Cropland<-terra::rast("raw_data/Global_cropland_3km_2019.tif")
    crs(HazardCrop)<-crs(Cropland)
    
    # Resample to hazard layer
    Cropland<-terra::resample(Cropland,HazardCrop)
    
    # Convert % to km2
    CroplandSize<-terra::cellSize(Cropland,unit="km")
    Cropland<-(Cropland/100)*CroplandSize
  
  
  ## Load mapspam #####
  spam_ssa_file<-"raw_data/SPAM/ssa_crop_vop15_intd15.tif"
  spam_global_file<-"raw_data/SPAM/global_crop_vop15_int15.tif"
  spam_ssa_file_rs<-"Data/SPAM/ssa_crop_vop15_intd15_rs.tif"
  spam_global_file_rs<-"Data/SPAM/global_crop_vop15_int15_rs.tif"
  
  spam_africa<-terra::rast(spam_ssa_file)
  spam_global<-terra::rast(spam_global_file)
  
  # Resmaple to hazards
  cs<-cellSize(spam_africa,unit="km")
  spam_africa_d<-spam_africa/cs
  spam_africa_d<-terra::resample(spam_africa_d,Hazard)
  spam_africa<-spam_africa_d*cellSize(spam_africa_d,unit="km")
  
  terra::writeRaster(spam_africa,filename=spam_ssa_file_rs,overwrite=T)
  
  # Resmaple to hazards
  cs<-cellSize(spam_global,unit="km")
  spam_global_d<-spam_global/cs
  spam_global_d<-terra::resample(spam_global_d,Hazard)
  spam_global<-spam_global_d*cellSize(spam_global_d,unit="km")
  
  terra::writeRaster(spam_global,filename=spam_global_file_rs,overwrite=T)
  
  ## extract hazard x vop by admin area ######
  haz_vop_file<-"Data/hazard_vop_admin.parquet"
  RegxCrop<-rbindlist(lapply(strsplit(SPAMbyRegion$Label,"-"),FUN=function(X){data.table(Region=X[1],Country=X[2],Crop=X[3],Rank=X[4])}))
  areas<-unique(RegxCrop[,.(Region,Country)])
  
  Data<-rbindlist(pbapply::pblapply(1:nrow(areas),FUN=function(i){
    REG<-areas$Region[i]
    country<-areas$Country[i]
    
    Crops<-RegxCrop[Region==REG & Country==country,as.character(Crop)]
    # Load crops for selected region
    if(REG %in% c("ESA","WCA")){
      SPAM<-spam_africa
    }else{
      SPAM<-spam_global
    }
    
    # Hazard descriptions table
    HazTab<-data.table(Code=1:9,Hazard=c('Drought (D)','Flood (F)','Climate variability (V)','D + V', 'Growing season reduction (R)', 'High growing season temperature (T)','F + T','V + T','Other combination'))
    
    # Crop & hazard by selected CGIAR region
    if(country!=""){
      Region<-cgiar_countries[cgiar_countries$CG_REG==REG & cgiar_countries$ADMIN==country,]
    }else{
      Region<-cgiar_regions[cgiar_regions$CG_REG==REG,]
    }
    
    SPAM <- terra::mask(terra::crop(SPAM,Region),Region)
    Haz<-terra::mask(terra::crop(Hazard,Region),Region)
    
    dff<-rbindlist(lapply(1:nlyr(SPAM),FUN=function(i){
      # Take one crop
      cr1 <- SPAM[[i]]
      
      ex<-data.table(VoP=as.numeric(values(cr1)),Code=as.numeric(values(Haz)))
      ex<-merge(ex,HazTab,by="Code",all.x=T)
      ex<-ex[,list(VoP=sum(VoP,na.rm = T)),by=list(Hazard)
      ][,Crop:=names(cr1)]
      ex
    }))
    
    dff[,Region:=REG][,Country:=country]
    
    dff
  }))
  
  Data<-merge(Data,RegxCrop,by=c("Region","Country","Crop"),all.x=T,sort=F)
  Data[,VoP_total_crop:=sum(VoP,na.rm=T),by=.(Region,Country,Crop)
  ][,VoP_total:=sum(VoP,na.rm=T),by=.(Region,Country)
  ][,VoP_perc_total:=round(100*VoP/VoP_total,2)
  ][,Label:=paste0(Region,"-",Country,"-",Crop,"-",Rank)
  ][,VoP:=round(VoP,0)
  ][,VoP_total_crop:=round(VoP_total_crop,0)
  ][,VoP_perc_total:=round(VoP_perc_total,0)]
  
  
  Data[is.na(Hazard),Hazard:="No hazard"]
  
  # Remove rows with zero value
  Data<-Data[VoP>0]
  # Add back label
  Data[,Label:=paste0(Region,"-",Country,"-",Crop,"-",Rank)][,Admin:=Country][Admin=="",Admin:=Region]
  
  arrow::write_parquet(Data,haz_vop_file)
  
  ## load population #####
  Population<-terra::rast("raw_data/cell5m_ruralpop_2020_v3.tif")
  TotalPop<-terra::resample(terra::crop(Population,HazardCrop),HazardCrop)
  names(TotalPop)<-"TotalPop"
  MaskedPop<-terra::mask(TotalPop,HazardCrop)
  names(MaskedPop)<-"AtRiskPop"

  ## Extract cropland hazards by country #####
      EX<-data.table(terra::extract(c(HazardCrop,Cropland,TotalPop,MaskedPop),cgiar_countries))
      
      EX<-merge(EX,data.table(ID=1:length(cgiar_countries),Country=cgiar_countries$ADMIN),by="ID",all.x=T)
      
      # Hazard descriptions table
      setnames(EX,"combi_haz_table_cg","Code")
      
      EX<-merge(EX,HazTab[,list(Code,ShortName)],by="Code",all.x=T)[,ID:=NULL][,Code:=NULL]
      setnames(EX,"ShortName","Hazard")
      
      EX[is.na(Hazard),Hazard:="No Hazard"]
      
      EX<-EX[,list(Cropland.Area=round(sum(Global_cropland_3km_2019)),
                   RuralPop=round(sum(TotalPop,na.rm=T)),
                   RuralPopRisk=round(sum(AtRiskPop,na.rm = T))),by=list(Country,Hazard)
      ][,Total.Area:=sum(Cropland.Area),by=Country
      ][,RuralPop:=sum(RuralPop),by=Country
      ][,RuralPopRisk:=sum(RuralPopRisk),by=Country
      ][,Cropland.Perc:=round(100*Cropland.Area/Total.Area,1)
      ][!is.na(Hazard)
      ][,Cropland.Area:=NULL
      ][,Total.Perc:=sum(Cropland.Perc),by=Country
      ][,RuralPopRiskPerc:=round(100*RuralPopRisk/RuralPop,1)
      ][,RuralPop:=round(RuralPop/10^6,1)
      ][,RuralPopRisk:=round(RuralPopRisk/10^6,1)]
      
      EX<-dcast(EX,Country+Total.Area+Total.Perc+RuralPop+RuralPopRiskPerc+RuralPopRisk~Hazard,value.var = "Cropland.Perc")[Total.Area!=0]
      
      arrow::write_parquet(EX,Haz_Ext_File)