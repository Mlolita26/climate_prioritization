if(!require("pacman", character.only = TRUE)){install.packages("pacman",dependencies = T)}

p_load(data.table,terra,wbstats,pbapply,install = T)

# Use the isciences version and not the CRAN version of exactextractr
if(!require("exactextractr")){
  remotes::install_github("isciences/exactextractr")
}

# 1) Prepare Boundaries ####
cgiar_countries <- terra::vect('raw_data/CGIAR_countries_simplified.geojson')
region_file<-"raw_data/regions.geojson"
country_file<-"raw_data/countries.geojson"

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

  ## 1.1) Set Base raster #####
  base_rast<-terra::rast("raw_data/haz_comb/haz_agg_rf.tif")
  ## 1.2) Resample SPAM #####
  spam_techs<-c("r","i","a")
  spam_ssa_path<-"raw_data/SPAM/ssa_crop_vop15_intd15_"
  spam_global_path<-"raw_data/SPAM/global_crop_vop15_int15_"
  spam_files<-apply(expand.grid(c(spam_ssa_path,spam_global_path),spam_techs,".tif"),1,FUN=paste,collapse="")

  overwrite<-F
  for(i in 1:length(spam_files)){
    cat(i,"/",length(spam_files),"\n")
    file_new<-gsub(".tif","_rs.tif",spam_files[i])
    if(!file.exists(file_new)|overwrite==T){
    spam_dat<-terra::rast(spam_files[i])
    
    cs<-cellSize(spam_dat,unit="km")
    spam_dat_d<-spam_dat/cs
    if(grepl("global",file_new)){
      spam_dat_d<-terra::resample(spam_dat_d,base_rast,method="near")
    }else{
      spam_dat_d<-terra::resample(spam_dat_d,base_rast,method="bilinear")
    }
    spam_dat<-spam_dat_d*cellSize(spam_dat_d,unit="km")
    
    terra::writeRaster(spam_dat,filename=file_new,overwrite=T)
    }
  }
  
  ## 1.3) Rasterize geographies #####
  Regions_rast<-lapply(1:length(cgiar_regions),FUN=function(i){terra::rasterize(cgiar_regions[i],base_rast,field="CG_REG")})
  names(Regions_rast)<-cgiar_regions$CG_REG
  Countries_rast<-terra::rasterize(cgiar_countries,base_rast,field="ADMIN")
  
  countries<-data.table(levels(Countries_rast)[[1]])
  country_list<-pblapply(1:nrow(countries),FUN=function(i){
    a<-Countries_rast
    value<-countries[i,ID]
    a[values(a)!=value]<-NA
    return(a)
  })
  names(country_list)<-countries$ADMIN
  
  ## 1.4) Combine spam global and spam africa #####
  overwrite<-F
  spam_comb_path<-gsub("global","comb",spam_global_path)
  for(tech in spam_techs){
    save_path<-paste0(spam_comb_path,tech,"_rs.tif")
    if(!file.exists(save_path)|overwrite){
    # Load mapspam data
    spam_africa<-terra::rast(paste0(spam_ssa_path,tech,"_rs.tif"))+0
    spam_global<-terra::rast(paste0(spam_global_path,tech,"_rs.tif"))+0
    if("BANA" %in% names(spam_global)){
      names(spam_global)<-ms_codes[match(names(spam_global),Code),Fullname]
    }
    
    crops<-names(spam_africa)
    crops<-crops[crops %in% names(spam_global)]
    
    spam_combined <-rast(lapply(crops, FUN = function(crop) {
        cat(crop,"\n")
        a <- spam_global[[crop]]
        b <- spam_africa[[crop]]
        a[!is.na(values(b))] <- values(b)[!is.na(values(b))]
        a
      }))

    writeRaster(spam_combined,save_path)
    }
  }
    
  if(F){
    # Check results

    a<-terra::rast(paste0(spam_comb_path,"a_rs.tif"))
    i<-terra::rast(paste0(spam_comb_path,"i_rs.tif"))
    r<-terra::rast(paste0(spam_comb_path,"r_rs.tif"))
    
    crop<-"arabica coffee"
    crop<-"soybean"
    
    a<-a[crop]
    i<-i[crop]
    r<-r[crop]
    
    # i + r should be virtually the same as a
    plot(a-sum(c(i,r),na.rm=T))
  
    admin<-"LAC"
    v<-cgiar_regions[cgiar_regions$CG_REG==admin,]
    a<-crop(mask(a,v),v)
    i<-crop(mask(i,v),v)
    r<-crop(mask(r,v),v)
    plot(c(a,sum(c(i,r),na.rm=T)))
  }

# 2) Extract VoP ####
for(tech in spam_techs){
  ## Load mapspam data #####
  spam_combined<-terra::rast(paste0(spam_comb_path,tech,"_rs.tif"))

  # Regions
  region_spam<-data.table(rbindlist(pblapply(1:length(Regions_rast),FUN=function(k){
    zonal(spam_combined,Regions_rast[[k]],fun="sum",na.rm=T)
  })))
  setnames(region_spam,"CG_REG","Region")
  region_spam<-data.table::melt(region_spam,id.vars = "Region",value.name = "VoP",variable.name = "Crop")

  # Countries
  country_spam<-data.table(zonal(spam_combined,Countries_rast,fun="sum",na.rm=T))
  setnames(country_spam,c("ADMIN"),c("Country"))
  country_spam<-data.table::melt(country_spam,id.vars = c("Country"),value.name = "VoP",variable.name = "Crop")
  country_spam<-merge(country_spam,data.frame(cgiar_countries),by.x="Country",by.y="ADMIN",all.x=T)
  # Save results
  spam<-rbindlist(list(region_spam,country_spam),fill=T)
  fwrite(spam,file=file.path("raw_data/SPAM",paste0("SPAMextracted_",tech,".csv")))
}

## 2.1) Check results #####
spam_i<-fread("raw_data/SPAM/SPAMextracted_i.csv")
spam_r<-fread("raw_data/SPAM/SPAMextracted_r.csv")
spam_a<-fread("raw_data/SPAM/SPAMextracted_a.csv")

crop<-"coffee"
crop<-"soybean"
crop<-"cassava"

admin<-"Africa"
spam_i[Region==admin & Country=="" & grepl(crop,Crop),sum(VoP)] +
spam_r[Region==admin & Country=="" & grepl(crop,Crop),sum(VoP)]
spam_a[Region==admin & Country=="" & grepl(crop,Crop),sum(VoP)]

admin<-"LAC"
spam_i[Region==admin & Country=="" & grepl(crop,Crop),sum(VoP)] +
spam_r[Region==admin & Country=="" & grepl(crop,Crop),sum(VoP)]
spam_a[Region==admin & Country=="" & grepl(crop,Crop),sum(VoP)]


# 3) Prepare Hazards Data ####
  hazard_layers<-data.table(file=c("raw_data/haz_comb/haz_full_rf.tif",
                                   "raw_data/haz_comb/haz_full_ir.tif"),
                            legend=c('raw_data/haz_comb/haz_full_rf.csv',
                                     'raw_data/haz_comb/haz_full_ir.csv'))
  hazard_layers[,file_n:=c("rf","ir")][,irrigated:=c(F,T)]
  
  
  for(i in 1:nrow(hazard_layers)){
    cat(i,"/",nrow(hazard_layers))
    haz_choice<-hazard_layers$file_n[i]
    irrigated<-hazard_layers[i,irrigated]
  
  haz_vop_file<-paste0("raw_data/haz_comb/hazard_",haz_choice,"_vop_admin.parquet")
  hazard_file<-hazard_layers[file_n==haz_choice,file]
  HazTab<-fread(hazard_layers[file_n==haz_choice,legend])
  
  ## Load hazard layer #####
  Hazard<- terra::rast(hazard_file)
  Hazard[is.na(Hazard[])]<-0
  HazardCrop<-terra::rast(hazard_file)
  
  ## Load mapspam #####
  if(irrigated==F){
    spam_file<-paste0(spam_comb_path,"r_rs.tif")
    SPAMbyRegion<-fread("raw_data/SPAM/SPAMextracted_r.csv")
  }else{
    spam_file<-paste0(spam_comb_path,"i_rs.tif")
    SPAMbyRegion<-fread("raw_data/SPAM/SPAMextracted_i.csv")
  }
  
  spam<-terra::rast(spam_file)

  ## extract hazard x vop by admin area ######

  spam_regions<-rbind(cgiar_countries,cgiar_regions)

  Data<-rbindlist(pbapply::pblapply(1:length(spam_regions),FUN=function(j){
    Region<-spam_regions[j]
    
    SPAM <- terra::mask(terra::crop(spam,Region),Region)
    Haz<-terra::mask(terra::crop(Hazard,Region),Region)
    
    dff<-rbindlist(lapply(1:nlyr(SPAM),FUN=function(k){
      # Take one crop
      cr1 <- SPAM[[k]]
      
      ex<-data.table(VoP=cr1[],Code=as.numeric(values(Haz)))
      colnames(ex)[1]<-"VoP"
      ex<-merge(ex,HazTab,by="Code",all.x=T)
      ex[is.na(Hazard),Hazard:="none"]
      ex<-ex[,list(VoP=sum(VoP,na.rm = T)),by=list(Hazard)
      ][,Crop:=names(cr1)]
      ex
    }))
    
    dff[,Region:=data.frame(Region)$CG_REG]
    cntry<-as.character(Region$ADMIN)
    dff[,Country:=cntry]
    
    dff
  }))
  
  Data[,VoP_total_crop:=sum(VoP,na.rm=T),by=.(Region,Country,Crop)
  ][,VoP_total:=sum(VoP,na.rm=T),by=.(Region,Country)
  ][,VoP_perc_total:=round(100*VoP/VoP_total,2)
  ][,VoP:=round(VoP,0)
  ][,VoP_total_crop:=round(VoP_total_crop,0)
  ][,VoP_perc_total:=round(VoP_perc_total,0)]
  
  Data[,Admin:=Country][Admin=="",Admin:=Region]
  
  arrow::write_parquet(Data,haz_vop_file)
  }
    
# 4) Combine Rainfed and Irrigated ####
  files<-list.files("raw_data/haz_comb","vop_admin",full.names = T)
  
  # All hazards rainfed
  a<-arrow::read_parquet(files[1])
  # All hazards irrigated
  b<-arrow::read_parquet(files[2])
  
  # Check data
  crop<-"rice"
  region<-"SEA"
  a[Region==region & is.na(Country) & Crop==crop,sum(VoP)]
  b[Region==region & is.na(Country) & Crop==crop,sum(VoP)]
  
  region<-"SA"
  a[Region==region & is.na(Country) & Crop==crop,sum(VoP)]
  b[Region==region & is.na(Country) & Crop==crop,sum(VoP)]
  
  ab<-rbind(a,b)[order(Region,Country,Crop,Hazard)]
  ab[Hazard=="No hazard",Hazard:="none"]
  ab<-ab[,.(VoP=sum(VoP)),by=.(Region,Country,Crop,Hazard,Admin)]
  ab[,VoP_tot_crop:=sum(VoP),by=.(Region,Country,Crop)]
  ab[,VoP_tot:=sum(VoP),by=.(Region,Country)]
  
  ab[,Admin:=Region][!is.na(Country),Admin:=Country]

  arrow::write_parquet(ab,"Data/hazard_vop_admin.parquet")
  
  

    