if(!require("pacman", character.only = TRUE)){install.packages("pacman",dependencies = T)}

p_load(data.table,terra,install = T)

# Use the isciences version and not the CRAN version of exactextractr
if(!require("exactextractr")){
  remotes::install_github("isciences/exactextractr")
}

# Load data ####
cgiar_countries <- vect('Data/CGIAR_countries_simplified.shp')
cgiar_regions<-aggregate(cgiar_countries,by="CG_REG")

spam_africa<-terra::rast("Data/SPAM/ssa_crop_vop15_intd15.tif")
spam_global<-terra::rast("Data/SPAM/global_crop_vop15_int15.tif")

# CGIAR regions ####

## NOT SSA #####
Regions<-cgiar_regions[!cgiar_regions$CG_REG %in% c("ESA","WCA"),]
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

## SSA #####
Regions<-cgiar_regions[cgiar_regions$CG_REG %in% c("ESA","WCA"),]
SPAMbyRegion<-data.table(exactextractr::exact_extract(spam_africa,sf::st_as_sf(Regions),fun="sum",append_cols=c("CG_REG")))
setnames(SPAMbyRegion,"CG_REG","Region")
SPAMbyRegion<-data.table::melt(SPAMbyRegion,id.vars = "Region",value.name = "VoP",variable.name = "Crop")
SPAMbyRegion[,Crop:=gsub("sum.","",Crop)]
# Add rank
SPAMbyRegion<-SPAMbyRegion[order(Region,VoP,decreasing = T)][,Rank:=1:.N,by=Region][VoP==0,Rank:=NA]
# Add labels
SPAMbyRegion[,Label:=paste0(Region,"-","","-",Crop,"-",Rank)]

data$ssa_regions<-SPAMbyRegion


# CGIAR countries ####
## NOT SSA #####
Regions<-cgiar_countries[!cgiar_regions$CG_REG %in% c("ESA","WCA"),]
SPAMbyRegion<-data.table(exactextractr::exact_extract(spam_global,sf::st_as_sf(Regions),fun="sum",append_cols=c("CG_REG","ADMIN","ADM0_A3")))
setnames(SPAMbyRegion,c("CG_REG","ADMIN","ADM0_A3"),c("Region","Country","iso3"))
SPAMbyRegion<-data.table::melt(SPAMbyRegion,id.vars = c("Region","Country","iso3"),value.name = "VoP",variable.name = "Crop")
SPAMbyRegion[,Crop:=gsub("sum.","",Crop)]
# Add rank
SPAMbyRegion<-SPAMbyRegion[order(Region,Country,VoP,decreasing = T)][,Rank:=1:.N,by=.(Region,Country)][VoP==0,Rank:=NA]
# Add labels
SPAMbyRegion[,Label:=paste0(Region,"-",Country,"-",Crop,"-",Rank)]

data$global_countries<-SPAMbyRegion

## SSA #####
Regions<-cgiar_countries[cgiar_regions$CG_REG %in% c("ESA","WCA"),]
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




