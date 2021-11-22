
## load palm distributional data 

# get shapefilenames 
filenames <- list.files("DATA/VAR/Palms/Shapefiles/")[grep(".shp", list.files("DATA/VAR/Palms/Shapefiles/"))]
# prune .shp.xml filetypes
filenames <- filenames[-grep(".xml", filenames)]
# recreate path
filepath <- paste0("DATA/VAR/Palms/Shapefiles/", filenames)

# remove shapes with errors
filenames <-filenames[-grep("Geonoma_arundinacea.shp", filenames)]
filenames <-filenames[-grep("TEMPLATE_WGS84.shp", filenames)]
filepath <- paste0("DATA/VAR/Palms/Shapefiles/", filenames)

# Make a list of shapefiles
PalmShapeList <- lapply(filepath, sf::read_sf)
# standardize crs
crs(neotropic)<- crs( sf::as_Spatial(PalmShapeList[[1]]))
# query which neotropical provicnes match 
ProvDistPalm <- lapply(1:length(filenames), function(x) over(sf::as_Spatial(PalmShapeList[[x]]), neotropic)$Province_1)
names(ProvDistPalm) <- stringr::str_split(filenames, "\\.", simplify = T)[,1]
ProvDistPalm <- reshape2::melt(ProvDistPalm)


## load mammal distributional data 


## select species with some degree of frugivory
frugTraits <- animalTrait[animalTrait$Diet.Fruit > 0,]
frugTraits <- frugTraits[!c(names(frugTraits) %in% c("Diet.Source","Diet.Certainty", "BodyMass.Source", "BodyMass.SpecLevel"))]


## load species distribution data

mammShape <- shapefile("AnimalDistribution/Mammal/TERRESTRIAL_MAMMALS.shp")
# get the shapes corresponding with animal frugivores with we have trait information data 
mammShape <- mammShape[mammShape@data$binomial %in% frugTraits$FRUGIVORE,]
# subset those shapes so it only contains america 
# get the lower limit (lat)
Corxmax <- sapply(1:length(mammShape$binomial), function(x) 
  extent(bbox(subset(mammShape, 
                     mammShape$binomial == mammShape$binomial[[x]])))@xmax)
# subset only for american species
AmMammal <- lapply(which(Corxmax <0), function(x)
  subset(mammShape, 
         mammShape$binomial == mammShape$binomial[[x]]))



# query which neotropical provicnes match 
ProvDistMammal <- lapply(1:length(AmMammal), function(x) over(AmMammal[[x]], neotropic)$Province_1)
names(ProvDistMammal) <-sapply(1:length(AmMammal), function(x) unique(AmMammal[[x]]$binomial))
ProvDistMammal <- reshape2::melt(ProvDistMammal)

################
## Calculate richness and add to the neotropic shapefile

PalmsPerProvince <- rowSums(table(ProvDistPalm$value, ProvDistPalm$L1))
MammalsPerProvince <- rowSums(table(ProvDistMammal$value, ProvDistMammal$L1))

neotropic$PalmRichnes <- PalmsPerProvince[match(neotropic$Province_1, names(PalmsPerProvince))]
neotropic$MammalRichnes <- MammalsPerProvince[match(neotropic$Province_1, names(MammalsPerProvince))]

# Order and match
ProvDistPalm$L1 <- stringr::str_replace(ProvDistPalm$L1, "_", " ")
ProvDistPalm$SH <- palmTraitsRel$MaxStemHeight_m[match( ProvDistPalm$L1,palmTraitsRel$SpecName)]
ProvDistPalm$FL <- palmTraitsRel$AverageFruitLength_cm[match( ProvDistPalm$L1,palmTraitsRel$SpecName)]



# Order and match

ProvDistMammal$BM <- frugTraits$BodyMass.Value[match( ProvDistMammal$L1,frugTraits$FRUGIVORE)]


scSp <- scores(dietRDA)$sites
ProvDistMammal$H_I <- scSp[,1][match( ProvDistMammal$L1,rownames(scSp))]
ProvDistMammal$FD <- scSp[,2][match( ProvDistMammal$L1,rownames(scSp))]





################
## Calculate mean and sd variation per relevant trait and to the neotropic shapefile

# palms 

SDPalm <- aggregate(log1p(ProvDistPalm$FL), list(ProvDistPalm$value), sd, na.rm = T)
SDPalm$SH <- aggregate(log1p(ProvDistPalm$SH), list(ProvDistPalm$value), sd, na.rm = T)$x
SDPalm$FL_mean <- aggregate(log1p(ProvDistPalm$FL), list(ProvDistPalm$value), mean, na.rm = T)$x
SDPalm$SH_mean <- aggregate(log1p(ProvDistPalm$SH), list(ProvDistPalm$value), mean, na.rm = T)$x


# mammals
SDMamm <- aggregate(log1p(ProvDistMammal$BM), list(ProvDistMammal$value), sd, na.rm = T)
SDMamm$BM_mean <- aggregate(log1p(ProvDistMammal$BM), list(ProvDistMammal$value), mean, na.rm = T)$x
SDMamm$H_I_x <- aggregate((ProvDistMammal$H_I), list(ProvDistMammal$value), mean, na.rm = T)$x
SDMamm$FD_x <- aggregate((ProvDistMammal$FD), list(ProvDistMammal$value), mean, na.rm = T)$x

SDMamm$H_I_sd <- aggregate((ProvDistMammal$H_I), list(ProvDistMammal$value), sd, na.rm = T)$x
SDMamm$FD_sd <- aggregate((ProvDistMammal$FD), list(ProvDistMammal$value), sd, na.rm = T)$x


#### add trait data to neotropic shapefile 
# palms
neotropic$FL_SD <- SDPalm$x[match(neotropic$Province_1, SDPalm$Group.1)]
neotropic$SH_SD <- SDPalm$SH[match(neotropic$Province_1, SDPalm$Group.1)]
neotropic$SH_x <- SDPalm$SH_mean[match(neotropic$Province_1, SDPalm$Group.1)]
neotropic$FL_x <- SDPalm$FL_mean[match(neotropic$Province_1, SDPalm$Group.1)]
# mammals

neotropic$BM_SD <- SDMamm$x[match(neotropic$Province_1, SDMamm$Group.1)]
neotropic$BM_x <- SDMamm$BM_mean[match(neotropic$Province_1, SDMamm$Group.1)]

neotropic$HI_SD <- SDMamm$H_I_sd[match(neotropic$Province_1, SDMamm$Group.1)]
neotropic$HI_x <- SDMamm$H_I_x[match(neotropic$Province_1, SDMamm$Group.1)]
neotropic$FD_SD <- SDMamm$FD_sd[match(neotropic$Province_1, SDMamm$Group.1)]
neotropic$FD_x <- SDMamm$FD_x[match(neotropic$Province_1, SDMamm$Group.1)]


