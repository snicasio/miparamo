library(raster)
library(rgeos)
library(rgdal)
library(geoR)
library(dplyr)
library(rgbif)
library(spThin)


    #   MOODELOS DE DISTRIBUCION

        ##  LECTURA DE DATOS Y DELIMITACION DEL AREA DE ANALISIS

            ### Lectura de coordenadas y descarga desde gbif

lw <- occ_data(scientificName="Leopardus wiedii",hasCoordinate=T)  #   Leopardus weidii de gbif
mr <- occ_data(scientificName="Mazama rufina",hasCoordinate=T)  #   Mazama rufina de gbif
to <- occ_data(scientificName="Tremarctos ornatus",hasCoordinate=T)  #   Tremarctos ornatus de gbif
    lw_coor <- as.data.frame(lw$data[,3:4]); names(lw_coor) <- c("Latitude","Longitude")
    lw_coor <- cbind(Species=rep("Leopardus weidii",nrow(lw_coor)),lw_coor)
    mr_coor <- as.data.frame(mr$data[,3:4]); names(mr_coor) <- c("Latitude","Longitude")
    mr_coor <- cbind(Species=rep("Mazama rufina",nrow(mr_coor)),mr_coor)
    to_coor <- as.data.frame(to$data[,3:4]); names(to_coor) <- c("Latitude","Longitude")
    to_coor <- cbind(Species=rep("Tremarctos ornatus",nrow(to_coor)),to_coor)

species_list <- list(Leopardus_wiedii=lw_coor,
                     Mazama_rufina=mr_coor,
                     Tremarctos_ornatus=to_coor)


            ### Eliminamos aquellas coordenadas duplicadas

duplicado <- list()
    for(i in seq(species_list))
    {
        duplicado[[i]] <- sapply(1:nrow(species_list[[i]]), function(x)
            if(duplicated(species_list[[i]][,2])[x]==TRUE & duplicated(species_list[[i]][,3])[x]==TRUE){
                T
            }else{
                F
            })
        species_list[[i]] <- species_list[[i]][duplicado[[i]]==F,]
        
                    ### Eliminamos los registros que tengan una distancia < 50 km entre si (para evitar sesgos en los modelos)
        
        thin(loc.data=species_list[[i]],
             lat.col="Latitude",
             long.col="Longitude",
             spec.col="Species",
             thin.par=50,reps=100,
             out.dir=getwd(),
             max.files=1,
             out.base=gsub(" ","_",names(species_list)[i]))
    }


            ### Limites del area de analisis

colombia <- getData("GADM",country="Colombia",level=1)
deptos <- bind(colombia[colombia$NAME_1 %in% c("Santander"),],
               colombia[colombia$NAME_1 %in% c("Norte de Santander"),])
deptos_cent <- gCentroid(colombia)


            ### Lectura de registros limpiados y conversion a shapefile

thinned_list <- lapply(list.files(pattern="thin1"),read.csv)
    names(thinned_list) <- gsub("_thin1.csv","",list.files(pattern="thin1"))

coord_list <- list()
    for(i in seq(thinned_list))
    {
        coord_list[[i]] <- SpatialPoints(thinned_list[[i]][,-1],proj4string=CRS("+proj=longlat"))
        coord_list[[i]] <- crop(coord_list[[i]],colombia)
    }
    names(coord_list) <- names(thinned_list)

par(mfrow=c(1,3))
    for(i in seq(coord_list))
    {
        plot(colombia)
        plot(coord_list[[i]],main=names(coord_list)[i],add=T)
    }


            ### Descarga de capas bioclimaticas y corte de la region de estudio

bioclim <- mask(crop(getData("worldclim",
                             var="bio",
                             res=0.5,
                             lon=as.data.frame(deptos_cent)[1,1],
                             lat=as.data.frame(deptos_cent)[1,2]),
                     colombia),
                colombia)


            ### Extraccion de valores climaticos por registro de cada especie

bio_data <- lapply(coord_list,extract,x=bioclim)
    for(i in seq(bio_data))
    {
        bio_data[[i]] <- as.data.frame(bio_data[[i]])
        names(bio_data[[i]]) <- gsub("_23","",names(bio_data[[i]]))
        bio_data[[i]] <- bio_data[[i]][!is.na(bio_data[[i]]$bio1)==T,]
    }


            ### Analisis de correlacion y de componentes principales para determiinar las variables con mayor importancia para la distribucion climatica

cor(bio_data[[1]][,-c(5:11,13,14,16:19)]) # Matriz de correlacion entre variables bioclimaticas (eliminamos aquellas con una correlacion > 0.7)
cor(bio_data[[2]][,-c(5:11,13,16:19)]) # Matriz de correlacion entre variables bioclimaticas (eliminamos aquellas con una correlacion > 0.7)
cor(bio_data[[3]][,-c(5:11,13,16:19)]) # Matriz de correlacion entre variables bioclimaticas (eliminamos aquellas con una correlacion > 0.7)


pca_lw <- prcomp(bio_data[[1]][,-c(5:11,13,14,16:19)]) #    Analisis de componentes principales
    pca_lw; summary(pca_lw); biplot(pca_lw)

pca_mr <- prcomp(bio_data[[2]][,-c(5:11,13,16:19)]) #    Analisis de componentes principales
    pca_mr; summary(pca_mr); biplot(pca_mr)

pca_to <- prcomp(bio_data[[3]][,-c(5:11,13,16:19)]) #    Analisis de componentes principales
    pca_to; summary(pca_to); biplot(pca_to)
    
list_pca <- list(pca_lw,pca_mr,pca_to); names(list_pca) <- names(bio_data)

bio_pca <- lapply(list_pca, function(x) as.numeric(gsub("bio","",rownames(as.data.frame(x$rotation))))) #   Esto me servira para guardar las capas de interes para los analisis de cada especie


            ### Guardamos los datos con el formato requerido para ejecutar el paquete "kuenm"

                ####    Puntos de entrenamiento, de prueba e independientes

independent <- read.csv("Especies.csv")
independent_list <- lapply(unique(independent$Species), function(x) independent %>% filter(Species==x))
    names(species_list) <- unique(independent$Species)


lista_csv <- rep(list(list()),length(coord_list))   #   Los archivos *.csv por especie
    names(lista_csv) <- names(coord_list)
    for(i in seq(coord_list))
    {
        lista_csv[[i]][[1]] <- data.frame(Species=rep(names(coord_list)[i],times=nrow(coord_list[[i]]@coords)),
                                          coord_list[[i]]@coords)
        lista_csv[[i]][[1]] <- lista_csv[[i]][[1]][as.numeric(rownames(bio_data[[i]])),]
        rownames(lista_csv[[i]][[1]]) <- seq(rownames(lista_csv[[i]][[1]]))
            
        lista_csv[[i]][[2]] <- lista_csv[[i]][[1]][sample(1:nrow(lista_csv[[i]][[1]]), size=round(nrow(lista_csv[[i]][[1]])*0.7), replace=F),]
        
        lista_csv[[i]][[3]] <- lista_csv[[i]][[1]][-as.numeric(rownames(lista_csv[[i]][[2]])),]
        
        lista_csv[[i]][[4]] <- independent_list[[i]]
        
        names(lista_csv[[i]]) <- c("Sp_joint","Sp_train","Sp_test","Sp_ind")
        
        for(j in seq(lista_csv[[i]]))
        {
            write.csv(lista_csv[[i]][[j]],paste(getwd(),names(lista_csv)[i],paste0(names(lista_csv[[i]])[j],".csv"),sep="/"),row.names=F)
        }
    }

            ### Funcion permutadora

#perm <- function(v) {
#    n <- length(v)
#    if (n == 1) v
#    else {
#        X <- NULL
#        for (i in 1:n) X <- rbind(X, cbind(v[i], perm(v[-i])))
#        X
#    }
#}

                ####    Archivos raster de las variables ambientales en formato ascii

lista_asc <- lapply(bio_pca, function(x) subset(bioclim,x))
    for(i in seq(lista_asc))
        for(j in seq(lista_asc[[i]]@layers))
        {
            writeRaster(lista_asc[[i]]@layers[[j]],
                        paste(getwd(),
                              names(lista_csv)[i],
                              "M_variables",
                              "Set_1",
                              paste0(gsub("_23","",names(lista_asc[[i]]))[j],
                                     ".asc"),
                              sep="/"),
                        format="ascii",
                        overwrite=T)
        }



        ##  ANALISIS MAXENT

            ### Installing and loading packages

#if(!require(devtools)){
#    install.packages("devtools")
#}
#
#if(!require(kuenm)){
#    devtools::install_github("marlonecobos/kuenm")
#}
    

library(kuenm)

            ### Leopardus wiedii

setwd(paste(getwd(),"Leopardus_wiedii",sep="/"))
dir()
file_name <- "lwie_enm_process"
kuenm_start(file.name = file_name)



setwd(paste(getwd(),"Mazama_rufina",sep="/"))
dir()
file_name <- "mruf_enm_process"
kuenm_start(file.name = file_name)



setwd(paste(getwd(),"Tremarctos_ornatus",sep="/"))
dir()
file_name <- "torn_enm_process"
kuenm_start(file.name = file_name)