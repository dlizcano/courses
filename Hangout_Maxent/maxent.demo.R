#Instalación paquetes
#Estas lineas se deben correr si es la primera instalación de Maxent
#install.packages("dismo") #Añadir maxent.jar en la carpeta java de dismo
#install.packages("raster")
#install.packages("rgdal")
#install.packages("rJava")

#Configuración inicial
options(java.parameters = "-Xmx1g" )
library(raster)
library(dismo)
setwd("C:/Users/Jorge/Dropbox/Presentaciones/Hangout_Maxent")

#Cargar datos
occs <- read.csv("./samples/bradypus.csv") #Datos de presencia
View(occs)
layers <- stack(list.files("./layers","*.asc$",full.names=T))
plot(layers)
plot(layers[[1]])
layers #Ver en consola las caracteristicas del stack
layers[[1]] #Ver en consola las características de un stack particular
points(occs[,2:3],pch=18,cex=0.6)

#Run model
me <- maxent(layers, occs[,2:3], factors='ecoreg', removeDuplicates=TRUE, 
             path="./outputR")

#Explorar directorio y el html, comparar con el obtenido mediante la consola

#Formato SWD (samples with data)
pres.covs<-extract(layers, occs[,2:3])
pres.covs<-na.omit(pres.covs)
pres.covs<-unique(pres.covs)
View(pres.covs)

#Una forma más "correcta" de eliminar duplicados de la misma celda
pres.covs<-extract(layers, occs[,2:3],cellnumbers=T)
View(pres.covs)
pres.covs<-na.omit(pres.covs)
pres.covs<-unique(pres.covs)
pres.covs<-pres.covs[,-1] #Elimina columna de número de celda

bkg.covs<-sampleRandom(layers,10000,cells=T)
bkg.covs<-unique(bkg.covs)
bkg.covs<-bkg.covs[,-1] #Elimina columna de número de celda

#Condensar todo en una sola tabla
env.values<-data.frame(rbind(pres.covs,bkg.covs)) #rbind une tablas verticalmente
env.values$ecoreg<-as.factor(env.values$ecoreg) #Importante cuando se usan variables categoricas

#Etiquetar presencias (1) y background (0)
y <- c(rep(1,nrow(pres.covs)), rep(0,nrow(bkg.covs)))
me <- maxent(env.values, y, args=c("addsamplestobackground=true"), 
             path="./outputR")

#Visualizar los resultados
map <- predict(me, layers, progress="text")
plot(map)

#Guardar los resultados
writeRaster(map,"./outputR/map.tif")
save(me,file="./outputR/mx_obj.RData")

#Evaluación usando kfold partitioning
#Ejemplo para 1 fold
fold <- kfold(pres.covs, k=5) #Genera un indice aleatorio de los folds
occtest <- pres.covs[fold == 1, ]
occtrain <- pres.covs[fold != 1, ]

env.values<-data.frame(rbind(pres.covs, bkg.covs))
env.values$ecoreg<-as.factor(env.values$ecoreg) #Importante cuando se usan variables categoricas

me <- maxent(env.values, y, args=c("addsamplestobackground=true"), path="./outputR")
e <- evaluate(me, p=data.frame(occtest), a=data.frame(bkg.covs))
str(e)
?ModelEvaluation-class

#Threshold value that maximizes Kappa
plot(e@t,e@kappa,type="l")
e@t[which.max(e@kappa)]

#Computing True Skill Statistic = TPR(Sensitivity)+TNR(Specificity)-1
tss <- e@TPR+e@TNR-1
plot(e@t,tss,type="l")
e@t[which.max(tss)]

#AUC Plot: X=1-Specificity, Y=Sensitivity
plot((1-e@TNR),e@TPR,type="l",xlab="Fractional Predicted Area",
     ylab="Sensitiviy")
e@auc

#Algunas lecturas sobre estadísticos de desempeño y validación de modelos
#Fielding, A.H. & Bell, J.F. (1997) A review of methods for the assessment of
#prediction errors in conservation presence/absence models. 
#Environmental Conservation, 24, 38–49.

#ALLOUCHE, O., TSOAR, A. and KADMON, R. (2006), Assessing the accuracy of
#species distribution models: prevalence, kappa and the true skill statistic
#(TSS). Journal of Applied Ecology, 43: 1223–1232. 


#Now, for all folds
auc<-rep(NA,5)
max.tss<-rep(NA,5)
for (i in 1:5){
  occtest <- pres.covs[fold == i, ]
  occtrain <- pres.covs[fold != i, ]
  env.values<-data.frame(rbind(pres.covs, bkg.covs))
  env.values$ecoreg<-as.factor(env.values$ecoreg) #Importante cuando se usan variables categoricas
  me <- maxent(env.values, y, args=c("addsamplestobackground=true"), path="./outputR")
  e<-evaluate(me, p=data.frame(occtest), a=data.frame(bkg.covs))
  auc[i]<-e@auc
  lines((1-e@TNR),e@TPR)
  tss<-e@TPR+e@TNR-1
  max.tss[i]<-e@t[which.max(tss)]
}
mean(auc)
sd(auc)
mean(max.tss)

#Gráficos de curvas de respuesta. Recuerde que me era el nombre del objeto
#maxent con el modelo con todas las presencias el cual fue sobreescrito en
#el anterior loop. Sin embargo, lo podemos restaurar ya que lo guardamos.
load("./outputR/mx_obj.RData")
response(me)
response(me,var=1) #Usando indicador de columna
response(me, var="ecoreg") #Usando nombre de la variable

#Proyección modelos a otras épocas
hotlayers <- stack(list.files("./hotlayers", "*.asc$", full.names=T))
map.future <- predict(me, hotlayers, progress="text")

par(mfrow=c(1,2)) #Dividir ventana de grafico para comparar modelos
plot(map)
plot(map.future)

#Comparar visualmente modelos aplicando un threshold
umbral.tss <- mean(max.tss)
plot(map >= umbral.tss)
plot(map.future >= umbral.tss)

map.tss <- (map>=umbral.tss) #Crear objeto raster con umbral
map.future.tss <- (map.future>=umbral.tss) #Crear objeto raster con umbral

#Calcular áreas
area.layers <- area(layers[[1]]) #Capas usadas son de 0.05 grados que corresponden en el ecuador a 6 km aprox
plot(area.layers)
plot(map.tss*area.layers)
cellStats(map.tss, sum)*(6*6) #Area presente en km2 sin corregir proyección
cellStats(map.tss * area.layers, sum) #Area presente en km2 corregida
cellStats(map.future.tss * area.layers, sum) #Area futuro en km2

#Configuración "personalizada" de Maxent
response(me)
mxnt.args=c("autofeature=FALSE",
            "linear=TRUE",
            "quadratic=TRUE",
            "product=FALSE",
            "hinge=FALSE",
            "threshold=FALSE") #Seleccionar features manualmente
me.mfeatures <- maxent(env.values, y, args=mxnt.args, path="./outputR")
response(me.mfeatures)
par(mfrow=c(1,2))
response(me,var="h_dem")
response(me.mfeatures,var="h_dem")

#Elith, J.et al. (2011), A statistical explanation of MaxEnt for ecologists. 
#Diversity and Distributions, 17: 43–57. 

#Argumentos para proyección
map2<-predict(me.mfeatures, layers, progress="text")
map.future2<-predict(me.mfeatures, hotlayers, 
                     args=c("extrapolate=FALSE", "doclamp=TRUE"), 
                     progress="text")
plot(map2)
plot(map.future2)

#Elith, J., Kearney, M. and Phillips, S. (2010), The art of modelling 
#range-shifting species. Methods in Ecology and Evolution, 1: 330–342. 

#Evaluar nuevo modelo
ev.stats <- evModel(fold, pres.covs, bkg.covs, factor.ind=3, mxnt.args,
                     path="./outputR")
mean(auc) #Estadisticas modelo con autofeature
mean(ev.stats$auc) #Estadisticas modelo con autofeature
#Comparación presentes
plot(map.tss)
plot(map2>mean(ev.stats$max.tss))

#Comparación futuros
plot(map.future.tss)
plot(map.future2>mean(ev.stats$max.tss))

#Dos casos especiales

###Selección de variables (Warren & Seifert 2011)
#Warren, D. L. and Seifert, S. N. (2011), Ecological niche modeling in 
#Maxent: the importance of model complexity and the performance of model
#selection criteria. Ecological Applications, 21: 335–342.

#An Introduction to Statistical Learning: with Applications in R (Capitulo 6)

betas=c(0.02, 0.1, 0.46, 1, 2.2, 4.6)
opt.lambda <- data.frame(beta=betas, nparams=NA, mean.auc=NA, mean.mtss=NA)

for(i in 1:length(betas)){
  mxnt.args=c("autofeature=FALSE",
              "linear=TRUE",
              "quadratic=TRUE",
              "product=FALSE",
              "hinge=FALSE",
              "threshold=FALSE",
              paste0("betamultiplier=",betas[i])) #Seleccionar features manualmente
  ev.stats <- evModel2(fold, pres.covs, bkg.covs, factor.ind=3,mxnt.args,
                       path="./outputR")
  opt.lambda[i,2:4]<-sapply(ev.stats,mean)
}

par(mfrow=c(1,1))
plot(opt.lambda$beta,opt.lambda$nparams,type="l",
     xlab="Beta multiplier",ylab="Numero de parámetros promedio")

plot(opt.lambda$beta,opt.lambda$mean.auc,type="l",
     xlab="Beta multiplier",ylab="AUC",col="blue")

opt.lambda

#Ver ayuda de Maxent para más flags

##Tomar background de M
#Anderson, R. P. and Raza, A. (2010), The effect of the extent of the study
#region on GIS models of species geographic distributions and estimates of
#niche evolution: preliminary tests with montane rodents (genus Nephelomys)
#in Venezuela. Journal of Biogeography, 37: 1378–1393. 

plot(layers[["ecoreg"]],col=rainbow(15))
points(occs[,2:3], pch=18,cex=0.6)
eco.mask<-layers[["ecoreg"]] %in% unique(pres.covs[,"ecoreg"])
eco.mask[eco.mask==0] <- NA
plot(eco.mask)
writeRaster(eco.mask, "./layers/ecomask.asc", overwrite=T)
layers.eco <- stack(list.files("./layers","*.asc$",full.names=T))

me.eco <- maxent(layers.eco, occs[,2:3], factors='ecoreg', removeDuplicates=TRUE, path="./outputR")
map.eco <- predict(layers.eco, me.eco, progress="text")
plot(map.eco)

eco.mask[!is.na(layers[[1]])] <- 1
plot(eco.mask)
writeRaster(!is.na(layers[["ecoreg"]]), "./layers/ecomask.asc", overwrite=T)
layers.eco <- stack(list.files("./layers","*.asc$",full.names=T))

map.proj <- predict(layers.eco, me.eco, progress="text")

par(mfrow=c(1,2))
plot(map.eco)
plot(map.proj)
