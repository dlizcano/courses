
library(raster)
options(java.parameters = "-Xmx1g" )
library(dismo)

occs <- read.csv("C:/Users/GIC 14/Dropbox/Presentaciones/Hangout_Maxent/samples/bradypus.csv")
layers<-stack(list.files("C:/Users/GIC 14/Dropbox/Presentaciones/Hangout_Maxent/layers","*.asc$",full.names=T))
#plot(layers)
plot(layers[[1]])
points(occs[,2:3])

#Run model
me <- maxent(layers, occs[,2:3], factors='ecoreg')
# Los resultados obtenidos con R y la consola son parecidos pero no iguales
# Vemos que los puntos de presencia son iguales pero el background no
# Tampoco tenemos el jpg del modelo
# Volvemos a la consola para escribir las predicciones del background

bkg <- read.csv("C:/Users/GIC 14/Dropbox/Presentaciones/Hangout_Maxent/outputs/bradypus_variegatus_backgroundPredictions.csv")
occs <- read.csv("C:/Users/GIC 14/Dropbox/Presentaciones/Hangout_Maxent/outputs/bradypus_variegatus_samplePredictions.csv")

pres.values<-extract(layers,occs[,1:2])
bkg.values<-extract(layers,bkg[,1:2])
env.values<-data.frame(rbind(pres.values,bkg.values))
env.values$ecoreg<-as.factor(env.values$ecoreg) #Importante cuando se usan variables categoricas
y<-c(rep(1,nrow(pres.values)),rep(0,nrow(bkg.values)))

me <- maxent(env.values, y, args=c("autofeature=true"))
