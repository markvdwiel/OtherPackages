library(NoWaves)
#setwd("C:\\VUData\\CGHdata\\NoWaves")
#setwd("C:\\ExternData\\AMC")
setwd("C:\\VUData\\CalibrationFiles")

#CGHTumor <- read.table("C:\\VUData\\Rebecca\\Corrected\\Corrected_Dutch_tumor.txt",header=TRUE, sep="\t", na.strings = c("NA","#N/A"),comment.char="%")
#OR load("Dutch_tumor.Rdata") #should contains object called CGHTumor 
#load("C:\\VUData\\Normals44K\\Dutch_clingenet.Rdata")
CGHNormal <- read.table("C:\\VUData\\CGHdata\\Gent\\HG_samples.txt",header=TRUE, sep="\t", na.strings = c("NA","#N/A"),comment.char="%")
CGHNormal <- read.table("C:\\VUData\\CGHdata\\Gent\\HG_samples.txt",header=TRUE, sep="\t", na.strings = c("NA","#N/A"),comment.char="%")
CGHNormal <- read.table("C:\\VUData\\CalibrationFiles\\Agilent180K.txt",header=TRUE, sep="\t", na.strings = c("NA","#N/A"),comment.char="%")
CGHNormal <- read.table("C:\\VUData\\CalibrationFiles\\NimbleGen135KNormals.txt",header=TRUE, sep="\t", na.strings = c("NA","#N/A"),comment.char="%")
CGHNormal <- data.frame(CGHNormal[,1:4],apply(CGHNormal[,-(1:4)],c(1,2),round,digits=3))
save(CGHNormal,file="C:\\VUData\\CalibrationFiles\\NimbleGen135K.Rdata")


NormalsSmooth <- SmoothNormals(CGHNormal,bandwidth=1)
#save(NormalsSmooth,file="C:\\VUData\\Normals44K\\NormalsSmooth44.Rdata")
#save(NormalsSmooth,file="C:\\VUData\\CGHdata\\Gent\\NormalsSmooth44_Gent.Rdata")
#save(NormalsSmooth,file="C:\\VUData\\CalibrationFiles\\AgilentSmooth180K.Rdata")
save(NormalsSmooth,file="C:\\VUData\\CalibrationFiles\\NimbleGenSmooth135K.Rdata")

#load("C:\\VUData\\CalibrationFiles\\AgilentSmooth44K.Rdata")
#load("C:\\VUData\\CGHdata\\Gent\\NormalsSmooth44_Gent.Rdata")
cormat <- CorTumNorm(CGHTumor,NormalsSmooth)
cormat
corrected <- CorrectTumors(CGHTumor,NormalsSmooth)
plotboth(samp=1,CGHTumor,corrected)
write.table(corrected, file="corrected_tumor.txt", sep="\t", quote=F,row.names=F)
#write.table(corrected, file="C:\\VUData\\Rebecca\\Corrected\\Corrected_Dutch_tumor.txt", sep="\t", quote=F,row.names=F)

#save(corrected, file="corrected_tumor.Rdata")
        
