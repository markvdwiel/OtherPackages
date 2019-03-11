library(CGHtest)
data(cghreg)
datainfo <- data.frame(chromosomes(cghreg), bpstart(cghreg), bpend(cghreg), nclone(cghreg), avedist(cghreg))
datacgh <- regions(cghreg)
gfr <- groupfreq(datacgh,group=c(7,30),groupnames=c("MSI+", "CIN+"),af=0.1)
pvs <- pvalstest(datacgh,datainfo,teststat = "Chi-square",group=c(7,30),groupnames=c("MSI+", "CIN+"),af=0.1,niter=500)
fdrs <- fdrperm(pvs,mtdirection="stepup")


library(CGHtestpar)
data(cghreg)
datainfo <- data.frame(chromosomes(cghreg), bpstart(cghreg), bpend(cghreg), nclone(cghreg), avedist(cghreg))
datacgh <- regions(cghreg)
gfr <- groupfreq(datacgh,group=c(7,30),groupnames=c("MSI+", "CIN+"),af=0.1)

#The most time-consuming function can use multiple cpus:
pmt <- proc.time()
pvs <- pvalstest(datacgh,datainfo,teststat = "Chi-square",group=c(7,30),groupnames=c("MSI+", "CIN+"),af=0.1,niter=200,ncpus=2)
proc.time()-pmt
fdrs <- fdrperm(pvs,mtdirection="stepup")
