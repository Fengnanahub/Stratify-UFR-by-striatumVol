
########################################UFR63 compute
# install.packages("bnlearn")
library(bnlearn)

setwd('')

AlldataforBayesianGHR63 <- read.table('', head = TRUE, sep=",")

# exclude NA and abnormal row;
AlldataforBayesianGHR63cl <- AlldataforBayesianGHR63[c(1:9, 11:56, 58:63),]

dataOrig<- AlldataforBayesianGHR63cl[, c(2:5)] # no Cov
dataOrig <- as.data.frame(scale(dataOrig, center = T, scale = T))

# regress out covariates for each var;
#(1) res1: regress age
roiRputamenvol.lm = lm(roiRputamenvol ~ age, data=AlldataforBayesianGHR63cl) # roi putamen
roiRputamenvol.res = resid(roiRputamenvol.lm)
cluCaudatevol.lm = lm(cluCaudatevol ~ age, data=AlldataforBayesianGHR63cl) # cluster caudate
cluCaudatevol.res = resid(cluCaudatevol.lm)
LESN.lm = lm(LESN ~ age, data=AlldataforBayesianGHR63cl) # LESN
LESN.res = resid(LESN.lm)
SIPStot.lm = lm(SIPStot ~ age, data=AlldataforBayesianGHR63cl) # SIPStot
SIPStot.res = resid(SIPStot.lm)

data.regressed <- cbind(roiRputamenvol.res, cluCaudatevol.res, LESN.res, SIPStot.res)
colnames(data.regressed) <- c("roiRputamenvol.res1", "cluCaudatevol.res1", "LESN.res1", "SIPStot.res1")
data.regressed <- as.data.frame(scale(data.regressed, center = T, scale = T))

#(2) res2: regress sex after regress age(through res1)
roiRputamenvol.lm = lm(roiRputamenvol.res ~ sex, data=AlldataforBayesianGHR63cl) # roi putamen
roiRputamenvol.res = resid(roiRputamenvol.lm)
cluCaudatevol.lm = lm(cluCaudatevol.res ~ sex, data=AlldataforBayesianGHR63cl) # cluster caudate
cluCaudatevol.res = resid(cluCaudatevol.lm)
LESN.lm = lm(LESN.res ~ sex, data=AlldataforBayesianGHR63cl) # LESN
LESN.res = resid(LESN.lm)
SIPStot.lm = lm(SIPStot.res ~ sex, data=AlldataforBayesianGHR63cl) # SIPStot
SIPStot.res = resid(SIPStot.lm)

data.regressed <- cbind(roiRputamenvol.res, cluCaudatevol.res, LESN.res, SIPStot.res)
colnames(data.regressed) <- c("roiRputamenvol.res2", "cluCaudatevol.res2", "LESN.res2", "SIPStot.res2")
data.regressed <- as.data.frame(scale(data.regressed, center = T, scale = T))

#(3) res3: regress edu after regress age，sex(through res1->res2)
roiRputamenvol.lm = lm(roiRputamenvol.res ~ edu, data=AlldataforBayesianGHR63cl) # roi putamen
roiRputamenvol.res = resid(roiRputamenvol.lm)
cluCaudatevol.lm = lm(cluCaudatevol.res ~ edu, data=AlldataforBayesianGHR63cl) # cluster caudate
cluCaudatevol.res = resid(cluCaudatevol.lm)
LESN.lm = lm(LESN.res ~ edu, data=AlldataforBayesianGHR63cl) # LESN
LESN.res = resid(LESN.lm)
SIPStot.lm = lm(SIPStot.res ~ edu, data=AlldataforBayesianGHR63cl) # SIPStot
SIPStot.res = resid(SIPStot.lm)

rm(data.regressed)
data.regressed <- cbind(roiRputamenvol.res, cluCaudatevol.res, LESN.res, SIPStot.res)
colnames(data.regressed) <- c("roiRputamenvol.res3", "cluCaudatevol.res3", "LESN.res3", "SIPStot.res3")
data.regressed <- as.data.frame(scale(data.regressed, center = T, scale = T))

#(4) res4: regress TIV after regress age，sex, edu(through res1->res2->res3);
roiRputamenvol.lm = lm(roiRputamenvol.res ~ TIV, data=AlldataforBayesianGHR63cl) # roi putamen
roiRputamenvol.res = resid(roiRputamenvol.lm)
cluCaudatevol.lm = lm(cluCaudatevol.res ~ TIV, data=AlldataforBayesianGHR63cl) # cluster caudate
cluCaudatevol.res = resid(cluCaudatevol.lm)
LESN.lm = lm(LESN.res ~ TIV, data=AlldataforBayesianGHR63cl) # LESN
LESN.res = resid(LESN.lm)
SIPStot.lm = lm(SIPStot.res ~ TIV, data=AlldataforBayesianGHR63cl) # SIPStot
SIPStot.res = resid(SIPStot.lm)

rm(data.regressed)
data.regressed <- cbind(roiRputamenvol.res, cluCaudatevol.res, LESN.res, SIPStot.res)
colnames(data.regressed) <- c("roiRputamenvol.res4", "cluCaudatevol.res4", "LESN.res4", "SIPStot.res4")
data.regressed <- as.data.frame(scale(data.regressed, center = T, scale = T))
# hc
bn.hc <- hc(data.regressed, blacklist = NULL, score = "bge")
plot(bn.hc)
# bootstrap
# ver 1: 1000 and 0.70
boot = boot.strength(data = data.regressed, R = 1000, m=nrow(data.regressed), algorithm = "hc", 
                     algorithm.args = list(score ="bge", iss = 10))

boot[(boot$strength > 0.7) & (boot$direction >= 0.5), ]
avg.boot = averaged.network(boot, threshold = 0.70)
plot(avg.boot)
# save(data.regressed, file = '~/finalres4_data.regressed.rda')



