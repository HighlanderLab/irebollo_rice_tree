
# Models 1-7 with PRM, SRM, and BRM
# I. Rebollo, D. Tolhurst

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
data <- read.table("Phenotype.txt", header = TRUE)
head(data)

# A = PRM with individuals in the same order as levels in data$LINE_ID
# G = SRM with individuals in the same order as levels in data$LINE_ID
# TG = BRM with individuals in the same order as levels in data$LINE_ID

##############################
# Model0: main effects only
# No YEAR_LOCATION means or blocks
# single overall mean
# single genetic and reisdual vars
##############################

# pedigree NRM
require(asreml)
asr0.A <- asreml(YIELD ~ 1,
                 random = ~ vm(LINE_ID, A),
                 residual = ~units,
                 data = data, 
                 na.action = na.method(x="include"),
                 workspace = 1e9)
asr0.A <- update(asr0.A)
# -186,223.7 
# -186,109.3 w/ non-additive geno effects - some difference

# Van Raden GRM
diag(G) <- diag(G) + 0.0000001
dim(G)
G <- G[levels(data$LINE_ID), levels(data$LINE_ID)]
asr0.G <- asreml(YIELD ~ 1,
                 random = ~ vm(LINE_ID, G),
                 residual = ~units,
                 data = data, 
                 na.action = na.method(x="include"),
                 workspace = 1e9)
asr0.G <- update(asr0.G)
# -186,017.4
# -186,010.2 w/ non-additive geno effects - no diff
# compare with pedigree NRM
cor(asr0.A$coefficients$random[grep("vm", rownames(asr0.A$coefficients$random))],
    asr0.G$coefficients$random[grep("vm", rownames(asr0.G$coefficients$random))])
# 0.95  main effects

# tree-based GRM
dim(TG)
# 936
asr0.T <- asreml(YIELD ~ 1,
                 random = ~ vm(LINE_ID, TG),
                 residual = ~units,
                 data = data,
                 na.action = na.method(x="include"),
                 workspace = 1e9)
asr0.T <- update(asr0.T)
# -186,009.5
# -186,009.4  w/ non-additive geno effects - no diff
# compare with Van Raden GRM
cor(asr0.G$coefficients$random[grep("vm", rownames(asr0.G$coefficients$random))],
    asr0.T$coefficients$random[grep("vm", rownames(asr0.T$coefficients$random))])
# 0.995  main effects

# compare with pedigree NRM
cor(asr0.A$coefficients$random[grep("vm", rownames(asr0.A$coefficients$random))],
    asr0.T$coefficients$random[grep("vm", rownames(asr0.T$coefficients$random))])
# 0.957

# look at some principal components
M <- svd(G)$u %*% diag(sqrt(svd(G)$d))
range(G - M %*% t(M))

TM <- svd(TG)$u %*% diag(sqrt(svd(TG)$d))
range(TG - TM %*% t(TM))

plot(svd(G)$u[,1], svd(TG)$u[,1])
plot(svd(TG)$d, svd(G)$d)
tt <- c()
for(i in 1:936){tt[i] <- abs(cor(M[,i], TM[,i]))}
svd(TG)$d/sum(svd(TG)$d)

plot(diag(M[,1:500] %*% t(M[,1:500])), diag(TM[,1:500] %*% t(TM[,1:500])))


##############################
# Model1: diagonal effects only
# add YEAR_LOCATION means and blocks
# single genetic var
# het block and residual vars
##############################

str(data)
with(data, table(YEAR_LOCATION))
nrow(data)
# 23,594
nlevels(data$YEAR_LOCATION)
# 49
nlevels(data$LINE_ID)
# 936
genos.YEAR_LOCATION <- with(data, tapply(Gkeep, YEAR_LOCATION, function(x) length(unique(x))))
with(data, tapply(TRIAL_ID, YEAR_LOCATION, function(x) length(unique(x))))
data$TRIAL_ID
data5 <- droplevels(data[data$YEAR_LOCATION %in% names(genos.YEAR_LOCATION)[genos.YEAR_LOCATION >= 20],])
# subset out for YEAR_LOCATIONs with more than 20 genos
nrow(data5)
# 23,009
nlevels(data5$YEAR_LOCATION)
# 40
nlevels(data5$LINE_ID)
# 936
# ok so we lost 9 YEAR_LOCATIONs, but no genos

# pedigree NRM
asr1.A <- asreml(YIELD ~ YEAR_LOCATION,
                 random = ~ diag(YEAR_LOCATION):vm(LINE_ID, A) + 
                            diag(YEAR_LOCATION):TRIAL_ID:BLOCK,
                 residual = ~dsum(~units|YEAR_LOCATION),
                 data = data5, 
                 na.action = na.method(x="include"),
                 workspace = 1e9)
asr1.A <- update(asr1.A)
# -171,666.2


# Van Raden GRM
asr1.G <- asreml(YIELD ~ YEAR_LOCATION,
                 random = ~ diag(YEAR_LOCATION):vm(LINE_ID, G) + 
                            diag(YEAR_LOCATION):TRIAL_ID:BLOCK,
                 residual = ~dsum(~units|YEAR_LOCATION),
                 data = data5, 
                 na.action = na.method(x="include"),
                 workspace = 1e9)
asr1.G <- update(asr1.G)
# -171,527.4
# compare with pedigree NRM
cor(asr1.A$coefficients$random[grep("vm", rownames(asr1.A$coefficients$random))],
    asr1.G$coefficients$random[grep("vm", rownames(asr1.G$coefficients$random))])
# 0.67  int. effects
require(ASExtras4)
met.asreml(asr1.G) 
plot(asr1.G)
# look at residuals for potential outliers


# tree-based GRM
asr1.T <- asreml(YIELD ~ YEAR_LOCATION,
                 random = ~ diag(YEAR_LOCATION):vm(LINE_ID, TG) + 
                            diag(YEAR_LOCATION):TRIAL_ID:BLOCK,
                 residual = ~dsum(~units|YEAR_LOCATION),
                 data = data5,
                 na.action = na.method(x="include"),
                 workspace = 1e9)
asr1.T <- update(asr1.T)
# -171,487.0
# compare with Van Raden GRM
cor(asr1.G$coefficients$random[grep("vm", rownames(asr1.G$coefficients$random))],
    asr1.T$coefficients$random[grep("vm", rownames(asr1.T$coefficients$random))])
# 0.99  int. effects

cor(asr1.A$coefficients$random[grep("vm", rownames(asr1.A$coefficients$random))],
    asr1.T$coefficients$random[grep("vm", rownames(asr1.T$coefficients$random))])
# 0.68  int. effects


##############################
# Model2: main effects only
# add YEAR_LOCATION means and blocks
# single genetic var
# het block and residual vars
##############################

# pedigree NRM
asr2.A <- asreml(YIELD ~ YEAR_LOCATION,
                   random = ~ vm(LINE_ID, A) +
                              diag(YEAR_LOCATION):TRIAL_ID:BLOCK,
                   residual = ~dsum(~units|YEAR_LOCATION),
                   data = data5, 
                   na.action = na.method(x="include"),
                   workspace = 1e9)
asr2.A <- update(asr2.A)
# -172,530.8  
# -172,498.7  w/ non-additive geno effects - small diff


# Van Raden GRM
asr2.G <- asreml(YIELD ~ YEAR_LOCATION,
                   random = ~ vm(LINE_ID, G) +
                              diag(YEAR_LOCATION):TRIAL_ID:BLOCK,
                   residual = ~dsum(~units|YEAR_LOCATION),
                   data = data5, 
                   na.action = na.method(x="include"),
                   workspace = 1e9)
asr2.G <- update(asr2.G)
# -172,408.4
# -172,398.3  w/ non-additive geno effects - no diff
# compare with pedigree NRM
cor(asr2.A$coefficients$random[grep("vm", rownames(asr2.A$coefficients$random))],
    asr2.G$coefficients$random[grep("vm", rownames(asr2.G$coefficients$random))])
# 0.940  main effects


# tree-based GRM
asr2.T <- asreml(YIELD ~ YEAR_LOCATION,
                   random = ~ vm(LINE_ID, TG) + 
                              diag(YEAR_LOCATION):TRIAL_ID:BLOCK,
                   residual = ~dsum(~units|YEAR_LOCATION),
                   data = data5,
                   na.action = na.method(x="include"),
                   workspace = 1e9)
asr2.T <- update(asr2.T)
# -172,399.4
# -172,399.4  w/ non-additive geno effects - no diff
# leave the non-add effects out here after
# compare with Van Raden GRM
cor(asr2.G$coefficients$random[grep("vm", rownames(asr2.G$coefficients$random))],
    asr2.T$coefficients$random[grep("vm", rownames(asr2.T$coefficients$random))])
# 0.993  main effects
cor(asr2.A$coefficients$random[grep("vm", rownames(asr2.A$coefficients$random))],
    asr2.T$coefficients$random[grep("vm", rownames(asr2.T$coefficients$random))])
# 0.957

##############################
# Model3: Compound symmetry
# add YEAR_LOCATION means and blocks
# single main and int. genetic vars
# het block and residual vars
##############################

# pedigree NRM
asr3.A <- asreml(YIELD ~ YEAR_LOCATION,
                 random = ~ vm(LINE_ID, A) + idv(YEAR_LOCATION):vm(LINE_ID, A) + 
                            diag(YEAR_LOCATION):TRIAL_ID:BLOCK,
                 residual = ~dsum(~units|YEAR_LOCATION),
                 data = data5, 
                 na.action = na.method(x="include"),
                 workspace = 1e9)
asr3.A <- update(asr3.A)
# -171,477.2 


# Van Raden GRM
asr3.G <- asreml(YIELD ~ YEAR_LOCATION,
                 random = ~ vm(LINE_ID, G) + idv(YEAR_LOCATION):vm(LINE_ID, G) + 
                            diag(YEAR_LOCATION):TRIAL_ID:BLOCK,
                 residual = ~dsum(~units|YEAR_LOCATION),
                 data = data5,
                 na.action = na.method(x="include"),
                 workspace = 1e9)
asr3.G <- update(asr3.G)
# -171,324.4
# compare with pedigree NRM
cor(asr3.A$coefficients$random[grep("vm", rownames(asr3.A$coefficients$random))],
    asr3.G$coefficients$random[grep("vm", rownames(asr3.G$coefficients$random))])
# 0.63  main effects
cor(asr3.A$coefficients$random[grep("YEAR_LOCATION.*vm", rownames(asr3.A$coefficients$random))],
    asr3.G$coefficients$random[grep("YEAR_LOCATION.*vm", rownames(asr3.G$coefficients$random))])
# 0.62  int. effects
cor(asr3.A$coefficients$random[grep("^vm", rownames(asr3.A$coefficients$random))] + asr3.A$coefficients$random[grep("YEAR_LOCATION.*vm", rownames(asr3.A$coefficients$random))],
    asr3.G$coefficients$random[grep("^vm", rownames(asr3.G$coefficients$random))] + asr3.G$coefficients$random[grep("YEAR_LOCATION.*vm", rownames(asr3.G$coefficients$random))])
# 0.81  main + int. effects


# tree-based GRM
asr3.T <- asreml(YIELD ~ YEAR_LOCATION,
                 random = ~ vm(LINE_ID, TG) + idv(YEAR_LOCATION):vm(LINE_ID, TG) + 
                            diag(YEAR_LOCATION):TRIAL_ID:BLOCK,
                 residual = ~dsum(~units|YEAR_LOCATION),
                 data = data5,
                 na.action = na.method(x="include"),
                 workspace = 1e9)
asr3.T <- update(asr3.T)
# -171,288.9
# compare with Van Raden GRM
cor(asr3.G$coefficients$random[grep("^vm", rownames(asr3.G$coefficients$random))],
    asr3.T$coefficients$random[grep("^vm", rownames(asr3.T$coefficients$random))])
# 0.99  main effects
cor(asr3.G$coefficients$random[grep("YEAR_LOCATION.*vm", rownames(asr3.G$coefficients$random))],
    asr3.T$coefficients$random[grep("YEAR_LOCATION.*vm", rownames(asr3.T$coefficients$random))])
# 0.99  int. effects
cor(asr3.G$coefficients$random[grep("^vm", rownames(asr3.G$coefficients$random))] + asr3.G$coefficients$random[grep("YEAR_LOCATION.*vm", rownames(asr3.G$coefficients$random))],
    asr3.T$coefficients$random[grep("^vm", rownames(asr3.T$coefficients$random))] + asr3.T$coefficients$random[grep("YEAR_LOCATION.*vm", rownames(asr3.T$coefficients$random))])
# 0.99  main + int. effects

cor(asr3.A$coefficients$random[grep("^vm", rownames(asr3.A$coefficients$random))],
    asr3.T$coefficients$random[grep("^vm", rownames(asr3.T$coefficients$random))])
# 0.91  main effects
cor(asr3.A$coefficients$random[grep("YEAR_LOCATION.*vm", rownames(asr3.A$coefficients$random))],
    asr3.T$coefficients$random[grep("YEAR_LOCATION.*vm", rownames(asr3.T$coefficients$random))])
# 0.63  int. effects
cor(asr3.A$coefficients$random[grep("^vm", rownames(asr3.A$coefficients$random))] + asr3.A$coefficients$random[grep("YEAR_LOCATION.*vm", rownames(asr3.A$coefficients$random))],
    asr3.T$coefficients$random[grep("^vm", rownames(asr3.T$coefficients$random))] + asr3.T$coefficients$random[grep("YEAR_LOCATION.*vm", rownames(asr3.T$coefficients$random))])
# 0.83  main + int. effects

##############################
# Model4: Main effect + diag
# add YEAR_LOCATION means and blocks
# single main and het int. genetic vars
# het block and residual vars
##############################

# pedigree NRM
asr4.A <- asreml(YIELD ~ YEAR_LOCATION,
                 random = ~ vm(LINE_ID, A) + diag(YEAR_LOCATION):vm(LINE_ID, A) + 
                            diag(YEAR_LOCATION):TRIAL_ID:BLOCK,
                 residual = ~dsum(~units|YEAR_LOCATION),
                 data = data5, 
                 na.action = na.method(x="include"),
                 workspace = 1e9)
asr4.A <- update(asr4.A)
# -171,380.8
max(summary(asr4.A)$varcom[,5], na.rm=T)


# Van Raden GRM
asr4.G <- asreml(YIELD ~ YEAR_LOCATION,
                    random = ~ vm(LINE_ID, G) + diag(YEAR_LOCATION):vm(LINE_ID, G) + 
                               diag(YEAR_LOCATION):TRIAL_ID:BLOCK,
                    residual = ~dsum(~units|YEAR_LOCATION),
                    data = data5,
                    na.action = na.method(x="include"),
                    workspace = 1e9)
asr4.G <- update(asr4.G)
# -171,240.6
max(summary(asr4.G)$varcom[,5], na.rm=T)
# compare with pedigree NRM
cor(asr4.A$coefficients$random[grep("vm", rownames(asr4.A$coefficients$random))],
    asr4.G$coefficients$random[grep("vm", rownames(asr4.G$coefficients$random))])
# 0.60  main effects
cor(asr4.A$coefficients$random[grep("YEAR_LOCATION.*vm", rownames(asr4.A$coefficients$random))],
    asr4.G$coefficients$random[grep("YEAR_LOCATION.*vm", rownames(asr4.G$coefficients$random))])
# 0.59  int. effects
cor(asr4.A$coefficients$random[grep("^vm", rownames(asr4.A$coefficients$random))] + asr4.A$coefficients$random[grep("YEAR_LOCATION.*vm", rownames(asr4.A$coefficients$random))],
    asr4.G$coefficients$random[grep("^vm", rownames(asr4.G$coefficients$random))] + asr4.G$coefficients$random[grep("YEAR_LOCATION.*vm", rownames(asr4.G$coefficients$random))])
# 0.81  main + int. effects


# tree GRM
asr4.T <- asreml(YIELD ~ YEAR_LOCATION,
                     random = ~ vm(LINE_ID, TG) + diag(YEAR_LOCATION):vm(LINE_ID, TG) + 
                                diag(YEAR_LOCATION):TRIAL_ID:BLOCK,
                     residual = ~dsum(~units|YEAR_LOCATION),
                     data = data5,
                     na.action = na.method(x="include"),
                     workspace = 1e9)
asr4.T <- update(asr4.T)
# -171,206.3
max(summary(asr4.T)$varcom[,5], na.rm=T)
# compare with Van Raden GRM
cor(asr4.G$coefficients$random[grep("^vm", rownames(asr4.G$coefficients$random))],
    asr4.T$coefficients$random[grep("^vm", rownames(asr4.T$coefficients$random))])
# 0.995  main effects
cor(asr4.G$coefficients$random[grep("YEAR_LOCATION.*vm", rownames(asr4.G$coefficients$random))],
    asr4.T$coefficients$random[grep("YEAR_LOCATION.*vm", rownames(asr4.T$coefficients$random))])
# 0.98  int. effects
cor(asr4.G$coefficients$random[grep("^vm", rownames(asr4.G$coefficients$random))] + asr4.G$coefficients$random[grep("YEAR_LOCATION.*vm", rownames(asr4.G$coefficients$random))],
    asr4.T$coefficients$random[grep("^vm", rownames(asr4.T$coefficients$random))] + asr4.T$coefficients$random[grep("YEAR_LOCATION.*vm", rownames(asr4.T$coefficients$random))])
# 0.992  main + int. effects

# compare with pedigree NRM
cor(asr4.A$coefficients$random[grep("^vm", rownames(asr4.A$coefficients$random))],
    asr4.T$coefficients$random[grep("^vm", rownames(asr4.T$coefficients$random))])
# 0.911  main effects
cor(asr4.A$coefficients$random[grep("YEAR_LOCATION.*vm", rownames(asr4.A$coefficients$random))],
    asr4.T$coefficients$random[grep("YEAR_LOCATION.*vm", rownames(asr4.T$coefficients$random))])
# 0.609  int. effects
cor(asr4.A$coefficients$random[grep("^vm", rownames(asr4.A$coefficients$random))] + asr4.A$coefficients$random[grep("YEAR_LOCATION.*vm", rownames(asr4.A$coefficients$random))],
    asr4.T$coefficients$random[grep("^vm", rownames(asr4.T$coefficients$random))] + asr4.T$coefficients$random[grep("YEAR_LOCATION.*vm", rownames(asr4.T$coefficients$random))])
# 0.823  main + int. effects


##############################
# Model5: rr1 + diag (FA1)
# add YEAR_LOCATION means and blocks
# single main and het int. genetic vars
# het block and residual vars
##############################

# pedigree NRM
asr5.A <- asreml(YIELD ~ YEAR_LOCATION,
                 random = ~ rr(YEAR_LOCATION,1):vm(LINE_ID, A) +  diag(YEAR_LOCATION):vm(LINE_ID, A) + 
                            diag(YEAR_LOCATION):TRIAL_ID:BLOCK,
                 residual = ~dsum(~units|YEAR_LOCATION),
                 data = data5,
                 na.action = na.method(x="include"),
                 workspace = 1e9)
asr5.A <- update(asr5.A)
# -171,334.7
max(summary(asr5.A)$varcom[,5], na.rm=T)


# Van Raden GRM
asr5.G <- asreml(YIELD ~ YEAR_LOCATION,
                 random = ~ rr(YEAR_LOCATION,1):vm(LINE_ID, G) +  diag(YEAR_LOCATION):vm(LINE_ID, G) + 
                            diag(YEAR_LOCATION):TRIAL_ID:BLOCK,
                 residual = ~dsum(~units|YEAR_LOCATION),
                 data = data5,
                 na.action = na.method(x="include"),
                 workspace = 1e9)
asr5.G <- update(asr5.G)
# -171,202.6
max(summary(asr5.G)$varcom[,5], na.rm=T)
# compare with pedigree NRM
cor(asr5.A$coefficients$random[grep("Comp", rownames(asr5.A$coefficients$random))],
    asr5.G$coefficients$random[grep("Comp", rownames(asr5.G$coefficients$random))])
# 0.91  scores (f)
cor(asr5.A$coefficients$random[grep("^rr.*-U", rownames(asr5.A$coefficients$random))],
    asr5.G$coefficients$random[grep("^rr.*-U", rownames(asr5.G$coefficients$random))])
# 0.89  regression BLUPs (Lamf)
cor(asr5.A$coefficients$random[grep("^YEAR_LOCATION.*vm", rownames(asr5.A$coefficients$random))],
    asr5.G$coefficients$random[grep("^YEAR_LOCATION.*vm", rownames(asr5.G$coefficients$random))])
# 0.57  specific BLUPs  (delta)
cor(asr5.A$coefficients$random[grep("^rr.*-U", rownames(asr5.A$coefficients$random))] + asr5.A$coefficients$random[grep("^YEAR_LOCATION.*vm", rownames(asr5.A$coefficients$random))],
    asr5.G$coefficients$random[grep("^rr.*-U", rownames(asr5.G$coefficients$random))] + asr5.G$coefficients$random[grep("^YEAR_LOCATION.*vm", rownames(asr5.G$coefficients$random))])
# 0.84  total BLUPs (Lamf + delta)

# tree-based GRM
asr5.T <- asreml(YIELD ~ YEAR_LOCATION,
                 random = ~ rr(YEAR_LOCATION,1):vm(LINE_ID, TG) + diag(YEAR_LOCATION):vm(LINE_ID, TG) + 
                            diag(YEAR_LOCATION):TRIAL_ID:BLOCK,
                 residual = ~dsum(~units|YEAR_LOCATION),
                 data = data5,
                 na.action = na.method(x="include"),
                 workspace = 1e9)
asr5.T <- update(asr5.T)
# -171,169.5
max(summary(asr5.T)$varcom[,5], na.rm=T)
# compare with Van Raden GRM
cor(asr5.G$coefficients$random[grep("Comp", rownames(asr5.G$coefficients$random))],
    asr5.T$coefficients$random[grep("Comp", rownames(asr5.T$coefficients$random))])
# 0.994  scores (f), negativE????
cor(asr5.G$coefficients$random[grep("^rr.*-U", rownames(asr5.G$coefficients$random))],
    asr5.T$coefficients$random[grep("^rr.*-U", rownames(asr5.T$coefficients$random))])
# 0.991  regression BLUPs (Lamf)
cor(asr5.G$coefficients$random[grep("^YEAR_LOCATION.*vm", rownames(asr5.G$coefficients$random))],
    asr5.T$coefficients$random[grep("^YEAR_LOCATION.*vm", rownames(asr5.T$coefficients$random))])
# 0.978  specific BLUPs  (delta)
cor(asr5.G$coefficients$random[grep("^rr.*-U", rownames(asr5.G$coefficients$random))] + asr5.G$coefficients$random[grep("^YEAR_LOCATION.*vm", rownames(asr5.G$coefficients$random))],
    asr5.T$coefficients$random[grep("^rr.*-U", rownames(asr5.T$coefficients$random))] + asr5.T$coefficients$random[grep("^YEAR_LOCATION.*vm", rownames(asr5.T$coefficients$random))])
# 0.991  total BLUPs (Lamf + delta)

# compare withpedigree NRM
cor(asr5.A$coefficients$random[grep("Comp", rownames(asr5.A$coefficients$random))],
    asr5.T$coefficients$random[grep("Comp", rownames(asr5.T$coefficients$random))])
# -0.923  scores (f)
cor(asr5.A$coefficients$random[grep("^rr.*-U", rownames(asr5.A$coefficients$random))],
    asr5.T$coefficients$random[grep("^rr.*-U", rownames(asr5.T$coefficients$random))])
# 0.907  regression BLUPs (Lamf)
cor(asr5.A$coefficients$random[grep("^YEAR_LOCATION.*vm", rownames(asr5.A$coefficients$random))],
    asr5.T$coefficients$random[grep("^YEAR_LOCATION.*vm", rownames(asr5.T$coefficients$random))])
# 0.583  specific BLUPs  (delta)
cor(asr5.A$coefficients$random[grep("^rr.*-U", rownames(asr5.A$coefficients$random))] + asr5.A$coefficients$random[grep("^YEAR_LOCATION.*vm", rownames(asr5.A$coefficients$random))],
    asr5.T$coefficients$random[grep("^rr.*-U", rownames(asr5.T$coefficients$random))] + asr5.T$coefficients$random[grep("^YEAR_LOCATION.*vm", rownames(asr5.T$coefficients$random))])
# 0.853  total BLUPs (Lamf + delta)

##############################
# Model6: rr2 + diag (FA2)
# add YEAR_LOCATION means and blocks
# single main and het int. genetic vars
# het block and residual vars
##############################

# pedigree NRM
asr6.A <- asreml(YIELD ~ YEAR_LOCATION,
                 random = ~ rr(YEAR_LOCATION,2):vm(LINE_ID, A) + diag(YEAR_LOCATION):vm(LINE_ID, A) + 
                            diag(YEAR_LOCATION):TRIAL_ID:BLOCK,
                 residual = ~dsum(~units|YEAR_LOCATION),
                 data = data5,
                 na.action = na.method(x="include"),
                 workspace = 1e9)
asr6.A <- update(asr6.A)
# DT: -173,160.0 # I got different, -171,245.4
asr6.A$loglik 
max(summary(asr6.A)$varcom[,5], na.rm=T)


# Van Raden GRM
asr6.G <- asreml(YIELD ~ YEAR_LOCATION,
                 random = ~ rr(YEAR_LOCATION,2):vm(LINE_ID, G) + diag(YEAR_LOCATION):vm(LINE_ID, G) + 
                            diag(YEAR_LOCATION):TRIAL_ID:BLOCK,
                 residual = ~dsum(~units|YEAR_LOCATION),
                 data = data5,
                 na.action = na.method(x="include"),
                 workspace = 1e9)
asr6.G <- update(asr6.G)
# DT: -173,160 , I got -171,145.3
asr6.G$loglik
max(summary(asr6.G)$varcom[,5], na.rm=T)
# compare with pedigree NRM
cor(asr6.A$coefficients$random[grep("Comp", rownames(asr6.A$coefficients$random))],
    asr6.G$coefficients$random[grep("Comp", rownames(asr6.G$coefficients$random))])
# 0.77  scores (f) I got 0.39
cor(asr6.A$coefficients$random[grep("^rr.*-U", rownames(asr6.A$coefficients$random))],
    asr6.G$coefficients$random[grep("^rr.*-U", rownames(asr6.G$coefficients$random))])
# 0.85  regression BLUPs (Lamf) I got 0.84
cor(asr6.A$coefficients$random[grep("^YEAR_LOCATION.*vm", rownames(asr6.A$coefficients$random))],
    asr6.G$coefficients$random[grep("^YEAR_LOCATION.*vm", rownames(asr6.G$coefficients$random))])
# 0.47  specific BLUPs  (delta) I get the same
cor(asr6.A$coefficients$random[grep("^rr.*-U", rownames(asr6.A$coefficients$random))] + asr6.A$coefficients$random[grep("^YEAR_LOCATION.*vm", rownames(asr6.A$coefficients$random))],
    asr6.G$coefficients$random[grep("^rr.*-U", rownames(asr6.G$coefficients$random))] + asr6.G$coefficients$random[grep("^YEAR_LOCATION.*vm", rownames(asr6.G$coefficients$random))])
# 0.85  total BLUPs (Lamf + delta) I get the same

# tree-based GRM
asr6.T <- asreml(YIELD ~ YEAR_LOCATION,
                   random = ~ rr(YEAR_LOCATION,2):vm(LINE_ID, TG) + diag(YEAR_LOCATION):vm(LINE_ID, TG) + 
                              diag(YEAR_LOCATION):TRIAL_ID:BLOCK,
                   residual = ~dsum(~units|YEAR_LOCATION),
                   data = data5,
                   na.action = na.method(x="include"),
                   workspace = 1e9)
asr6.T <- update(asr6.T)
# -173,160.0 # I got 171,113.8
asr6.T$loglik
max(summary(asr6.T)$varcom[,5], na.rm=T)
# compare with Van Raden GRM
cor(rep(c(1,-1), each = 936) * asr6.G$coefficients$random[grep("Comp", rownames(asr6.G$coefficients$random))],
    asr6.T$coefficients$random[grep("Comp", rownames(asr6.T$coefficients$random))])
# DT: 0.99  scores (f) I get 0.85
cor(asr6.G$coefficients$random[grep("^rr.*-U", rownames(asr6.G$coefficients$random))],
    asr6.T$coefficients$random[grep("^rr.*-U", rownames(asr6.T$coefficients$random))])
# DT: 0.99  regression BLUPs (Lamf) I get 0.98
cor(asr6.G$coefficients$random[grep("^YEAR_LOCATION.*vm", rownames(asr6.G$coefficients$random))],
    asr6.T$coefficients$random[grep("^YEAR_LOCATION.*vm", rownames(asr6.T$coefficients$random))])
# DT: 0.95  specific BLUPs  (delta) I get 0.88
cor(asr6.G$coefficients$random[grep("^rr.*-U", rownames(asr6.G$coefficients$random))] + asr6.G$coefficients$random[grep("^YEAR_LOCATION.*vm", rownames(asr6.G$coefficients$random))],
    asr6.T$coefficients$random[grep("^rr.*-U", rownames(asr6.T$coefficients$random))] + asr6.T$coefficients$random[grep("^YEAR_LOCATION.*vm", rownames(asr6.T$coefficients$random))])
# DT: 0.99  total BLUPs (Lamf + delta) I get same

# compare with pedigree NRM
cor(rep(c(1,-1), each = 936) * asr6.A$coefficients$random[grep("Comp", rownames(asr6.A$coefficients$random))],
    asr6.T$coefficients$random[grep("Comp", rownames(asr6.T$coefficients$random))])
# -0.04  scores (f)
cor(asr6.A$coefficients$random[grep("^rr.*-U", rownames(asr6.A$coefficients$random))],
    asr6.T$coefficients$random[grep("^rr.*-U", rownames(asr6.T$coefficients$random))])
# 0.85=  regression BLUPs (Lamf)
cor(asr6.A$coefficients$random[grep("^YEAR_LOCATION.*vm", rownames(asr6.A$coefficients$random))],
    asr6.T$coefficients$random[grep("^YEAR_LOCATION.*vm", rownames(asr6.T$coefficients$random))])
# 0.46  specific BLUPs  (delta)
cor(asr6.A$coefficients$random[grep("^rr.*-U", rownames(asr6.A$coefficients$random))] + asr6.A$coefficients$random[grep("^YEAR_LOCATION.*vm", rownames(asr6.A$coefficients$random))],
    asr6.T$coefficients$random[grep("^rr.*-U", rownames(asr6.T$coefficients$random))] + asr6.T$coefficients$random[grep("^YEAR_LOCATION.*vm", rownames(asr6.T$coefficients$random))])
# 0.86  total BLUPs (Lamf + delta)

##############################
# Model7: rr2 + diag (FA3)
# add YEAR_LOCATION means and blocks
# single main and het int. genetic vars
# het block and residual vars
##############################

# pedigree NRM
asr7.A <- asreml(YIELD ~ YEAR_LOCATION,
                 random = ~ rr(YEAR_LOCATION,3):vm(LINE_ID, A) + diag(YEAR_LOCATION):vm(LINE_ID, A) + 
                     diag(YEAR_LOCATION):TRIAL_ID:BLOCK,
                 residual = ~dsum(~units|YEAR_LOCATION),
                 data = data5,
                 na.action = na.method(x="include"),
                 workspace = 1e9)
asr7.A <- update(asr7.A)
# -171204.9
asr7.A$loglik 
max(summary(asr7.A)$varcom[,5], na.rm=T)


# Van Raden GRM
asr7.G <- asreml(YIELD ~ YEAR_LOCATION,
                 random = ~ rr(YEAR_LOCATION,3):vm(LINE_ID, G) + diag(YEAR_LOCATION):vm(LINE_ID, G) + 
                     diag(YEAR_LOCATION):TRIAL_ID:BLOCK,
                 residual = ~dsum(~units|YEAR_LOCATION),
                 data = data5,
                 na.action = na.method(x="include"),
                 workspace = 1e9)
asr7.G <- update(asr7.G)
# -171145.3
asr7.G$loglik
max(summary(asr7.G)$varcom[,5], na.rm=T)
# compare with pedigree NRM
cor(asr7.A$coefficients$random[grep("Comp", rownames(asr7.A$coefficients$random))],
    asr7.G$coefficients$random[grep("Comp", rownames(asr7.G$coefficients$random))])
# 0.05  scores (f) 
cor(asr7.A$coefficients$random[grep("^rr.*-U", rownames(asr7.A$coefficients$random))],
    asr7.G$coefficients$random[grep("^rr.*-U", rownames(asr7.G$coefficients$random))])
# 0.03  regression BLUPs (Lamf) 
cor(asr7.A$coefficients$random[grep("^YEAR_LOCATION.*vm", rownames(asr7.A$coefficients$random))],
    asr7.G$coefficients$random[grep("^YEAR_LOCATION.*vm", rownames(asr7.G$coefficients$random))])
# 0.04  specific BLUPs  (delta) 
cor(asr7.A$coefficients$random[grep("^rr.*-U", rownames(asr7.A$coefficients$random))] + asr7.A$coefficients$random[grep("^YEAR_LOCATION.*vm", rownames(asr7.A$coefficients$random))],
    asr7.G$coefficients$random[grep("^rr.*-U", rownames(asr7.G$coefficients$random))] + asr7.G$coefficients$random[grep("^YEAR_LOCATION.*vm", rownames(asr7.G$coefficients$random))])
# 0.04 total BLUPs (Lamf + delta) 


# tree-based GRM
asr7.T <- asreml(YIELD ~ YEAR_LOCATION,
                 random = ~ rr(YEAR_LOCATION,3):vm(LINE_ID, TG) + diag(YEAR_LOCATION):vm(LINE_ID, TG) + 
                     diag(YEAR_LOCATION):TRIAL_ID:BLOCK,
                 residual = ~dsum(~units|YEAR_LOCATION),
                 data = data5,
                 na.action = na.method(x="include"),
                 workspace = 1e9)
asr7.T <- update(asr7.T)
# -171069.6
asr7.T$loglik
max(summary(asr7.T)$varcom[,5], na.rm=T)
# compare with Van Raden GRM
cor(rep(c(-1, 1, -1), each = 936) *  # ???? why?
    asr7.G$coefficients$random[grep("Comp", rownames(asr7.G$coefficients$random))],
    asr7.T$coefficients$random[grep("Comp", rownames(asr7.T$coefficients$random))])
# 0.30  scores (f) 
cor(asr7.G$coefficients$random[grep("^rr.*-U", rownames(asr7.G$coefficients$random))],
    asr7.T$coefficients$random[grep("^rr.*-U", rownames(asr7.T$coefficients$random))])
# 0.90  regression BLUPs (Lamf) 
cor(asr7.G$coefficients$random[grep("^YEAR_LOCATION.*vm", rownames(asr7.G$coefficients$random))],
    asr7.T$coefficients$random[grep("^YEAR_LOCATION.*vm", rownames(asr7.T$coefficients$random))])
# 0.51  specific BLUPs  (delta) 
cor(asr7.G$coefficients$random[grep("^rr.*-U", rownames(asr7.G$coefficients$random))] + asr7.G$coefficients$random[grep("^YEAR_LOCATION.*vm", rownames(asr7.G$coefficients$random))],
    asr7.T$coefficients$random[grep("^rr.*-U", rownames(asr7.T$coefficients$random))] + asr7.T$coefficients$random[grep("^YEAR_LOCATION.*vm", rownames(asr7.T$coefficients$random))])
# 0.94  total BLUPs (Lamf + delta) 

# compare with pedigree NRM
cor(rep(c(-1, -1, 1), each = 936) * asr7.A$coefficients$random[grep("Comp", rownames(asr7.A$coefficients$random))],
    asr7.T$coefficients$random[grep("Comp", rownames(asr7.T$coefficients$random))])
# 0.03  scores (f)
cor(asr7.A$coefficients$random[grep("^rr.*-U", rownames(asr7.A$coefficients$random))],
    asr7.T$coefficients$random[grep("^rr.*-U", rownames(asr7.T$coefficients$random))])
# 0.04  regression BLUPs (Lamf)
cor(asr7.A$coefficients$random[grep("^YEAR_LOCATION.*vm", rownames(asr7.A$coefficients$random))],
    asr7.T$coefficients$random[grep("^YEAR_LOCATION.*vm", rownames(asr7.T$coefficients$random))])
# 0.06  specific BLUPs  (delta)
cor(asr7.A$coefficients$random[grep("^rr.*-U", rownames(asr7.A$coefficients$random))] + asr7.A$coefficients$random[grep("^YEAR_LOCATION.*vm", rownames(asr7.A$coefficients$random))],
    asr7.T$coefficients$random[grep("^rr.*-U", rownames(asr7.T$coefficients$random))] + asr7.T$coefficients$random[grep("^YEAR_LOCATION.*vm", rownames(asr7.T$coefficients$random))])
# 0.04 total BLUPs (Lamf + delta)

u_main <- c(asr2.A$coefficients$random[grep("^vm", rownames(asr2.A$coefficients$random))],
            asr2.G$coefficients$random[grep("^vm", rownames(asr2.G$coefficients$random))],
            asr2.T$coefficients$random[grep("^vm", rownames(asr2.T$coefficients$random))],
            asr4.A$coefficients$random[grep("^vm", rownames(asr4.A$coefficients$random))], 
            asr4.G$coefficients$random[grep("^vm", rownames(asr4.G$coefficients$random))],
            asr4.T$coefficients$random[grep("^vm", rownames(asr4.T$coefficients$random))])

length(u_main)
936 * 6 #perfect!
summary(u_main)

u_int <- c(asr4.A$coefficients$random[grep("^YEAR_LOCATION.*vm", rownames(asr4.A$coefficients$random))],
           asr4.G$coefficients$random[grep("^YEAR_LOCATION.*vm", rownames(asr4.G$coefficients$random))],
           asr4.T$coefficients$random[grep("^YEAR_LOCATION.*vm", rownames(asr4.T$coefficients$random))])
length(u_int)
# model 2 doesn't have interaction effects

saveRDS(u_main, "mainBV_models2&4.rds")
saveRDS(u_int, "intBV_models2&4.rds")

# Get SE from models ------------------------------------------------------

asr2.G_pred <- predict(asr2.G, classify = "LINE_ID", only = "vm(LINE_ID, G)", maxit = 2, vcov = T, pworkspace = 6e8)
asr2.T_pred <- predict(asr2.T, classify = "LINE_ID", only = "vm(LINE_ID, TG)", maxit = 2, vcov = T, pworkspace = 6e8)
asr4.G_pred <- predict(asr4.G, classify = "LINE_ID", only = "vm(LINE_ID, G)", maxit = 2, vcov = T, pworkspace = 6e8)
asr4.T_pred <- predict(asr4.T, classify = "LINE_ID", only = "vm(LINE_ID, TG)", maxit = 2, vcov = T, pworkspace = 6e8)
# You may need to increase the prediction workspace, you can do this by adding pworkspace=6e8  to the predict call. 

# In asreml(fixed = YIELD ~ YEAR_LOCATION, random = ~vm(LINE_ID, G) + diag(YEAR_LOCATION):TRIAL_ID:BLOCK,  :
#             Log-likelihood not converged
# add more space or is it something else?

#The PEV matrix will be output in the vcov slot of the list. 
PEVs <- list()
PEVs[["2_G"]] <- asr2.G_pred$vcov
PEVs[["2_T"]] <- asr2.T_pred$vcov
PEVs[["4_G"]] <- asr4.G_pred$vcov
PEVs[["4_T"]] <- asr4.T_pred$vcov

#get sigma_g^2 from asreml summary object (variance of BVs)
summary(asr2.G)$varcom[1,] #531764.1
summary(asr2.T)$varcom[1,] #39409468
summary(asr4.G)$varcom[1,] #421298.9
summary(asr4.T)$varcom[1,] #29878180

# Variance explained ------------------------------------------------------

# Model 1 and 2: NA
#Model 3: sigma_g^2/(sigma_g^2 + sigma_ge^2), where sigma_g^2 is the main effect variance component and sigma_ge^2 is the interaction variance component

model3_VAF <- function(model) {
    v_g <- summary(model)$varcom[grep("^vm\\(LINE_ID,", rownames(summary(model)$varcom)), "component"]
    v_ge <- summary(model)$varcom[grep("^YEAR_LOCATION:vm\\(LINE_ID, ", rownames(summary(model)$varcom)), "component"]
    # v_g <- summary(model)$varcom[1, "component"]
    # v_ge <- summary(model)$varcom[2, "component"]
    round(v_g / (v_g + v_ge) * 100, 2)
}
model3_VAF(asr3.A) # 43.67
model3_VAF(asr3.G) # 54.82
model3_VAF(asr3.T) # 47.57

# Model 4: sigma_g^2/(sigma_g^2 + sum(sigma_ge_j^2)/p) where sigma_ge_j^2 is the interaction variance component for environment j (the diag variance parameters) and p is the number of environments

model4_VAF <- function(model) {
    # model <- asr4.A
    v_g <- summary(model)$varcom[grep("^vm\\(LINE_ID,", rownames(summary(model)$varcom)), "component"]
    v_ges <- summary(model)$varcom[grep("^YEAR_LOCATION:vm\\(LINE_ID, ", rownames(summary(model)$varcom)), "component"]
    # v_g <- summary(model)$varcom[1, "component"]
    # v_ge <- summary(model)$varcom[2, "component"]
    round(v_g / (v_g + sum(v_ges)/length(v_ges)) * 100, 2)
}
model4_VAF(asr4.A) # 35.45
model4_VAF(asr4.G) # 50.65
model4_VAF(asr4.T) # 42.71


# Models 5 and 6: sum(diag(LamLam’))/sum(diag(LamLam’ + Psi)) where Lam is a p x k matrix of environment loadings (the rr..fa1 and rr..fa2 variance parameters) and Psi is a p x p diagonal matrix with diagonal elements given by the sigma_ge_j for the p environments

FAmodel_VAF <- function(model, k) {
    # model <- asr5.A
    # model <- asr6.A
    Lam <- matrix(summary(model)$varcom[grep("rr.*fa.$", rownames(summary(model)$varcom)), "component"], byrow = FALSE, ncol = k)
    Psi <- diag(summary(model)$varcom[grep("^YEAR_LOCATION:vm\\(LINE_ID, ", rownames(summary(model)$varcom)), "component"])
    p <- nlevels(data5$YEAR_LOCATION) # 40
    # v_g <- summary(model)$varcom[1, "component"]
    # v_ge <- summary(model)$varcom[2, "component"]
    round(sum(diag(Lam %*% t(Lam))) / sum(diag(Lam %*% t(Lam) + Psi)) * 100, 2)
}

FAmodel_VAF(asr5.A, 1) # 49.57
FAmodel_VAF(asr5.G, 1) # 55.23
FAmodel_VAF(asr5.T, 1) # 46.33

FAmodel_VAF(asr6.A, 2) # 67.51
FAmodel_VAF(asr6.G, 2) # 68.86
FAmodel_VAF(asr6.T, 2) # 64.46

FAmodel_VAF(asr7.A, 3) # 81.79
FAmodel_VAF(asr7.G, 3) # 78.41
FAmodel_VAF(asr7.T, 3) # 80.21
