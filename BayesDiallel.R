rm(list = ls())
library("BayesDiallel")

########################################################################################
## Diallel Example Analyses
##

################################################################################
## Load in R packages
##
##   These R packages must be installed from CRAN to use BayesDiallel to its fullest.
try(require(lattice, quietly=TRUE, warn.conflicts=FALSE), silent=TRUE);
try(require(coda, quietly=TRUE, warn.conflicts=FALSE), silent=TRUE);
try(require(R.methodsS3, quietly=TRUE, warn.conflicts=FALSE), silent=TRUE);
try(require(R.oo, quietly=TRUE, warn.conflicts=FALSE), silent=TRUE);
try(require(corpcor, quietly=TRUE, warn.conflicts=FALSE), silent=TRUE);
require(BayesDiallel, quietly = TRUE, warn.conflicts=FALSE);


################################################################################
##  First example, examination of "Piximis Mouse Data"
##
##  The data set csv files for inspection of phenotypes for some Diallel data from
## crosses of Collaborative Cross founder strains are within the "Examples"
## directory of this package.  Loading "data(PiximusData)"  automatically loads
## those files into the R terminal.  The dataset is divided into the data file "Piximus.Data"
## the Models File "Examples.Piximus.Models" and the Bayes priors file "Examples.Piximus.tau.Prior.Info".
##
##
##  Length and weight phenotypes can be analyzed using default methods for 
## function "DiallelAnalyzer()"
##

###############################################################################
## First Example Piximus Diallel Mouse Data.
data(PiximusData);
##  Data is in "Piximus.Data"
##  Models are in "Examples.Piximus.Models"
##  tau.Prior is in "Examples.Piximus.tau.Prior.Info"
####################################
###### Run the Algorithm->
AFD = DiallelAnalyzer(data = Piximus.Data, father.strain="father.strain.name",
                      mother.strain="mother.strain.name", phenotype="MouseWeight", is.female="is.female",
                      sep="", FixedEffects = NULL,
                      RandomEffects = NULL, 
                      Models=Example.Piximus.Models[1], sigmasq.start = 1,  numChains = 5,  lengthChains=2500,
                      burnin = 1,
                      DIC.Only=FALSE,  tauPriorFile = Example.Piximus.tau.Prior.Info,
                      SaveAFDFile = "SaveAFDBackUpPiximusData.RData",
                      LogTransform=TRUE);
save(AFD=AFD, this=AFD, file="SaveAFDBackUpPiximusData.RData");

AFDPreCentered = DiallelAnalyzer(data = Piximus.Data, father.strain="father.strain.name",
                                 mother.strain="mother.strain.name", phenotype="MouseWeight", is.female="is.female",
                                 sep="", FixedEffects = NULL,
                                 RandomEffects = NULL, 
                                 Models=Example.Piximus.Models[1:2], sigmasq.start = 1,  numChains = 5,  lengthChains=2500,
                                 burnin = 1,
                                 DIC.Only=FALSE,  tauPriorFile = Example.Piximus.tau.Prior.Info,
                                 SaveAFDFile = "SaveAFDBackUpPiximusDataCentered.RData",
                                 LogTransform=TRUE, DoFirstCenter =TRUE, Verbose = 1);

## Load pregenerated AFD from backup
if (FALSE) {
  library(R.oo);  library(BayesDiallel);
  load("SaveAFDBackUpPiximusData.RData");
  
  library(R.oo);  library(BayesDiallel);
  load("SaveAFDBackUpPiximusDataCentered.RData");
}

##List of Elements available to AFD
AFD$ListElements

##List of Functions available to AFD
AFD$ListFunctions;

## Get a summary of the Coda Chains for first Diallel Model
summary(AFD$AllDiallelObs[[1]]$raw.chains)
plot(AFD$AllDiallelObs[[1]]$raw.chains);
plot.hpd(AFD$AllDiallelObs[[1]]$raw.chains);

## Now make those TwoDiallelPlots!
## Plot Females, observed versus fitted
AFD$AllDiallelObs[[1]]$TwoDiallelPlot(FemalePlot = TRUE, MaleAgainstFemale = FALSE);

## PlotStrawPlot illustrates effects
AFD$AllDiallelObs[[1]]$PlotStrawPlot();
library(BayesDiallel);
data(PiximusData);
dir.create("BS");
AFDPreCenteredBMCwithWeight = DiallelAnalyzer(data = Piximus.Data, father.strain="father.strain.name",
                                              mother.strain="mother.strain.name", phenotype="BoneMineralContent", is.female="is.female",
                                              sep="", FixedEffects = "MouseWeight",
                                              RandomEffects = NULL, 
                                              Models=Example.Piximus.Models[1], sigmasq.start = 1,  numChains = 3,  lengthChains=900,
                                              burnin = 1,
                                              DIC.Only=FALSE,  tauPriorFile = Example.Piximus.tau.Prior.Info,
                                              SaveAFDFile = "SaveAFDBackUpPiximusDataWeightCentered.RData",
                                              DoBayesSpike=TRUE, BSSaveDir ="BS",
                                              LogTransform=FALSE, DoFirstCenter =TRUE, PriorProbFixed = c(-1,-1,-1,-1,1.0), 
                                              Verbose = 1, BSLengthChains = 400, BSNumChains =2);
AFDPreCenteredBMCMaybeWeight = DiallelAnalyzer(data = Piximus.Data, father.strain="father.strain.name",
                                               mother.strain="mother.strain.name", phenotype="BoneMineralContent", is.female="is.female",
                                               sep="", FixedEffects = "MouseWeight",
                                               RandomEffects = NULL, 
                                               Models=Example.Piximus.Models[1:2], sigmasq.start = 1,  numChains = 2,  lengthChains=500,
                                               burnin = 1,
                                               DIC.Only=FALSE,  tauPriorFile = Example.Piximus.tau.Prior.Info,
                                               SaveAFDFile = "SaveAFDBackUpPiximusDataCentered.RData",
                                               DoBayesSpike=TRUE,  BSSaveDir ="BS",
                                               LogTransform=FALSE, DoFirstCenter =TRUE, PriorProbFixed = c(-1,-1,-1,-1,-1), Verbose = 1,
                                               BSLengthChains = 500, BSNumChains =2);


if (FALSE) {
  library(BayesDiallel);
  load("SaveAFDBackUpPiximusDataCentered.RData");
}
AFD <-  AFDPreCenteredBMCwithWeight;
AFD$AllDiallelObs[[1]]$PosteriorPSqHPDTable
AFD$AllDiallelObs[[1]]$PosteriorPSqHPDSubTable
AFD$AllDiallelObs[[1]]$PosteriorDSqHPDTable
AFD$AllDiallelObs[[1]]$PosteriorDSqHPDSubTable

AFD$.AllFixedVariables;
AFD$.AllRandomVariables;
BS1 <-AFDPreCenteredBMCwithWeight$BS
BS1$CodaList[[1]][1:5, 165:181]
BS2 <- AFDPreCenteredBMCMaybeWeight$BS 

ART <- BS2$XtResid - t(BS2$X) 

AXX <- t(BS2$X) 
sum(abs(AXX - BS2$XtX))

BS1$CodaList[[1]][1:10, BS1$OtherNameCodaList == "Prob:Fixed:FixedEffect:1"]   
BS2$CodaList[[1]][1:10, BS2$OtherNameCodaList == "Prob:Fixed:FixedEffect:1"]  


###############################################################################
## Simulate Fake data and experiment with fitting
##############################################################################

TrueModel = c("BSasmvb");
WalkAndReadCrossList(TrueModel)  
## Set the model shape to simulate
ModelVector = c(1,1,2,1,3,1);
ModelVector = WalkAndReadCrossList(TrueModel);
ajModel = ModelVector[1];
MotherModel = ModelVector[2];
CrossModel = ModelVector[3];
BetaInbredLevel = ModelVector[4];   ## 09062014 we rename "BetaHybridLevel"
##  to "BetaInbredLevel" for clarity
SexModel = ModelVector[5];
BetaInbredTimesSex = ModelVector[6];
IsFixed = 1; IsRandom = 1;
NumReps = 2 ## Imagine dataset of size 2 times full diallel
numj = 8 ## 8 strains of mice?


print("Sampling NonDiallel Random and Fake Effects");
InFixedEffects=NULL; InRandomEffects =NULL;
if (IsFixed || IsRandom) {
  WantLen = (numj *(numj+1)/2 ) *NumReps ;
  InRandomEffects = sample(size=WantLen, 1:5, replace=TRUE);
  InFixedEffects = rnorm(WantLen);
  if (IsFixed == 0) {InFixedEffects=NULL;}
  if (IsRandom == 0) {InRandomEffects=NULL;}
  mRandomEffects = 20;
  dfRandomEffects = 5; tauFixedEffects = 20;
  NumRepetitions=1;
}  else {
  NumRepetitions = NumReps;
}

## Simulate a data set given model of choice
print("Making first Simulation");
### Allocate AFDS Diallel Simulation Object
AFDS <- FullDiallelSimulate(numj = numj, InFixedEffects = InFixedEffects,
                            InRandomEffects = InRandomEffects, mRandomEffects = 20,
                            dfRandomEffects = 5, tauFixedEffects = 20, YesSex = TRUE, NumRepetitions = NumRepetitions);

## Sample the parameters for Diallel Simulation Object
defaultSimulateParam(AFDS);
if (is.null(InRandomEffects) && is.null(InFixedEffects)) {
  SetupFilledSimulatedListjk(AFDS, NumRepetitions = NumRepetitions);
}


## Set up the model for Diallel Simulation Object
SetupSimulateMethod(AFDS, ajModel = 1, MotherModel = MotherModel, CrossModel =
                      CrossModel, SexModel = SexModel, BetaInbredLevel = BetaInbredLevel, 
                    BetaInbredTimesSex = BetaInbredTimesSex);
## SimValues now simulates a random Beta vector, tau vector, and observed data Y      
SimValues(AFDS, Sigma = 1);

#############################################################################
###  Begin Analysis
###
### From Simulated Diallel in AFDS, we pull observed data to put in an 
###  Analysis Diallel AFD

Y = AFDS$Y; ListjkInput = AFDS$Listjk;
SexInputVec = AFDS$SexVector

## Pull Out the RandomEffects
ARandomEffects = rep(0, length(AFDS$RandomEffects[,1]));
for (ii in 1:length(AFDS$RandomEffects[1,])) {
  ARandomEffects[AFDS$RandomEffects[,ii] == 1] = ii;
}
if (is.null(InRandomEffects)) { ARandomEfects = NULL; }
##TargetLength=TARGET;
SampleBeta = AFDS$ADiallelOb$Beta;
FileDir = AllDir2;
quantile = .95;
INTIN = 0;

## These are the true simulation parameters from AFDS we would like experiment to learn  
TrueParametersList = c(AFDS$ADiallelOb$Beta, AFDS$ADiallelOb$tau,
                       AFDS$ADiallelOb$Sigma)

Models =  c(TrueModel, "BSadf3s", "Bvswabsaas");


########################################################
##  Call DiallelAnalyze
##
MyDataFileTable = cbind(Y, ListjkInput, SexInputVec, AFDS$FixedEffects, ARandomEffects);
colnames(MyDataFileTable) = c("phenotype", "mother.strain", "father.strain", 
                              "is.female", colnames(AFDS$FixedEffects), colnames(ARandomEffects));
AFD = DiallelAnalyzer(DataFile = MyDataFileTable, father.strain="father.strain",
                      mother.strain="mother.strain", phenotype="phenotype", is.female="is.female",
                      FixedEffects = colnames(AFDS$FixedEffects),
                      RandomEffects = ARandomEffects, 
                      Models=Models, sigmasq.start = 1;  numChains = 5;  lengthChains=2500,burnin = 1,
                      DIC.Only=FALSE,  TauFile = NULL);

## List of Elements Available to AFD
AFD$ListElements;
## List of Functions Available to AFD
AFD$ListFunctions;

ADO = AFD$AllDiallelObs[[1]];  ## This is the DiallelOb 
## corresponding to our true model

plot(ADO$raw.chains);  ## Plot the chains to understand
summary(ADO$cent.chains);  ## Get Summary information on centered parameters
plot(ADO$cent.chains);  ## Plot the chains to understand
summary(ADO$raw.chains);  ## Get Summary information on a parameters

plot.hpd(ADO$raw.chains); ## Using Library_diallel_helper, plot chains.

plot.corr(crosscorr(ADO$raw.chains)) ## Look at cross correlations of coda chains


CalcDICC(ADO, AFD) ## DIC value for this Coda Chain  

## This creates a summary table of applicable elements
MySum = summary.table(ADO, burnin = 400, AFD=AFD);

#MyPostSumGivesValues for all the diallel squares of every gender
MyPostSum = PosteriorPredSummary(ADO, AFD, burnin = 400, AFD=AFD, keep = TRUE);

SimmedGrid = GiveGrid(AFDS);
par(mfrow=c(1,2));
draw.diallel(data =SimmedGrid[MyPostSum$"Female/Male" == 1,], phenotype = "Mean.Simulated")
draw.diallel(data = MyPostSum[MyPostSum$"Female/Male" == 1,], phenotype = "Mean Posterior")



