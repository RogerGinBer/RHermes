## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, warning=FALSE, message = FALSE)
knitr::opts_chunk$set(collapse = TRUE)
library(plotly)

## ---- eval = FALSE------------------------------------------------------------
#  if(!requireNamespace("devtools", quietly = TRUE))
#      install.packages("devtools")
#  devtools::install_github("RogerGinBer/RHermes")

## ---- message=FALSE-----------------------------------------------------------
library(RHermes)
myHermes <- RHermesExp()

## -----------------------------------------------------------------------------
myHermes <- setExpParam(myHermes, params = ExpParam(ppm = 4, res = 120000,
                                                    ion = "+"))
myHermes <- setExpParam(myHermes, template = "orbi-pos")

## -----------------------------------------------------------------------------
myHermes <- setDB(myHermes, db = "hmdb") #Loads a small HMDB subset

## -----------------------------------------------------------------------------
myHermes <- setDB(myHermes, db = "hmdb", admult = 2, adcharge = 2)

## -----------------------------------------------------------------------------
#We will manually add the "M+2H" adduct
myHermes <- addAd(myHermes, name = "M+2H", deltam = 2*1.00727600, ch = 2,
                    mult = 1, toadd = "H2")

#For instance, remove adducts of unused solvents
myHermes <- remAd(myHermes, c("M+DMSO+H", "M+IsoProp+H"))


## ----eval=FALSE---------------------------------------------------------------
#  adlist(myHermes)

## ----echo=FALSE---------------------------------------------------------------
knitr::kable(adlist(myHermes))

## ---- message = TRUE----------------------------------------------------------
myHermes #Shows a summary of all set parameters and object info

## ---- echo = FALSE------------------------------------------------------------
#For speed with a small dataset
myHermes <- setCluster(myHermes, BiocParallel::SerialParam()) 

## -----------------------------------------------------------------------------
myHermes <- processMS1(myHermes,
                        system.file("ExtData/MS1TestData.mzML",
                                    package = "RHermes"))

## ---- message = TRUE----------------------------------------------------------
#You could do either of these, but PL it's more direct if you just want to check
#one file:
# myHermes
PL(myHermes, 1)

## -----------------------------------------------------------------------------
s <- getSOIpar("double")
myHermes <- findSOI(myHermes, s, 1)

## ----eval=FALSE---------------------------------------------------------------
#  #Generate a SOI list of file 2 without blank subtraction and another one
#  #using file 1 as blank
#  myHermes <- findSOI(myHermes, s, c(2,2), c(0,1))

## -----------------------------------------------------------------------------
myHermes <- filterSOI(myHermes, id = 1, minint = 20000, isofidelity = TRUE)

## ---- message = TRUE----------------------------------------------------------
#Only uses SOIs of adducts in "ad". In this case, only M+H adducts
myHermes <- generateIL(myHermes, id = 1, par = ILParam(priorization = "only", ad = "M+H"))

#If there's an M+H adduct annotation, don't include other adducts.
#Only works if you've performed filterSOI first, since it requires adduct 
#annotation grouping
myHermes <- generateIL(myHermes, id = 1, par = ILParam(priorization = "yes", ad = "M+H"))

#Use all SOIs in the SOI list
myHermes <- generateIL(myHermes, id = 1, par = ILParam(priorization = "all"))

## -----------------------------------------------------------------------------
#This removes all IL entries below 50000 intensity in the 0s-150s RT interval
myHermes <- filterIL(myHermes, 1, rts = c(0,150), minint = 5e4)

## ---- eval = FALSE------------------------------------------------------------
#  exportIL(myHermes, id = 1, file = "InclusionList", maxOver = 5, sepFiles = TRUE)

## ---- eval = FALSE------------------------------------------------------------
#  myHermes <- processMS2(myHermes, id = 1,
#                         MS2files = c("./file1.mzML", "./file2.mzML",
#                                      "./file3.mzML", "./file4.mzML"),
#                         sstype = "regular", useDB = FALSE)

## -----------------------------------------------------------------------------
myHermes <- readRDS(system.file("extdata", "exampleObject.rds",
                                package = "RHermes"))

## -----------------------------------------------------------------------------
Ident(myHermes, 1)

## ---- eval = FALSE------------------------------------------------------------
#  saveRDS(myHermes, "testRHermes.rds")

## ---- eval = FALSE------------------------------------------------------------
#  myHermes <- readRDS("testRHermes.rds")

