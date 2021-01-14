---
title: "1st Andrea Brunner Workflow"
output:
  html_notebook: default
  html_document:
    df_print: paged
    theme: lumen
    toc: true
    toc_float: true
    toc_depth: 4
  pdf_document: default
vignette: >
  %\VignetteIndexEntry{1st Andrea Brunner Workflow}
  %\VignetteEngine{knitr::rmarkdown}
  \usepackage[utf8]{inputenc}
---

# Context

This is the second workflow corresponding to Andrea Brunner's analysis.
We will use the new RHermes S4 classes and methods to illustrate its use.

# MS1 Processing

## Setting the parameters

The first step is to generate an RHermesExp object and set the MS1 experimental 
parameters. To ease the job for the end-user, 4 different templates are
available for both Q-TOF and Orbitrap instruments. The user can also specify 
custom parameters using ExpParam() instead of the template.

New templat




```r
# setwd("../") #For unknitted use of the chunk
myHermes <- RHermesExp()
myHermes <- setExpParam(myHermes,
                        params = ExpParam(ppm = 3.5, res = 120000,
                                          instr = "Orbitrap", minmz = 50,
                                          maxmz = 1200))

negHermes <- RHermesExp()
negHermes <- setExpParam(negHermes,
                         params = ExpParam(ppm = 3.5, res = 120000,
                                          instr = "Orbitrap", minmz = 50,
                                          maxmz = 1200, ion = "-"))
```

Processing the RHermesExp object will set timestamps of the performed 
operations. You can use readTime(object) to see them.

```r
readTime(myHermes)
```

```
## System info: Windows/10 x64/build 18363/DESKTOP-P7QILPA/x86-64/Roger/Roger/Roger
```

Now we will set the formulas and adducts to consider in our experiment. The
default database is the Human Metabolome Database (HMDB), but NORMAN is also 
available and is particularly useful for environmental pollutant profiling.

The IDs to use for the two databases are, respectively, "hmdb", and "norman".

RHermes directly checks the polarity used in the experiment (set in the 
previous step) and selects the corresponding adducts.

The user can also select the characteristics of the adducts considered, in 
particular one can select the maximum charge (+, 2+, etc.) and the maximume
multiplicity of the adducts (M, 2M, 3M, etc.). The defaults are 1 for both
charge and multiplicity.


```r
myHermes <- setDB(myHermes, db = "custom", filename = "D:/norman.xls")
myHermes <- remAd(myHermes, myHermes@metadata@ExpParam@adlist$adduct[-c(1:5)])
```

```
## This is the new adduct list:
```

```r
negHermes <- setDB(negHermes, "custom", filename = "D:/norman.xls")
```

```
## Parsing the custom formula database
```

```r
negHermes <- remAd(negHermes, negHermes@metadata@ExpParam@adlist$adduct[-c(1,7)])
readTime(myHermes)
```

## File processing into PLs

Once the basic parameters have been set it is time to process some files. To do
this we just need to use the function FileProc with the names of the files we 
want processed. There is the possibility to process labelled data so that all 
possible C isotopes are searched but, in that case, the labelled files must be 
processed in a separate function call (or just process everything in _labelled 
mode_, just bear in mind it is significantly slower).


```r
dir <- "D:/ABrunner"
fil <- list.files(dir, pattern = ".*pos.*.mzML", full.names = TRUE)[c(1,4,7)]
myHermes <- FileProc(myHermes, files = fil)
```

```
## Preprocessing steps, calculating possible ionic formulas and isotopic distributions
```

```
##   |                                                                                                                                   |                                                                                                                           |   0%  |                                                                                                                                   |===============================                                                                                            |  25%  |                                                                                                                                   |==============================================================                                                             |  50%  |                                                                                                                                   |============================================================================================                               |  75%  |                                                                                                                                   |===========================================================================================================================| 100%
```

```
##  done.
```

```
##   |                                                                                                                                   |                                                                                                                           |   0%  |                                                                                                                                   |===============================                                                                                            |  25%  |                                                                                                                                   |==============================================================                                                             |  50%  |                                                                                                                                   |============================================================================================                               |  75%  |                                                                                                                                   |===========================================================================================================================| 100%
```

```
##   |                                                                                                                                   |                                                                                                                           |   0%  |                                                                                                                                   |===============================                                                                                            |  25%  |                                                                                                                                   |==============================================================                                                             |  50%  |                                                                                                                                   |============================================================================================                               |  75%  |                                                                                                                                   |===========================================================================================================================| 100%
```

```
##   |                                                                                                                                   |                                                                                                                           |   0%  |                                                                                                                                   |===============================                                                                                            |  25%  |                                                                                                                                   |==============================================================                                                             |  50%  |                                                                                                                                   |============================================================================================                               |  75%  |                                                                                                                                   |===========================================================================================================================| 100%
```

```
##   |                                                                                                                                   |                                                                                                                           |   0%  |                                                                                                                                   |===============================                                                                                            |  25%  |                                                                                                                                   |==============================================================                                                             |  50%  |                                                                                                                                   |============================================================================================                               |  75%  |                                                                                                                                   |===========================================================================================================================| 100%
```

```r
fil <- list.files(dir, pattern = ".*neg.*.mzML", full.names = TRUE)[c(1,4,7)]
negHermes <- FileProc(negHermes, files = fil)
```

```
##   |                                                                                                                                   |                                                                                                                           |   0%  |                                                                                                                                   |===============================                                                                                            |  25%  |                                                                                                                                   |==============================================================                                                             |  50%  |                                                                                                                                   |============================================================================================                               |  75%  |                                                                                                                                   |===========================================================================================================================| 100%
```

```
##  done.
```

```
##   |                                                                                                                                   |                                                                                                                           |   0%  |                                                                                                                                   |===============================                                                                                            |  25%  |                                                                                                                                   |==============================================================                                                             |  50%  |                                                                                                                                   |============================================================================================                               |  75%  |                                                                                                                                   |===========================================================================================================================| 100%
```

```
##   |                                                                                                                                   |                                                                                                                           |   0%  |                                                                                                                                   |===============================                                                                                            |  25%  |                                                                                                                                   |==============================================================                                                             |  50%  |                                                                                                                                   |============================================================================================                               |  75%  |                                                                                                                                   |===========================================================================================================================| 100%
```

```
##   |                                                                                                                                   |                                                                                                                           |   0%  |                                                                                                                                   |===============================                                                                                            |  25%  |                                                                                                                                   |==============================================================                                                             |  50%  |                                                                                                                                   |============================================================================================                               |  75%  |                                                                                                                                   |===========================================================================================================================| 100%
```

```
##   |                                                                                                                                   |                                                                                                                           |   0%  |                                                                                                                                   |===============================                                                                                            |  25%  |                                                                                                                                   |==============================================================                                                             |  50%  |                                                                                                                                   |============================================================================================                               |  75%  |                                                                                                                                   |===========================================================================================================================| 100%
```
## SOI generation

Once every file has been processed and all signals that match your F/A
combinations have been registered into _Peak Lists_ the next step will be SOI
detection.

As described in the paper, RHermes detects SOI (_Scans of Interest_) by applying
one or multiple filters to the data points registered in a PL. To provide an 
almost hands-off experience, we have designed multiple filter templates that can
be easily accessed and customized. They are:

* _single_ and _single-x_
* _double_ and _double-x_
* _triple_ and _triple-x_

Where the "-x" suffix means "extended" and is to be used when dealing with long
chromatographic runs where compounds can elute during a long time. Double and
triple filter modes use staggered, partially-overlapping time segments to avoid
signal splitting.

In case you need more filters or if you just want to customize them you can find
the "SOITemplates.csv" file at /app/www in the RHermes package folder. Just add
more lines with your desired parameters using the same ID in the first column
and it will be available to getSOIpar().

It is likely that you may want to test some filtering parameters before
commiting to perform a set of MSMS runs. With this in mind, SOIfinder() allows
the input of multiple parameter objects in the form of a list (or c()). One 
 thing to keep in mind is that, in the case the length of parameters is
shorter than the number of SOI lists to generate, **only** the last parameter of
the list will be recycled.

Another feature that RHermes offers is the possibility to select a PL as
blank and _substract_ those signals when processing another PL. This step is 
optional but really advisable, since it will reduce the number of MSMS 
injections needed afterwards. It works like this:

 * You can input a numeric vector telling the program which PLs you want to use
 as blank. 
 * If you don't enter anything it doesn't perform any substraction by
 default, so don't worry.
 * 0 means _No blank substraction_, and it is important because both the 
 _fileID_ and _against_ vectors will be matched **one-to-one**.
 



```r
setwd("../") #For unknitted use of the chunk
s <- getSOIpar("double")
myHermes <- SOIfinder(myHermes, params = s, fileID = c(1,2,2,3,3,3),
                      against = c(0,0,1,0,1,2))
```

```
## Starting density filtering:
```

```
## Running Density Filter
```

```
##   |                                                                                                                                   |                                                                                                                           |   0%  |                                                                                                                                   |===============================                                                                                            |  25%  |                                                                                                                                   |==============================================================                                                             |  50%  |                                                                                                                                   |============================================================================================                               |  75%  |                                                                                                                                   |===========================================================================================================================| 100%
```

```
## Running Density Interpreter
```

```
##   |                                                                                                                                   |                                                                                                                           |   0%  |                                                                                                                                   |===============================                                                                                            |  25%  |                                                                                                                                   |==============================================================                                                             |  50%  |                                                                                                                                   |============================================================================================                               |  75%  |                                                                                                                                   |===========================================================================================================================| 100%
```

```
## Running Density Filter
```

```
##   |                                                                                                                                   |                                                                                                                           |   0%  |                                                                                                                                   |===============================                                                                                            |  25%  |                                                                                                                                   |==============================================================                                                             |  50%  |                                                                                                                                   |============================================================================================                               |  75%  |                                                                                                                                   |===========================================================================================================================| 100%
```

```
## Running Density Interpreter
```

```
##   |                                                                                                                                   |                                                                                                                           |   0%  |                                                                                                                                   |===============================                                                                                            |  25%  |                                                                                                                                   |==============================================================                                                             |  50%  |                                                                                                                                   |============================================================================================                               |  75%  |                                                                                                                                   |===========================================================================================================================| 100%
```

```
## Now Grouping:
```

```
## Merging filter 2 out of 2 with previous SOI list
```

```
##   |                                                                                                                                   |                                                                                                                           |   0%  |                                                                                                                                   |===============================                                                                                            |  25%  |                                                                                                                                   |==============================================================                                                             |  50%  |                                                                                                                                   |============================================================================================                               |  75%  |                                                                                                                                   |===========================================================================================================================| 100%
```

```
## Adding entries unique from the first DF
```

```
## Adding entries unique from the second DF
```

```
## Solving inner conflicts, Round: 1
```

```
## Solving inner conflicts, Round: 2
```

```
## Solving inner conflicts, Round: 3
```

```
## Initial peak retrieval:
```

```
##   |                                                                                                                                   |                                                                                                                           |   0%  |                                                                                                                                   |===============================                                                                                            |  25%  |                                                                                                                                   |==============================================================                                                             |  50%  |                                                                                                                                   |============================================================================================                               |  75%  |                                                                                                                                   |===========================================================================================================================| 100%
```

```
## 
```

```
## Starting group characterization:
```

```
## Mass calculation:
```

```
## Shortening and selecting long groups:
```

```
##   |                                                                                                                                   |                                                                                                                           |   0%  |                                                                                                                                   |===============================                                                                                            |  25%  |                                                                                                                                   |==============================================================                                                             |  50%  |                                                                                                                                   |============================================================================================                               |  75%  |                                                                                                                                   |===========================================================================================================================| 100%
```

```
## 
```

```
## Width calculation:
```

```
## Maximum intensity calculation:
```

```
## Number of scans:
```

```
## Converting from ionic formula to F/A combinations:
```

```
## Generating peaklist for plotting:
```

```
##   |                                                                                                                                   |                                                                                                                           |   0%  |                                                                                                                                   |===============================                                                                                            |  25%  |                                                                                                                                   |==============================================================                                                             |  50%  |                                                                                                                                   |============================================================================================                               |  75%  |                                                                                                                                   |===========================================================================================================================| 100%
```

```
## Starting density filtering:
```

```
## Running Density Filter
```

```
##   |                                                                                                                                   |                                                                                                                           |   0%  |                                                                                                                                   |===============================                                                                                            |  25%  |                                                                                                                                   |==============================================================                                                             |  50%  |                                                                                                                                   |============================================================================================                               |  75%  |                                                                                                                                   |===========================================================================================================================| 100%
```

```
## Running Density Interpreter
```

```
##   |                                                                                                                                   |                                                                                                                           |   0%  |                                                                                                                                   |===============================                                                                                            |  25%  |                                                                                                                                   |==============================================================                                                             |  50%  |                                                                                                                                   |============================================================================================                               |  75%  |                                                                                                                                   |===========================================================================================================================| 100%
```

```
## Running Density Filter
```

```
##   |                                                                                                                                   |                                                                                                                           |   0%  |                                                                                                                                   |===============================                                                                                            |  25%  |                                                                                                                                   |==============================================================                                                             |  50%  |                                                                                                                                   |============================================================================================                               |  75%  |                                                                                                                                   |===========================================================================================================================| 100%
```

```
## Running Density Interpreter
```

```
##   |                                                                                                                                   |                                                                                                                           |   0%  |                                                                                                                                   |===============================                                                                                            |  25%  |                                                                                                                                   |==============================================================                                                             |  50%  |                                                                                                                                   |============================================================================================                               |  75%  |                                                                                                                                   |===========================================================================================================================| 100%
```

```
## Now Grouping:
```

```
## Merging filter 2 out of 2 with previous SOI list
```

```
##   |                                                                                                                                   |                                                                                                                           |   0%  |                                                                                                                                   |===============================                                                                                            |  25%  |                                                                                                                                   |==============================================================                                                             |  50%  |                                                                                                                                   |============================================================================================                               |  75%  |                                                                                                                                   |===========================================================================================================================| 100%
```

```
## Adding entries unique from the first DF
```

```
## Adding entries unique from the second DF
```

```
## Solving inner conflicts, Round: 1
```

```
## Solving inner conflicts, Round: 2
```

```
## Solving inner conflicts, Round: 3
```

```
## Initial peak retrieval:
```

```
##   |                                                                                                                                   |                                                                                                                           |   0%  |                                                                                                                                   |===============================                                                                                            |  25%  |                                                                                                                                   |==============================================================                                                             |  50%  |                                                                                                                                   |============================================================================================                               |  75%  |                                                                                                                                   |===========================================================================================================================| 100%
```

```
## 
```

```
## Starting group characterization:
```

```
## Mass calculation:
```

```
## Shortening and selecting long groups:
```

```
##   |                                                                                                                                   |                                                                                                                           |   0%  |                                                                                                                                   |===============================                                                                                            |  25%  |                                                                                                                                   |==============================================================                                                             |  50%  |                                                                                                                                   |============================================================================================                               |  75%  |                                                                                                                                   |===========================================================================================================================| 100%
```

```
## 
```

```
## Width calculation:
```

```
## Maximum intensity calculation:
```

```
## Number of scans:
```

```
## Converting from ionic formula to F/A combinations:
```

```
## Generating peaklist for plotting:
```

```
##   |                                                                                                                                   |                                                                                                                           |   0%  |                                                                                                                                   |===============================                                                                                            |  25%  |                                                                                                                                   |==============================================================                                                             |  50%  |                                                                                                                                   |============================================================================================                               |  75%  |                                                                                                                                   |===========================================================================================================================| 100%
```

```
## Starting density filtering:
```

```
## Running Density Filter
```

```
##   |                                                                                                                                   |                                                                                                                           |   0%  |                                                                                                                                   |===============================                                                                                            |  25%  |                                                                                                                                   |==============================================================                                                             |  50%  |                                                                                                                                   |============================================================================================                               |  75%  |                                                                                                                                   |===========================================================================================================================| 100%
```

```
## Running Density Interpreter
```

```
##   |                                                                                                                                   |                                                                                                                           |   0%  |                                                                                                                                   |===============================                                                                                            |  25%  |                                                                                                                                   |==============================================================                                                             |  50%  |                                                                                                                                   |============================================================================================                               |  75%  |                                                                                                                                   |===========================================================================================================================| 100%
```

```
## Running Density Filter
```

```
##   |                                                                                                                                   |                                                                                                                           |   0%  |                                                                                                                                   |===============================                                                                                            |  25%  |                                                                                                                                   |==============================================================                                                             |  50%  |                                                                                                                                   |============================================================================================                               |  75%  |                                                                                                                                   |===========================================================================================================================| 100%
```

```
## Running Density Interpreter
```

```
##   |                                                                                                                                   |                                                                                                                           |   0%  |                                                                                                                                   |===============================                                                                                            |  25%  |                                                                                                                                   |==============================================================                                                             |  50%  |                                                                                                                                   |============================================================================================                               |  75%  |                                                                                                                                   |===========================================================================================================================| 100%
```

```
## Now Grouping:
```

```
## Merging filter 2 out of 2 with previous SOI list
```

```
##   |                                                                                                                                   |                                                                                                                           |   0%  |                                                                                                                                   |===============================                                                                                            |  25%  |                                                                                                                                   |==============================================================                                                             |  50%  |                                                                                                                                   |============================================================================================                               |  75%  |                                                                                                                                   |===========================================================================================================================| 100%
```

```
## Adding entries unique from the first DF
```

```
## Adding entries unique from the second DF
```

```
## Solving inner conflicts, Round: 1
```

```
## Solving inner conflicts, Round: 2
```

```
## Solving inner conflicts, Round: 3
```

```
## Initial peak retrieval:
```

```
##   |                                                                                                                                   |                                                                                                                           |   0%  |                                                                                                                                   |===============================                                                                                            |  25%  |                                                                                                                                   |==============================================================                                                             |  50%  |                                                                                                                                   |============================================================================================                               |  75%  |                                                                                                                                   |===========================================================================================================================| 100%
```

```
## 
```

```
## Starting group characterization:
```

```
## Mass calculation:
```

```
## Shortening and selecting long groups:
```

```
##   |                                                                                                                                   |                                                                                                                           |   0%  |                                                                                                                                   |===============================                                                                                            |  25%  |                                                                                                                                   |==============================================================                                                             |  50%  |                                                                                                                                   |============================================================================================                               |  75%  |                                                                                                                                   |===========================================================================================================================| 100%
```

```
## Blank substraction:
```

```
##   |                                                                                                                                   |                                                                                                                           |   0%  |                                                                                                                                   |===============================                                                                                            |  25%  |                                                                                                                                   |==============================================================                                                             |  50%  |                                                                                                                                   |============================================================================================                               |  75%  |                                                                                                                                   |===========================================================================================================================| 100%
## 
##   |                                                                                                                                   |                                                                                                                           |   0%  |                                                                                                                                   |===============================                                                                                            |  25%  |                                                                                                                                   |==============================================================                                                             |  50%  |                                                                                                                                   |============================================================================================                               |  75%  |                                                                                                                                   |===========================================================================================================================| 100%
```

```
## 
```

```
## Width calculation:
```

```
## Maximum intensity calculation:
```

```
## Number of scans:
```

```
## Converting from ionic formula to F/A combinations:
```

```
## Generating peaklist for plotting:
```

```
##   |                                                                                                                                   |                                                                                                                           |   0%  |                                                                                                                                   |===============================                                                                                            |  25%  |                                                                                                                                   |==============================================================                                                             |  50%  |                                                                                                                                   |============================================================================================                               |  75%  |                                                                                                                                   |===========================================================================================================================| 100%
```

```
## Starting density filtering:
```

```
## Running Density Filter
```

```
##   |                                                                                                                                   |                                                                                                                           |   0%  |                                                                                                                                   |===============================                                                                                            |  25%  |                                                                                                                                   |==============================================================                                                             |  50%  |                                                                                                                                   |============================================================================================                               |  75%  |                                                                                                                                   |===========================================================================================================================| 100%
```

```
## Running Density Interpreter
```

```
##   |                                                                                                                                   |                                                                                                                           |   0%  |                                                                                                                                   |===============================                                                                                            |  25%  |                                                                                                                                   |==============================================================                                                             |  50%  |                                                                                                                                   |============================================================================================                               |  75%  |                                                                                                                                   |===========================================================================================================================| 100%
```

```
## Running Density Filter
```

```
##   |                                                                                                                                   |                                                                                                                           |   0%  |                                                                                                                                   |===============================                                                                                            |  25%  |                                                                                                                                   |==============================================================                                                             |  50%  |                                                                                                                                   |============================================================================================                               |  75%  |                                                                                                                                   |===========================================================================================================================| 100%
```

```
## Running Density Interpreter
```

```
##   |                                                                                                                                   |                                                                                                                           |   0%  |                                                                                                                                   |===============================                                                                                            |  25%  |                                                                                                                                   |==============================================================                                                             |  50%  |                                                                                                                                   |============================================================================================                               |  75%  |                                                                                                                                   |===========================================================================================================================| 100%
```

```
## Now Grouping:
```

```
## Merging filter 2 out of 2 with previous SOI list
```

```
##   |                                                                                                                                   |                                                                                                                           |   0%  |                                                                                                                                   |===============================                                                                                            |  25%  |                                                                                                                                   |==============================================================                                                             |  50%  |                                                                                                                                   |============================================================================================                               |  75%  |                                                                                                                                   |===========================================================================================================================| 100%
```

```
## Adding entries unique from the first DF
```

```
## Adding entries unique from the second DF
```

```
## Solving inner conflicts, Round: 1
```

```
## Solving inner conflicts, Round: 2
```

```
## Solving inner conflicts, Round: 3
```

```
## Initial peak retrieval:
```

```
##   |                                                                                                                                   |                                                                                                                           |   0%  |                                                                                                                                   |===============================                                                                                            |  25%  |                                                                                                                                   |==============================================================                                                             |  50%  |                                                                                                                                   |============================================================================================                               |  75%  |                                                                                                                                   |===========================================================================================================================| 100%
```

```
## 
```

```
## Starting group characterization:
```

```
## Mass calculation:
```

```
## Shortening and selecting long groups:
```

```
##   |                                                                                                                                   |                                                                                                                           |   0%  |                                                                                                                                   |===============================                                                                                            |  25%  |                                                                                                                                   |==============================================================                                                             |  50%  |                                                                                                                                   |============================================================================================                               |  75%  |                                                                                                                                   |===========================================================================================================================| 100%
```

```
## 
```

```
## Width calculation:
```

```
## Maximum intensity calculation:
```

```
## Number of scans:
```

```
## Converting from ionic formula to F/A combinations:
```

```
## Generating peaklist for plotting:
```

```
##   |                                                                                                                                   |                                                                                                                           |   0%  |                                                                                                                                   |===============================                                                                                            |  25%  |                                                                                                                                   |==============================================================                                                             |  50%  |                                                                                                                                   |============================================================================================                               |  75%  |                                                                                                                                   |===========================================================================================================================| 100%
```

```
## Starting density filtering:
```

```
## Running Density Filter
```

```
##   |                                                                                                                                   |                                                                                                                           |   0%  |                                                                                                                                   |===============================                                                                                            |  25%  |                                                                                                                                   |==============================================================                                                             |  50%  |                                                                                                                                   |============================================================================================                               |  75%  |                                                                                                                                   |===========================================================================================================================| 100%
```

```
## Running Density Interpreter
```

```
##   |                                                                                                                                   |                                                                                                                           |   0%  |                                                                                                                                   |===============================                                                                                            |  25%  |                                                                                                                                   |==============================================================                                                             |  50%  |                                                                                                                                   |============================================================================================                               |  75%  |                                                                                                                                   |===========================================================================================================================| 100%
```

```
## Running Density Filter
```

```
##   |                                                                                                                                   |                                                                                                                           |   0%  |                                                                                                                                   |===============================                                                                                            |  25%  |                                                                                                                                   |==============================================================                                                             |  50%  |                                                                                                                                   |============================================================================================                               |  75%  |                                                                                                                                   |===========================================================================================================================| 100%
```

```
## Running Density Interpreter
```

```
##   |                                                                                                                                   |                                                                                                                           |   0%  |                                                                                                                                   |===============================                                                                                            |  25%  |                                                                                                                                   |==============================================================                                                             |  50%  |                                                                                                                                   |============================================================================================                               |  75%  |                                                                                                                                   |===========================================================================================================================| 100%
```

```
## Now Grouping:
```

```
## Merging filter 2 out of 2 with previous SOI list
```

```
##   |                                                                                                                                   |                                                                                                                           |   0%  |                                                                                                                                   |===============================                                                                                            |  25%  |                                                                                                                                   |==============================================================                                                             |  50%  |                                                                                                                                   |============================================================================================                               |  75%  |                                                                                                                                   |===========================================================================================================================| 100%
```

```
## Adding entries unique from the first DF
```

```
## Adding entries unique from the second DF
```

```
## Solving inner conflicts, Round: 1
```

```
## Solving inner conflicts, Round: 2
```

```
## Solving inner conflicts, Round: 3
```

```
## Initial peak retrieval:
```

```
##   |                                                                                                                                   |                                                                                                                           |   0%  |                                                                                                                                   |===============================                                                                                            |  25%  |                                                                                                                                   |==============================================================                                                             |  50%  |                                                                                                                                   |============================================================================================                               |  75%  |                                                                                                                                   |===========================================================================================================================| 100%
```

```
## 
```

```
## Starting group characterization:
```

```
## Mass calculation:
```

```
## Shortening and selecting long groups:
```

```
##   |                                                                                                                                   |                                                                                                                           |   0%  |                                                                                                                                   |===============================                                                                                            |  25%  |                                                                                                                                   |==============================================================                                                             |  50%  |                                                                                                                                   |============================================================================================                               |  75%  |                                                                                                                                   |===========================================================================================================================| 100%
```

```
## Blank substraction:
```

```
##   |                                                                                                                                   |                                                                                                                           |   0%  |                                                                                                                                   |===============================                                                                                            |  25%  |                                                                                                                                   |==============================================================                                                             |  50%  |                                                                                                                                   |============================================================================================                               |  75%  |                                                                                                                                   |===========================================================================================================================| 100%
## 
##   |                                                                                                                                   |                                                                                                                           |   0%  |                                                                                                                                   |===============================                                                                                            |  25%  |                                                                                                                                   |==============================================================                                                             |  50%  |                                                                                                                                   |============================================================================================                               |  75%  |                                                                                                                                   |===========================================================================================================================| 100%
```

```
## 
```

```
## Width calculation:
```

```
## Maximum intensity calculation:
```

```
## Number of scans:
```

```
## Converting from ionic formula to F/A combinations:
```

```
## Generating peaklist for plotting:
```

```
##   |                                                                                                                                   |                                                                                                                           |   0%  |                                                                                                                                   |===============================                                                                                            |  25%  |                                                                                                                                   |==============================================================                                                             |  50%  |                                                                                                                                   |============================================================================================                               |  75%  |                                                                                                                                   |===========================================================================================================================| 100%
```

```
## Starting density filtering:
```

```
## Running Density Filter
```

```
##   |                                                                                                                                   |                                                                                                                           |   0%  |                                                                                                                                   |===============================                                                                                            |  25%  |                                                                                                                                   |==============================================================                                                             |  50%  |                                                                                                                                   |============================================================================================                               |  75%  |                                                                                                                                   |===========================================================================================================================| 100%
```

```
## Running Density Interpreter
```

```
##   |                                                                                                                                   |                                                                                                                           |   0%  |                                                                                                                                   |===============================                                                                                            |  25%  |                                                                                                                                   |==============================================================                                                             |  50%  |                                                                                                                                   |============================================================================================                               |  75%  |                                                                                                                                   |===========================================================================================================================| 100%
```

```
## Running Density Filter
```

```
##   |                                                                                                                                   |                                                                                                                           |   0%  |                                                                                                                                   |===============================                                                                                            |  25%  |                                                                                                                                   |==============================================================                                                             |  50%  |                                                                                                                                   |============================================================================================                               |  75%  |                                                                                                                                   |===========================================================================================================================| 100%
```

```
## Running Density Interpreter
```

```
##   |                                                                                                                                   |                                                                                                                           |   0%  |                                                                                                                                   |===============================                                                                                            |  25%  |                                                                                                                                   |==============================================================                                                             |  50%  |                                                                                                                                   |============================================================================================                               |  75%  |                                                                                                                                   |===========================================================================================================================| 100%
```

```
## Now Grouping:
```

```
## Merging filter 2 out of 2 with previous SOI list
```

```
##   |                                                                                                                                   |                                                                                                                           |   0%  |                                                                                                                                   |===============================                                                                                            |  25%  |                                                                                                                                   |==============================================================                                                             |  50%  |                                                                                                                                   |============================================================================================                               |  75%  |                                                                                                                                   |===========================================================================================================================| 100%
```

```
## Adding entries unique from the first DF
```

```
## Adding entries unique from the second DF
```

```
## Solving inner conflicts, Round: 1
```

```
## Solving inner conflicts, Round: 2
```

```
## Solving inner conflicts, Round: 3
```

```
## Initial peak retrieval:
```

```
##   |                                                                                                                                   |                                                                                                                           |   0%  |                                                                                                                                   |===============================                                                                                            |  25%  |                                                                                                                                   |==============================================================                                                             |  50%  |                                                                                                                                   |============================================================================================                               |  75%  |                                                                                                                                   |===========================================================================================================================| 100%
```

```
## 
```

```
## Starting group characterization:
```

```
## Mass calculation:
```

```
## Shortening and selecting long groups:
```

```
##   |                                                                                                                                   |                                                                                                                           |   0%  |                                                                                                                                   |===============================                                                                                            |  25%  |                                                                                                                                   |==============================================================                                                             |  50%  |                                                                                                                                   |============================================================================================                               |  75%  |                                                                                                                                   |===========================================================================================================================| 100%
```

```
## Blank substraction:
```

```
##   |                                                                                                                                   |                                                                                                                           |   0%  |                                                                                                                                   |===============================                                                                                            |  25%  |                                                                                                                                   |==============================================================                                                             |  50%  |                                                                                                                                   |============================================================================================                               |  75%  |                                                                                                                                   |===========================================================================================================================| 100%
## 
##   |                                                                                                                                   |                                                                                                                           |   0%  |                                                                                                                                   |===============================                                                                                            |  25%  |                                                                                                                                   |==============================================================                                                             |  50%  |                                                                                                                                   |============================================================================================                               |  75%  |                                                                                                                                   |===========================================================================================================================| 100%
```

```
## 
```

```
## Width calculation:
```

```
## Maximum intensity calculation:
```

```
## Number of scans:
```

```
## Converting from ionic formula to F/A combinations:
```

```
## Generating peaklist for plotting:
```

```
##   |                                                                                                                                   |                                                                                                                           |   0%  |                                                                                                                                   |===============================                                                                                            |  25%  |                                                                                                                                   |==============================================================                                                             |  50%  |                                                                                                                                   |============================================================================================                               |  75%  |                                                                                                                                   |===========================================================================================================================| 100%
```

```r
saveRDS(myHermes, file = "D:/HermesResults/SurfaceWater/PosSOI.rds")

for(i in 1:6){
  myHermes <- SOIcleaner(myHermes, soiid = i, minint = 10000, TRUE)
  myHermes <- ISFproc(myHermes, id = i)
}
```

```
## Computing isotopic elution similarity:
```

```
## Calculating isotope similarity:
```

```
##   |                                                                                                                                   |                                                                                                                           |   0%  |                                                                                                                                   |===============================                                                                                            |  25%  |                                                                                                                                   |==============================================================                                                             |  50%  |                                                                                                                                   |============================================================================================                               |  75%  |                                                                                                                                   |===========================================================================================================================| 100%
```

```
## Calculating isotopic fidelity metrics:
```

```
##   |                                                                                                                                   |                                                                                                                           |   0%  |                                                                                                                                   |=                                                                                                                          |   0%  |                                                                                                                                   |=                                                                                                                          |   1%  |                                                                                                                                   |==                                                                                                                         |   1%  |                                                                                                                                   |==                                                                                                                         |   2%  |                                                                                                                                   |===                                                                                                                        |   2%  |                                                                                                                                   |===                                                                                                                        |   3%  |                                                                                                                                   |====                                                                                                                       |   3%  |                                                                                                                                   |====                                                                                                                       |   4%  |                                                                                                                                   |=====                                                                                                                      |   4%  |                                                                                                                                   |======                                                                                                                     |   4%  |                                                                                                                                   |======                                                                                                                     |   5%  |                                                                                                                                   |=======                                                                                                                    |   5%  |                                                                                                                                   |=======                                                                                                                    |   6%  |                                                                                                                                   |========                                                                                                                   |   6%  |                                                                                                                                   |========                                                                                                                   |   7%  |                                                                                                                                   |=========                                                                                                                  |   7%  |                                                                                                                                   |=========                                                                                                                  |   8%  |                                                                                                                                   |==========                                                                                                                 |   8%  |                                                                                                                                   |==========                                                                                                                 |   9%  |                                                                                                                                   |===========                                                                                                                |   9%  |                                                                                                                                   |============                                                                                                               |   9%  |                                                                                                                                   |============                                                                                                               |  10%  |                                                                                                                                   |=============                                                                                                              |  10%  |                                                                                                                                   |=============                                                                                                              |  11%  |                                                                                                                                   |==============                                                                                                             |  11%  |                                                                                                                                   |==============                                                                                                             |  12%  |                                                                                                                                   |===============                                                                                                            |  12%  |                                                                                                                                   |===============                                                                                                            |  13%  |                                                                                                                                   |================                                                                                                           |  13%  |                                                                                                                                   |=================                                                                                                          |  13%  |                                                                                                                                   |=================                                                                                                          |  14%  |                                                                                                                                   |==================                                                                                                         |  14%  |                                                                                                                                   |==================                                                                                                         |  15%  |                                                                                                                                   |===================                                                                                                        |  15%  |                                                                                                                                   |===================                                                                                                        |  16%  |                                                                                                                                   |====================                                                                                                       |  16%  |                                                                                                                                   |====================                                                                                                       |  17%  |                                                                                                                                   |=====================                                                                                                      |  17%  |                                                                                                                                   |======================                                                                                                     |  17%  |                                                                                                                                   |======================                                                                                                     |  18%  |                                                                                                                                   |=======================                                                                                                    |  18%  |                                                                                                                                   |=======================                                                                                                    |  19%  |                                                                                                                                   |========================                                                                                                   |  19%  |                                                                                                                                   |========================                                                                                                   |  20%  |                                                                                                                                   |=========================                                                                                                  |  20%  |                                                                                                                                   |=========================                                                                                                  |  21%  |                                                                                                                                   |==========================                                                                                                 |  21%  |                                                                                                                                   |==========================                                                                                                 |  22%  |                                                                                                                                   |===========================                                                                                                |  22%  |                                                                                                                                   |============================                                                                                               |  22%  |                                                                                                                                   |============================                                                                                               |  23%  |                                                                                                                                   |=============================                                                                                              |  23%  |                                                                                                                                   |=============================                                                                                              |  24%  |                                                                                                                                   |==============================                                                                                             |  24%  |                                                                                                                                   |==============================                                                                                             |  25%  |                                                                                                                                   |===============================                                                                                            |  25%  |                                                                                                                                   |===============================                                                                                            |  26%  |                                                                                                                                   |================================                                                                                           |  26%  |                                                                                                                                   |=================================                                                                                          |  26%  |                                                                                                                                   |=================================                                                                                          |  27%  |                                                                                                                                   |==================================                                                                                         |  27%  |                                                                                                                                   |==================================                                                                                         |  28%  |                                                                                                                                   |===================================                                                                                        |  28%  |                                                                                                                                   |===================================                                                                                        |  29%  |                                                                                                                                   |====================================                                                                                       |  29%  |                                                                                                                                   |====================================                                                                                       |  30%  |                                                                                                                                   |=====================================                                                                                      |  30%  |                                                                                                                                   |======================================                                                                                     |  30%  |                                                                                                                                   |======================================                                                                                     |  31%  |                                                                                                                                   |=======================================                                                                                    |  31%  |                                                                                                                                   |=======================================                                                                                    |  32%  |                                                                                                                                   |========================================                                                                                   |  32%  |                                                                                                                                   |========================================                                                                                   |  33%  |                                                                                                                                   |=========================================                                                                                  |  33%  |                                                                                                                                   |=========================================                                                                                  |  34%  |                                                                                                                                   |==========================================                                                                                 |  34%  |                                                                                                                                   |==========================================                                                                                 |  35%  |                                                                                                                                   |===========================================                                                                                |  35%  |                                                                                                                                   |============================================                                                                               |  35%  |                                                                                                                                   |============================================                                                                               |  36%  |                                                                                                                                   |=============================================                                                                              |  36%  |                                                                                                                                   |=============================================                                                                              |  37%  |                                                                                                                                   |==============================================                                                                             |  37%  |                                                                                                                                   |==============================================                                                                             |  38%  |                                                                                                                                   |===============================================                                                                            |  38%  |                                                                                                                                   |===============================================                                                                            |  39%  |                                                                                                                                   |================================================                                                                           |  39%  |                                                                                                                                   |=================================================                                                                          |  39%  |                                                                                                                                   |=================================================                                                                          |  40%  |                                                                                                                                   |==================================================                                                                         |  40%  |                                                                                                                                   |==================================================                                                                         |  41%  |                                                                                                                                   |===================================================                                                                        |  41%  |                                                                                                                                   |===================================================                                                                        |  42%  |                                                                                                                                   |====================================================                                                                       |  42%  |                                                                                                                                   |====================================================                                                                       |  43%  |                                                                                                                                   |=====================================================                                                                      |  43%  |                                                                                                                                   |======================================================                                                                     |  44%  |                                                                                                                                   |=======================================================                                                                    |  44%  |                                                                                                                                   |=======================================================                                                                    |  45%  |                                                                                                                                   |========================================================                                                                   |  45%  |                                                                                                                                   |========================================================                                                                   |  46%  |                                                                                                                                   |=========================================================                                                                  |  46%  |                                                                                                                                   |=========================================================                                                                  |  47%  |                                                                                                                                   |==========================================================                                                                 |  47%  |                                                                                                                                   |==========================================================                                                                 |  48%  |                                                                                                                                   |===========================================================                                                                |  48%  |                                                                                                                                   |============================================================                                                               |  48%  |                                                                                                                                   |============================================================                                                               |  49%  |                                                                                                                                   |=============================================================                                                              |  49%  |                                                                                                                                   |=============================================================                                                              |  50%  |                                                                                                                                   |==============================================================                                                             |  50%  |                                                                                                                                   |==============================================================                                                             |  51%  |                                                                                                                                   |===============================================================                                                            |  51%  |                                                                                                                                   |===============================================================                                                            |  52%  |                                                                                                                                   |================================================================                                                           |  52%  |                                                                                                                                   |=================================================================                                                          |  52%  |                                                                                                                                   |=================================================================                                                          |  53%  |                                                                                                                                   |==================================================================                                                         |  53%  |                                                                                                                                   |==================================================================                                                         |  54%  |                                                                                                                                   |===================================================================                                                        |  54%  |                                                                                                                                   |===================================================================                                                        |  55%  |                                                                                                                                   |====================================================================                                                       |  55%  |                                                                                                                                   |====================================================================                                                       |  56%  |                                                                                                                                   |=====================================================================                                                      |  56%  |                                                                                                                                   |======================================================================                                                     |  57%  |                                                                                                                                   |=======================================================================                                                    |  57%  |                                                                                                                                   |=======================================================================                                                    |  58%  |                                                                                                                                   |========================================================================                                                   |  58%  |                                                                                                                                   |========================================================================                                                   |  59%  |                                                                                                                                   |=========================================================================                                                  |  59%  |                                                                                                                                   |=========================================================================                                                  |  60%  |                                                                                                                                   |==========================================================================                                                 |  60%  |                                                                                                                                   |==========================================================================                                                 |  61%  |                                                                                                                                   |===========================================================================                                                |  61%  |                                                                                                                                   |============================================================================                                               |  61%  |                                                                                                                                   |============================================================================                                               |  62%  |                                                                                                                                   |=============================================================================                                              |  62%  |                                                                                                                                   |=============================================================================                                              |  63%  |                                                                                                                                   |==============================================================================                                             |  63%  |                                                                                                                                   |==============================================================================                                             |  64%  |                                                                                                                                   |===============================================================================                                            |  64%  |                                                                                                                                   |===============================================================================                                            |  65%  |                                                                                                                                   |================================================================================                                           |  65%  |                                                                                                                                   |=================================================================================                                          |  65%  |                                                                                                                                   |=================================================================================                                          |  66%  |                                                                                                                                   |==================================================================================                                         |  66%  |                                                                                                                                   |==================================================================================                                         |  67%  |                                                                                                                                   |===================================================================================                                        |  67%  |                                                                                                                                   |===================================================================================                                        |  68%  |                                                                                                                                   |====================================================================================                                       |  68%  |                                                                                                                                   |====================================================================================                                       |  69%  |                                                                                                                                   |=====================================================================================                                      |  69%  |                                                                                                                                   |=====================================================================================                                      |  70%  |                                                                                                                                   |======================================================================================                                     |  70%  |                                                                                                                                   |=======================================================================================                                    |  70%  |                                                                                                                                   |=======================================================================================                                    |  71%  |                                                                                                                                   |========================================================================================                                   |  71%  |                                                                                                                                   |========================================================================================                                   |  72%  |                                                                                                                                   |=========================================================================================                                  |  72%  |                                                                                                                                   |=========================================================================================                                  |  73%  |                                                                                                                                   |==========================================================================================                                 |  73%  |                                                                                                                                   |==========================================================================================                                 |  74%  |                                                                                                                                   |===========================================================================================                                |  74%  |                                                                                                                                   |============================================================================================                               |  74%  |                                                                                                                                   |============================================================================================                               |  75%  |                                                                                                                                   |=============================================================================================                              |  75%  |                                                                                                                                   |=============================================================================================                              |  76%  |                                                                                                                                   |==============================================================================================                             |  76%  |                                                                                                                                   |==============================================================================================                             |  77%  |                                                                                                                                   |===============================================================================================                            |  77%  |                                                                                                                                   |===============================================================================================                            |  78%  |                                                                                                                                   |================================================================================================                           |  78%  |                                                                                                                                   |=================================================================================================                          |  78%  |                                                                                                                                   |=================================================================================================                          |  79%  |                                                                                                                                   |==================================================================================================                         |  79%  |                                                                                                                                   |==================================================================================================                         |  80%  |                                                                                                                                   |===================================================================================================                        |  80%  |                                                                                                                                   |===================================================================================================                        |  81%  |                                                                                                                                   |====================================================================================================                       |  81%  |                                                                                                                                   |====================================================================================================                       |  82%  |                                                                                                                                   |=====================================================================================================                      |  82%  |                                                                                                                                   |=====================================================================================================                      |  83%  |                                                                                                                                   |======================================================================================================                     |  83%  |                                                                                                                                   |=======================================================================================================                    |  83%  |                                                                                                                                   |=======================================================================================================                    |  84%  |                                                                                                                                   |========================================================================================================                   |  84%  |                                                                                                                                   |========================================================================================================                   |  85%  |                                                                                                                                   |=========================================================================================================                  |  85%  |                                                                                                                                   |=========================================================================================================                  |  86%  |                                                                                                                                   |==========================================================================================================                 |  86%  |                                                                                                                                   |==========================================================================================================                 |  87%  |                                                                                                                                   |===========================================================================================================                |  87%  |                                                                                                                                   |============================================================================================================               |  87%  |                                                                                                                                   |============================================================================================================               |  88%  |                                                                                                                                   |=============================================================================================================              |  88%  |                                                                                                                                   |=============================================================================================================              |  89%  |                                                                                                                                   |==============================================================================================================             |  89%  |                                                                                                                                   |==============================================================================================================             |  90%  |                                                                                                                                   |===============================================================================================================            |  90%  |                                                                                                                                   |===============================================================================================================            |  91%  |                                                                                                                                   |================================================================================================================           |  91%  |                                                                                                                                   |=================================================================================================================          |  91%  |                                                                                                                                   |=================================================================================================================          |  92%  |                                                                                                                                   |==================================================================================================================         |  92%  |                                                                                                                                   |==================================================================================================================         |  93%  |                                                                                                                                   |===================================================================================================================        |  93%  |                                                                                                                                   |===================================================================================================================        |  94%  |                                                                                                                                   |====================================================================================================================       |  94%  |                                                                                                                                   |====================================================================================================================       |  95%  |                                                                                                                                   |=====================================================================================================================      |  95%  |                                                                                                                                   |=====================================================================================================================      |  96%  |                                                                                                                                   |======================================================================================================================     |  96%  |                                                                                                                                   |=======================================================================================================================    |  96%  |                                                                                                                                   |=======================================================================================================================    |  97%  |                                                                                                                                   |========================================================================================================================   |  97%  |                                                                                                                                   |========================================================================================================================   |  98%  |                                                                                                                                   |=========================================================================================================================  |  98%  |                                                                                                                                   |=========================================================================================================================  |  99%  |                                                                                                                                   |========================================================================================================================== |  99%  |                                                                                                                                   |========================================================================================================================== | 100%  |                                                                                                                                   |===========================================================================================================================| 100%
```

```
## Removing confirmed isotope entries:
```

```
## Recalculating peaklist for plotting:
```

```
##   |                                                                                                                                   |                                                                                                                           |   0%  |                                                                                                                                   |===============================                                                                                            |  25%  |                                                                                                                                   |==============================================================                                                             |  50%  |                                                                                                                                   |============================================================================================                               |  75%  |                                                                                                                                   |===========================================================================================================================| 100%
```

```
## Calculating adduct similarity:
```

```
##   |                                                                                                                                   |                                                                                                                           |   0%  |                                                                                                                                   |===============================                                                                                            |  25%  |                                                                                                                                   |==============================================================                                                             |  50%  |                                                                                                                                   |============================================================================================                               |  75%  |                                                                                                                                   |===========================================================================================================================| 100%
```

```
## 100 out of 14162
```

```
## 200 out of 14162
```

```
## 300 out of 14162
```

```
## 400 out of 14162
```

```
## 500 out of 14162
```

```
## 600 out of 14162
```

```
## 700 out of 14162
```

```
## 800 out of 14162
```

```
## 900 out of 14162
```

```
## 1000 out of 14162
```

```
## 1100 out of 14162
```

```
## 1200 out of 14162
```

```
## 1300 out of 14162
```

```
## 1400 out of 14162
```

```
## 1500 out of 14162
```

```
## 1600 out of 14162
```

```
## 1700 out of 14162
```

```
## 1800 out of 14162
```

```
## 1900 out of 14162
```

```
## 2000 out of 14162
```

```
## 2100 out of 14162
```

```
## 2200 out of 14162
```

```
## 2300 out of 14162
```

```
## 2400 out of 14162
```

```
## 2500 out of 14162
```

```
## 2600 out of 14162
```

```
## 2700 out of 14162
```

```
## 2800 out of 14162
```

```
## 2900 out of 14162
```

```
## 3000 out of 14162
```

```
## 3100 out of 14162
```

```
## 3200 out of 14162
```

```
## 3300 out of 14162
```

```
## 3400 out of 14162
```

```
## 3500 out of 14162
```

```
## 3600 out of 14162
```

```
## 3700 out of 14162
```

```
## 3800 out of 14162
```

```
## 3900 out of 14162
```

```
## 4000 out of 14162
```

```
## 4100 out of 14162
```

```
## 4200 out of 14162
```

```
## 4300 out of 14162
```

```
## 4400 out of 14162
```

```
## 4500 out of 14162
```

```
## 4600 out of 14162
```

```
## 4700 out of 14162
```

```
## 4800 out of 14162
```

```
## 4900 out of 14162
```

```
## 5000 out of 14162
```

```
## 5100 out of 14162
```

```
## 5200 out of 14162
```

```
## 5300 out of 14162
```

```
## 5400 out of 14162
```

```
## 5500 out of 14162
```

```
## 5600 out of 14162
```

```
## 5700 out of 14162
```

```
## 5800 out of 14162
```

```
## 5900 out of 14162
```

```
## 6000 out of 14162
```

```
## 6100 out of 14162
```

```
## 6200 out of 14162
```

```
## 6300 out of 14162
```

```
## 6400 out of 14162
```

```
## 6500 out of 14162
```

```
## 6600 out of 14162
```

```
## 6700 out of 14162
```

```
## 6800 out of 14162
```

```
## 6900 out of 14162
```

```
## 7000 out of 14162
```

```
## 7100 out of 14162
```

```
## 7200 out of 14162
```

```
## 7300 out of 14162
```

```
## 7400 out of 14162
```

```
## 7500 out of 14162
```

```
## 7600 out of 14162
```

```
## 7700 out of 14162
```

```
## 7800 out of 14162
```

```
## 7900 out of 14162
```

```
## 8000 out of 14162
```

```
## 8100 out of 14162
```

```
## 8200 out of 14162
```

```
## 8300 out of 14162
```

```
## 8400 out of 14162
```

```
## 8500 out of 14162
```

```
## 8600 out of 14162
```

```
## 8700 out of 14162
```

```
## 8800 out of 14162
```

```
## 8900 out of 14162
```

```
## 9000 out of 14162
```

```
## 9100 out of 14162
```

```
## 9200 out of 14162
```

```
## 9300 out of 14162
```

```
## 9400 out of 14162
```

```
## 9500 out of 14162
```

```
## 9600 out of 14162
```

```
## 9700 out of 14162
```

```
## 9800 out of 14162
```

```
## 9900 out of 14162
```

```
## 10000 out of 14162
```

```
## 10100 out of 14162
```

```
## 10200 out of 14162
```

```
## 10300 out of 14162
```

```
## 10400 out of 14162
```

```
## 10500 out of 14162
```

```
## 10600 out of 14162
```

```
## 10700 out of 14162
```

```
## 10800 out of 14162
```

```
## 10900 out of 14162
```

```
## 11000 out of 14162
```

```
## 11100 out of 14162
```

```
## 11200 out of 14162
```

```
## 11300 out of 14162
```

```
## 11400 out of 14162
```

```
## 11500 out of 14162
```

```
## 11600 out of 14162
```

```
## 11700 out of 14162
```

```
## 11800 out of 14162
```

```
## 11900 out of 14162
```

```
## 12000 out of 14162
```

```
## 12100 out of 14162
```

```
## 12200 out of 14162
```

```
## 12300 out of 14162
```

```
## 12400 out of 14162
```

```
## 12500 out of 14162
```

```
## 12600 out of 14162
```

```
## 12700 out of 14162
```

```
## 12800 out of 14162
```

```
## 12900 out of 14162
```

```
## 13000 out of 14162
```

```
## 13100 out of 14162
```

```
## 13200 out of 14162
```

```
## 13300 out of 14162
```

```
## 13400 out of 14162
```

```
## 13500 out of 14162
```

```
## 13600 out of 14162
```

```
## 13700 out of 14162
```

```
## 13800 out of 14162
```

```
## 13900 out of 14162
```

```
## 14000 out of 14162
```

```
## 14100 out of 14162
```

```
##   |                                                                                                                                   |                                                                                                                           |   0%  |                                                                                                                                   |===============================                                                                                            |  25%  |                                                                                                                                   |==============================================================                                                             |  50%  |                                                                                                                                   |============================================================================================                               |  75%  |                                                                                                                                   |===========================================================================================================================| 100%
```

```
## Computing isotopic elution similarity:
```

```
## Calculating isotope similarity:
```

```
##   |                                                                                                                                   |                                                                                                                           |   0%  |                                                                                                                                   |===============================                                                                                            |  25%  |                                                                                                                                   |==============================================================                                                             |  50%  |                                                                                                                                   |============================================================================================                               |  75%  |                                                                                                                                   |===========================================================================================================================| 100%
```

```
## Calculating isotopic fidelity metrics:
```

```
##   |                                                                                                                                   |                                                                                                                           |   0%  |                                                                                                                                   |=                                                                                                                          |   0%  |                                                                                                                                   |=                                                                                                                          |   1%  |                                                                                                                                   |==                                                                                                                         |   1%  |                                                                                                                                   |==                                                                                                                         |   2%  |                                                                                                                                   |===                                                                                                                        |   2%  |                                                                                                                                   |===                                                                                                                        |   3%  |                                                                                                                                   |====                                                                                                                       |   3%  |                                                                                                                                   |====                                                                                                                       |   4%  |                                                                                                                                   |=====                                                                                                                      |   4%  |                                                                                                                                   |======                                                                                                                     |   4%  |                                                                                                                                   |======                                                                                                                     |   5%  |                                                                                                                                   |=======                                                                                                                    |   5%  |                                                                                                                                   |=======                                                                                                                    |   6%  |                                                                                                                                   |========                                                                                                                   |   6%  |                                                                                                                                   |========                                                                                                                   |   7%  |                                                                                                                                   |=========                                                                                                                  |   7%  |                                                                                                                                   |=========                                                                                                                  |   8%  |                                                                                                                                   |==========                                                                                                                 |   8%  |                                                                                                                                   |==========                                                                                                                 |   9%  |                                                                                                                                   |===========                                                                                                                |   9%  |                                                                                                                                   |============                                                                                                               |   9%  |                                                                                                                                   |============                                                                                                               |  10%  |                                                                                                                                   |=============                                                                                                              |  10%  |                                                                                                                                   |=============                                                                                                              |  11%  |                                                                                                                                   |==============                                                                                                             |  11%  |                                                                                                                                   |==============                                                                                                             |  12%  |                                                                                                                                   |===============                                                                                                            |  12%  |                                                                                                                                   |===============                                                                                                            |  13%  |                                                                                                                                   |================                                                                                                           |  13%  |                                                                                                                                   |=================                                                                                                          |  13%  |                                                                                                                                   |=================                                                                                                          |  14%  |                                                                                                                                   |==================                                                                                                         |  14%  |                                                                                                                                   |==================                                                                                                         |  15%  |                                                                                                                                   |===================                                                                                                        |  15%  |                                                                                                                                   |===================                                                                                                        |  16%  |                                                                                                                                   |====================                                                                                                       |  16%  |                                                                                                                                   |====================                                                                                                       |  17%  |                                                                                                                                   |=====================                                                                                                      |  17%  |                                                                                                                                   |======================                                                                                                     |  17%  |                                                                                                                                   |======================                                                                                                     |  18%  |                                                                                                                                   |=======================                                                                                                    |  18%  |                                                                                                                                   |=======================                                                                                                    |  19%  |                                                                                                                                   |========================                                                                                                   |  19%  |                                                                                                                                   |========================                                                                                                   |  20%  |                                                                                                                                   |=========================                                                                                                  |  20%  |                                                                                                                                   |=========================                                                                                                  |  21%  |                                                                                                                                   |==========================                                                                                                 |  21%  |                                                                                                                                   |==========================                                                                                                 |  22%  |                                                                                                                                   |===========================                                                                                                |  22%  |                                                                                                                                   |============================                                                                                               |  22%  |                                                                                                                                   |============================                                                                                               |  23%  |                                                                                                                                   |=============================                                                                                              |  23%  |                                                                                                                                   |=============================                                                                                              |  24%  |                                                                                                                                   |==============================                                                                                             |  24%  |                                                                                                                                   |==============================                                                                                             |  25%  |                                                                                                                                   |===============================                                                                                            |  25%  |                                                                                                                                   |===============================                                                                                            |  26%  |                                                                                                                                   |================================                                                                                           |  26%  |                                                                                                                                   |=================================                                                                                          |  26%  |                                                                                                                                   |=================================                                                                                          |  27%  |                                                                                                                                   |==================================                                                                                         |  27%  |                                                                                                                                   |==================================                                                                                         |  28%  |                                                                                                                                   |===================================                                                                                        |  28%  |                                                                                                                                   |===================================                                                                                        |  29%  |                                                                                                                                   |====================================                                                                                       |  29%  |                                                                                                                                   |====================================                                                                                       |  30%  |                                                                                                                                   |=====================================                                                                                      |  30%  |                                                                                                                                   |======================================                                                                                     |  30%  |                                                                                                                                   |======================================                                                                                     |  31%  |                                                                                                                                   |=======================================                                                                                    |  31%  |                                                                                                                                   |=======================================                                                                                    |  32%  |                                                                                                                                   |========================================                                                                                   |  32%  |                                                                                                                                   |========================================                                                                                   |  33%  |                                                                                                                                   |=========================================                                                                                  |  33%  |                                                                                                                                   |=========================================                                                                                  |  34%  |                                                                                                                                   |==========================================                                                                                 |  34%  |                                                                                                                                   |==========================================                                                                                 |  35%  |                                                                                                                                   |===========================================                                                                                |  35%  |                                                                                                                                   |============================================                                                                               |  35%  |                                                                                                                                   |============================================                                                                               |  36%  |                                                                                                                                   |=============================================                                                                              |  36%  |                                                                                                                                   |=============================================                                                                              |  37%  |                                                                                                                                   |==============================================                                                                             |  37%  |                                                                                                                                   |==============================================                                                                             |  38%  |                                                                                                                                   |===============================================                                                                            |  38%  |                                                                                                                                   |===============================================                                                                            |  39%  |                                                                                                                                   |================================================                                                                           |  39%  |                                                                                                                                   |=================================================                                                                          |  39%  |                                                                                                                                   |=================================================                                                                          |  40%  |                                                                                                                                   |==================================================                                                                         |  40%  |                                                                                                                                   |==================================================                                                                         |  41%  |                                                                                                                                   |===================================================                                                                        |  41%  |                                                                                                                                   |===================================================                                                                        |  42%  |                                                                                                                                   |====================================================                                                                       |  42%  |                                                                                                                                   |====================================================                                                                       |  43%  |                                                                                                                                   |=====================================================                                                                      |  43%  |                                                                                                                                   |======================================================                                                                     |  43%  |                                                                                                                                   |======================================================                                                                     |  44%  |                                                                                                                                   |=======================================================                                                                    |  44%  |                                                                                                                                   |=======================================================                                                                    |  45%  |                                                                                                                                   |========================================================                                                                   |  45%  |                                                                                                                                   |========================================================                                                                   |  46%  |                                                                                                                                   |=========================================================                                                                  |  46%  |                                                                                                                                   |=========================================================                                                                  |  47%  |                                                                                                                                   |==========================================================                                                                 |  47%  |                                                                                                                                   |==========================================================                                                                 |  48%  |                                                                                                                                   |===========================================================                                                                |  48%  |                                                                                                                                   |============================================================                                                               |  48%  |                                                                                                                                   |============================================================                                                               |  49%  |                                                                                                                                   |=============================================================                                                              |  49%  |                                                                                                                                   |=============================================================                                                              |  50%  |                                                                                                                                   |==============================================================                                                             |  50%  |                                                                                                                                   |==============================================================                                                             |  51%  |                                                                                                                                   |===============================================================                                                            |  51%  |                                                                                                                                   |===============================================================                                                            |  52%  |                                                                                                                                   |================================================================                                                           |  52%  |                                                                                                                                   |=================================================================                                                          |  52%  |                                                                                                                                   |=================================================================                                                          |  53%  |                                                                                                                                   |==================================================================                                                         |  53%  |                                                                                                                                   |==================================================================                                                         |  54%  |                                                                                                                                   |===================================================================                                                        |  54%  |                                                                                                                                   |===================================================================                                                        |  55%  |                                                                                                                                   |====================================================================                                                       |  55%  |                                                                                                                                   |====================================================================                                                       |  56%  |                                                                                                                                   |=====================================================================                                                      |  56%  |                                                                                                                                   |=====================================================================                                                      |  57%  |                                                                                                                                   |======================================================================                                                     |  57%  |                                                                                                                                   |=======================================================================                                                    |  57%  |                                                                                                                                   |=======================================================================                                                    |  58%  |                                                                                                                                   |========================================================================                                                   |  58%  |                                                                                                                                   |========================================================================                                                   |  59%  |                                                                                                                                   |=========================================================================                                                  |  59%  |                                                                                                                                   |=========================================================================                                                  |  60%  |                                                                                                                                   |==========================================================================                                                 |  60%  |                                                                                                                                   |==========================================================================                                                 |  61%  |                                                                                                                                   |===========================================================================                                                |  61%  |                                                                                                                                   |============================================================================                                               |  61%  |                                                                                                                                   |============================================================================                                               |  62%  |                                                                                                                                   |=============================================================================                                              |  62%  |                                                                                                                                   |=============================================================================                                              |  63%  |                                                                                                                                   |==============================================================================                                             |  63%  |                                                                                                                                   |==============================================================================                                             |  64%  |                                                                                                                                   |===============================================================================                                            |  64%  |                                                                                                                                   |===============================================================================                                            |  65%  |                                                                                                                                   |================================================================================                                           |  65%  |                                                                                                                                   |=================================================================================                                          |  65%  |                                                                                                                                   |=================================================================================                                          |  66%  |                                                                                                                                   |==================================================================================                                         |  66%  |                                                                                                                                   |==================================================================================                                         |  67%  |                                                                                                                                   |===================================================================================                                        |  67%  |                                                                                                                                   |===================================================================================                                        |  68%  |                                                                                                                                   |====================================================================================                                       |  68%  |                                                                                                                                   |====================================================================================                                       |  69%  |                                                                                                                                   |=====================================================================================                                      |  69%  |                                                                                                                                   |=====================================================================================                                      |  70%  |                                                                                                                                   |======================================================================================                                     |  70%  |                                                                                                                                   |=======================================================================================                                    |  70%  |                                                                                                                                   |=======================================================================================                                    |  71%  |                                                                                                                                   |========================================================================================                                   |  71%  |                                                                                                                                   |========================================================================================                                   |  72%  |                                                                                                                                   |=========================================================================================                                  |  72%  |                                                                                                                                   |=========================================================================================                                  |  73%  |                                                                                                                                   |==========================================================================================                                 |  73%  |                                                                                                                                   |==========================================================================================                                 |  74%  |                                                                                                                                   |===========================================================================================                                |  74%  |                                                                                                                                   |============================================================================================                               |  74%  |                                                                                                                                   |============================================================================================                               |  75%  |                                                                                                                                   |=============================================================================================                              |  75%  |                                                                                                                                   |=============================================================================================                              |  76%  |                                                                                                                                   |==============================================================================================                             |  76%  |                                                                                                                                   |==============================================================================================                             |  77%  |                                                                                                                                   |===============================================================================================                            |  77%  |                                                                                                                                   |===============================================================================================                            |  78%  |                                                                                                                                   |================================================================================================                           |  78%  |                                                                                                                                   |=================================================================================================                          |  78%  |                                                                                                                                   |=================================================================================================                          |  79%  |                                                                                                                                   |==================================================================================================                         |  79%  |                                                                                                                                   |==================================================================================================                         |  80%  |                                                                                                                                   |===================================================================================================                        |  80%  |                                                                                                                                   |===================================================================================================                        |  81%  |                                                                                                                                   |====================================================================================================                       |  81%  |                                                                                                                                   |====================================================================================================                       |  82%  |                                                                                                                                   |=====================================================================================================                      |  82%  |                                                                                                                                   |=====================================================================================================                      |  83%  |                                                                                                                                   |======================================================================================================                     |  83%  |                                                                                                                                   |=======================================================================================================                    |  83%  |                                                                                                                                   |=======================================================================================================                    |  84%  |                                                                                                                                   |========================================================================================================                   |  84%  |                                                                                                                                   |========================================================================================================                   |  85%  |                                                                                                                                   |=========================================================================================================                  |  85%  |                                                                                                                                   |=========================================================================================================                  |  86%  |                                                                                                                                   |==========================================================================================================                 |  86%  |                                                                                                                                   |==========================================================================================================                 |  87%  |                                                                                                                                   |===========================================================================================================                |  87%  |                                                                                                                                   |============================================================================================================               |  87%  |                                                                                                                                   |============================================================================================================               |  88%  |                                                                                                                                   |=============================================================================================================              |  88%  |                                                                                                                                   |=============================================================================================================              |  89%  |                                                                                                                                   |==============================================================================================================             |  89%  |                                                                                                                                   |==============================================================================================================             |  90%  |                                                                                                                                   |===============================================================================================================            |  90%  |                                                                                                                                   |===============================================================================================================            |  91%  |                                                                                                                                   |================================================================================================================           |  91%  |                                                                                                                                   |=================================================================================================================          |  91%  |                                                                                                                                   |=================================================================================================================          |  92%  |                                                                                                                                   |==================================================================================================================         |  92%  |                                                                                                                                   |==================================================================================================================         |  93%  |                                                                                                                                   |===================================================================================================================        |  93%  |                                                                                                                                   |===================================================================================================================        |  94%  |                                                                                                                                   |====================================================================================================================       |  94%  |                                                                                                                                   |====================================================================================================================       |  95%  |                                                                                                                                   |=====================================================================================================================      |  95%  |                                                                                                                                   |=====================================================================================================================      |  96%  |                                                                                                                                   |======================================================================================================================     |  96%  |                                                                                                                                   |=======================================================================================================================    |  96%  |                                                                                                                                   |=======================================================================================================================    |  97%  |                                                                                                                                   |========================================================================================================================   |  97%  |                                                                                                                                   |========================================================================================================================   |  98%  |                                                                                                                                   |=========================================================================================================================  |  98%  |                                                                                                                                   |=========================================================================================================================  |  99%  |                                                                                                                                   |========================================================================================================================== |  99%  |                                                                                                                                   |========================================================================================================================== | 100%  |                                                                                                                                   |===========================================================================================================================| 100%
```

```
## Removing confirmed isotope entries:
```

```
## Recalculating peaklist for plotting:
```

```
##   |                                                                                                                                   |                                                                                                                           |   0%  |                                                                                                                                   |===============================                                                                                            |  25%  |                                                                                                                                   |==============================================================                                                             |  50%  |                                                                                                                                   |============================================================================================                               |  75%  |                                                                                                                                   |===========================================================================================================================| 100%
```

```
## Calculating adduct similarity:
```

```
##   |                                                                                                                                   |                                                                                                                           |   0%  |                                                                                                                                   |===============================                                                                                            |  25%  |                                                                                                                                   |==============================================================                                                             |  50%  |                                                                                                                                   |============================================================================================                               |  75%  |                                                                                                                                   |===========================================================================================================================| 100%
```

```
## 100 out of 14826
```

```
## 200 out of 14826
```

```
## 300 out of 14826
```

```
## 400 out of 14826
```

```
## 500 out of 14826
```

```
## 600 out of 14826
```

```
## 700 out of 14826
```

```
## 800 out of 14826
```

```
## 900 out of 14826
```

```
## 1000 out of 14826
```

```
## 1100 out of 14826
```

```
## 1200 out of 14826
```

```
## 1300 out of 14826
```

```
## 1400 out of 14826
```

```
## 1500 out of 14826
```

```
## 1600 out of 14826
```

```
## 1700 out of 14826
```

```
## 1800 out of 14826
```

```
## 1900 out of 14826
```

```
## 2000 out of 14826
```

```
## 2100 out of 14826
```

```
## 2200 out of 14826
```

```
## 2300 out of 14826
```

```
## 2400 out of 14826
```

```
## 2500 out of 14826
```

```
## 2600 out of 14826
```

```
## 2700 out of 14826
```

```
## 2800 out of 14826
```

```
## 2900 out of 14826
```

```
## 3000 out of 14826
```

```
## 3100 out of 14826
```

```
## 3200 out of 14826
```

```
## 3300 out of 14826
```

```
## 3400 out of 14826
```

```
## 3500 out of 14826
```

```
## 3600 out of 14826
```

```
## 3700 out of 14826
```

```
## 3800 out of 14826
```

```
## 3900 out of 14826
```

```
## 4000 out of 14826
```

```
## 4100 out of 14826
```

```
## 4200 out of 14826
```

```
## 4300 out of 14826
```

```
## 4400 out of 14826
```

```
## 4500 out of 14826
```

```
## 4600 out of 14826
```

```
## 4700 out of 14826
```

```
## 4800 out of 14826
```

```
## 4900 out of 14826
```

```
## 5000 out of 14826
```

```
## 5100 out of 14826
```

```
## 5200 out of 14826
```

```
## 5300 out of 14826
```

```
## 5400 out of 14826
```

```
## 5500 out of 14826
```

```
## 5600 out of 14826
```

```
## 5700 out of 14826
```

```
## 5800 out of 14826
```

```
## 5900 out of 14826
```

```
## 6000 out of 14826
```

```
## 6100 out of 14826
```

```
## 6200 out of 14826
```

```
## 6300 out of 14826
```

```
## 6400 out of 14826
```

```
## 6500 out of 14826
```

```
## 6600 out of 14826
```

```
## 6700 out of 14826
```

```
## 6800 out of 14826
```

```
## 6900 out of 14826
```

```
## 7000 out of 14826
```

```
## 7100 out of 14826
```

```
## 7200 out of 14826
```

```
## 7300 out of 14826
```

```
## 7400 out of 14826
```

```
## 7500 out of 14826
```

```
## 7600 out of 14826
```

```
## 7700 out of 14826
```

```
## 7800 out of 14826
```

```
## 7900 out of 14826
```

```
## 8000 out of 14826
```

```
## 8100 out of 14826
```

```
## 8200 out of 14826
```

```
## 8300 out of 14826
```

```
## 8400 out of 14826
```

```
## 8500 out of 14826
```

```
## 8600 out of 14826
```

```
## 8700 out of 14826
```

```
## 8800 out of 14826
```

```
## 8900 out of 14826
```

```
## 9000 out of 14826
```

```
## 9100 out of 14826
```

```
## 9200 out of 14826
```

```
## 9300 out of 14826
```

```
## 9400 out of 14826
```

```
## 9500 out of 14826
```

```
## 9600 out of 14826
```

```
## 9700 out of 14826
```

```
## 9800 out of 14826
```

```
## 9900 out of 14826
```

```
## 10000 out of 14826
```

```
## 10100 out of 14826
```

```
## 10200 out of 14826
```

```
## 10300 out of 14826
```

```
## 10400 out of 14826
```

```
## 10500 out of 14826
```

```
## 10600 out of 14826
```

```
## 10700 out of 14826
```

```
## 10800 out of 14826
```

```
## 10900 out of 14826
```

```
## 11000 out of 14826
```

```
## 11100 out of 14826
```

```
## 11200 out of 14826
```

```
## 11300 out of 14826
```

```
## 11400 out of 14826
```

```
## 11500 out of 14826
```

```
## 11600 out of 14826
```

```
## 11700 out of 14826
```

```
## 11800 out of 14826
```

```
## 11900 out of 14826
```

```
## 12000 out of 14826
```

```
## 12100 out of 14826
```

```
## 12200 out of 14826
```

```
## 12300 out of 14826
```

```
## 12400 out of 14826
```

```
## 12500 out of 14826
```

```
## 12600 out of 14826
```

```
## 12700 out of 14826
```

```
## 12800 out of 14826
```

```
## 12900 out of 14826
```

```
## 13000 out of 14826
```

```
## 13100 out of 14826
```

```
## 13200 out of 14826
```

```
## 13300 out of 14826
```

```
## 13400 out of 14826
```

```
## 13500 out of 14826
```

```
## 13600 out of 14826
```

```
## 13700 out of 14826
```

```
## 13800 out of 14826
```

```
## 13900 out of 14826
```

```
## 14000 out of 14826
```

```
## 14100 out of 14826
```

```
## 14200 out of 14826
```

```
## 14300 out of 14826
```

```
## 14400 out of 14826
```

```
## 14500 out of 14826
```

```
## 14600 out of 14826
```

```
## 14700 out of 14826
```

```
## 14800 out of 14826
```

```
##   |                                                                                                                                   |                                                                                                                           |   0%  |                                                                                                                                   |===============================                                                                                            |  25%  |                                                                                                                                   |==============================================================                                                             |  50%  |                                                                                                                                   |============================================================================================                               |  75%  |                                                                                                                                   |===========================================================================================================================| 100%
```

```
## Computing isotopic elution similarity:
```

```
## Calculating isotope similarity:
```

```
##   |                                                                                                                                   |                                                                                                                           |   0%  |                                                                                                                                   |===============================                                                                                            |  25%  |                                                                                                                                   |==============================================================                                                             |  50%  |                                                                                                                                   |============================================================================================                               |  75%  |                                                                                                                                   |===========================================================================================================================| 100%
```

```
## Calculating isotopic fidelity metrics:
```

```
##   |                                                                                                                                   |                                                                                                                           |   0%  |                                                                                                                                   |=                                                                                                                          |   0%  |                                                                                                                                   |=                                                                                                                          |   1%  |                                                                                                                                   |==                                                                                                                         |   1%  |                                                                                                                                   |==                                                                                                                         |   2%  |                                                                                                                                   |===                                                                                                                        |   2%  |                                                                                                                                   |===                                                                                                                        |   3%  |                                                                                                                                   |====                                                                                                                       |   3%  |                                                                                                                                   |====                                                                                                                       |   4%  |                                                                                                                                   |=====                                                                                                                      |   4%  |                                                                                                                                   |======                                                                                                                     |   4%  |                                                                                                                                   |======                                                                                                                     |   5%  |                                                                                                                                   |=======                                                                                                                    |   5%  |                                                                                                                                   |=======                                                                                                                    |   6%  |                                                                                                                                   |========                                                                                                                   |   6%  |                                                                                                                                   |========                                                                                                                   |   7%  |                                                                                                                                   |=========                                                                                                                  |   7%  |                                                                                                                                   |=========                                                                                                                  |   8%  |                                                                                                                                   |==========                                                                                                                 |   8%  |                                                                                                                                   |==========                                                                                                                 |   9%  |                                                                                                                                   |===========                                                                                                                |   9%  |                                                                                                                                   |============                                                                                                               |   9%  |                                                                                                                                   |============                                                                                                               |  10%  |                                                                                                                                   |=============                                                                                                              |  10%  |                                                                                                                                   |=============                                                                                                              |  11%  |                                                                                                                                   |==============                                                                                                             |  11%  |                                                                                                                                   |==============                                                                                                             |  12%  |                                                                                                                                   |===============                                                                                                            |  12%  |                                                                                                                                   |===============                                                                                                            |  13%  |                                                                                                                                   |================                                                                                                           |  13%  |                                                                                                                                   |=================                                                                                                          |  13%  |                                                                                                                                   |=================                                                                                                          |  14%  |                                                                                                                                   |==================                                                                                                         |  14%  |                                                                                                                                   |==================                                                                                                         |  15%  |                                                                                                                                   |===================                                                                                                        |  15%  |                                                                                                                                   |===================                                                                                                        |  16%  |                                                                                                                                   |====================                                                                                                       |  16%  |                                                                                                                                   |====================                                                                                                       |  17%  |                                                                                                                                   |=====================                                                                                                      |  17%  |                                                                                                                                   |======================                                                                                                     |  17%  |                                                                                                                                   |======================                                                                                                     |  18%  |                                                                                                                                   |=======================                                                                                                    |  18%  |                                                                                                                                   |=======================                                                                                                    |  19%  |                                                                                                                                   |========================                                                                                                   |  19%  |                                                                                                                                   |========================                                                                                                   |  20%  |                                                                                                                                   |=========================                                                                                                  |  20%  |                                                                                                                                   |=========================                                                                                                  |  21%  |                                                                                                                                   |==========================                                                                                                 |  21%  |                                                                                                                                   |==========================                                                                                                 |  22%  |                                                                                                                                   |===========================                                                                                                |  22%  |                                                                                                                                   |============================                                                                                               |  22%  |                                                                                                                                   |============================                                                                                               |  23%  |                                                                                                                                   |=============================                                                                                              |  23%  |                                                                                                                                   |=============================                                                                                              |  24%  |                                                                                                                                   |==============================                                                                                             |  24%  |                                                                                                                                   |==============================                                                                                             |  25%  |                                                                                                                                   |===============================                                                                                            |  25%  |                                                                                                                                   |===============================                                                                                            |  26%  |                                                                                                                                   |================================                                                                                           |  26%  |                                                                                                                                   |=================================                                                                                          |  26%  |                                                                                                                                   |=================================                                                                                          |  27%  |                                                                                                                                   |==================================                                                                                         |  27%  |                                                                                                                                   |==================================                                                                                         |  28%  |                                                                                                                                   |===================================                                                                                        |  28%  |                                                                                                                                   |===================================                                                                                        |  29%  |                                                                                                                                   |====================================                                                                                       |  29%  |                                                                                                                                   |====================================                                                                                       |  30%  |                                                                                                                                   |=====================================                                                                                      |  30%  |                                                                                                                                   |======================================                                                                                     |  30%  |                                                                                                                                   |======================================                                                                                     |  31%  |                                                                                                                                   |=======================================                                                                                    |  31%  |                                                                                                                                   |=======================================                                                                                    |  32%  |                                                                                                                                   |========================================                                                                                   |  32%  |                                                                                                                                   |========================================                                                                                   |  33%  |                                                                                                                                   |=========================================                                                                                  |  33%  |                                                                                                                                   |=========================================                                                                                  |  34%  |                                                                                                                                   |==========================================                                                                                 |  34%  |                                                                                                                                   |==========================================                                                                                 |  35%  |                                                                                                                                   |===========================================                                                                                |  35%  |                                                                                                                                   |============================================                                                                               |  35%  |                                                                                                                                   |============================================                                                                               |  36%  |                                                                                                                                   |=============================================                                                                              |  36%  |                                                                                                                                   |=============================================                                                                              |  37%  |                                                                                                                                   |==============================================                                                                             |  37%  |                                                                                                                                   |==============================================                                                                             |  38%  |                                                                                                                                   |===============================================                                                                            |  38%  |                                                                                                                                   |===============================================                                                                            |  39%  |                                                                                                                                   |================================================                                                                           |  39%  |                                                                                                                                   |=================================================                                                                          |  39%  |                                                                                                                                   |=================================================                                                                          |  40%  |                                                                                                                                   |==================================================                                                                         |  40%  |                                                                                                                                   |==================================================                                                                         |  41%  |                                                                                                                                   |===================================================                                                                        |  41%  |                                                                                                                                   |===================================================                                                                        |  42%  |                                                                                                                                   |====================================================                                                                       |  42%  |                                                                                                                                   |====================================================                                                                       |  43%  |                                                                                                                                   |=====================================================                                                                      |  43%  |                                                                                                                                   |======================================================                                                                     |  44%  |                                                                                                                                   |=======================================================                                                                    |  44%  |                                                                                                                                   |=======================================================                                                                    |  45%  |                                                                                                                                   |========================================================                                                                   |  45%  |                                                                                                                                   |========================================================                                                                   |  46%  |                                                                                                                                   |=========================================================                                                                  |  46%  |                                                                                                                                   |=========================================================                                                                  |  47%  |                                                                                                                                   |==========================================================                                                                 |  47%  |                                                                                                                                   |==========================================================                                                                 |  48%  |                                                                                                                                   |===========================================================                                                                |  48%  |                                                                                                                                   |============================================================                                                               |  48%  |                                                                                                                                   |============================================================                                                               |  49%  |                                                                                                                                   |=============================================================                                                              |  49%  |                                                                                                                                   |=============================================================                                                              |  50%  |                                                                                                                                   |==============================================================                                                             |  50%  |                                                                                                                                   |==============================================================                                                             |  51%  |                                                                                                                                   |===============================================================                                                            |  51%  |                                                                                                                                   |===============================================================                                                            |  52%  |                                                                                                                                   |================================================================                                                           |  52%  |                                                                                                                                   |=================================================================                                                          |  52%  |                                                                                                                                   |=================================================================                                                          |  53%  |                                                                                                                                   |==================================================================                                                         |  53%  |                                                                                                                                   |==================================================================                                                         |  54%  |                                                                                                                                   |===================================================================                                                        |  54%  |                                                                                                                                   |===================================================================                                                        |  55%  |                                                                                                                                   |====================================================================                                                       |  55%  |                                                                                                                                   |====================================================================                                                       |  56%  |                                                                                                                                   |=====================================================================                                                      |  56%  |                                                                                                                                   |======================================================================                                                     |  57%  |                                                                                                                                   |=======================================================================                                                    |  57%  |                                                                                                                                   |=======================================================================                                                    |  58%  |                                                                                                                                   |========================================================================                                                   |  58%  |                                                                                                                                   |========================================================================                                                   |  59%  |                                                                                                                                   |=========================================================================                                                  |  59%  |                                                                                                                                   |=========================================================================                                                  |  60%  |                                                                                                                                   |==========================================================================                                                 |  60%  |                                                                                                                                   |==========================================================================                                                 |  61%  |                                                                                                                                   |===========================================================================                                                |  61%  |                                                                                                                                   |============================================================================                                               |  61%  |                                                                                                                                   |============================================================================                                               |  62%  |                                                                                                                                   |=============================================================================                                              |  62%  |                                                                                                                                   |=============================================================================                                              |  63%  |                                                                                                                                   |==============================================================================                                             |  63%  |                                                                                                                                   |==============================================================================                                             |  64%  |                                                                                                                                   |===============================================================================                                            |  64%  |                                                                                                                                   |===============================================================================                                            |  65%  |                                                                                                                                   |================================================================================                                           |  65%  |                                                                                                                                   |=================================================================================                                          |  65%  |                                                                                                                                   |=================================================================================                                          |  66%  |                                                                                                                                   |==================================================================================                                         |  66%  |                                                                                                                                   |==================================================================================                                         |  67%  |                                                                                                                                   |===================================================================================                                        |  67%  |                                                                                                                                   |===================================================================================                                        |  68%  |                                                                                                                                   |====================================================================================                                       |  68%  |                                                                                                                                   |====================================================================================                                       |  69%  |                                                                                                                                   |=====================================================================================                                      |  69%  |                                                                                                                                   |=====================================================================================                                      |  70%  |                                                                                                                                   |======================================================================================                                     |  70%  |                                                                                                                                   |=======================================================================================                                    |  70%  |                                                                                                                                   |=======================================================================================                                    |  71%  |                                                                                                                                   |========================================================================================                                   |  71%  |                                                                                                                                   |========================================================================================                                   |  72%  |                                                                                                                                   |=========================================================================================                                  |  72%  |                                                                                                                                   |=========================================================================================                                  |  73%  |                                                                                                                                   |==========================================================================================                                 |  73%  |                                                                                                                                   |==========================================================================================                                 |  74%  |                                                                                                                                   |===========================================================================================                                |  74%  |                                                                                                                                   |============================================================================================                               |  74%  |                                                                                                                                   |============================================================================================                               |  75%  |                                                                                                                                   |=============================================================================================                              |  75%  |                                                                                                                                   |=============================================================================================                              |  76%  |                                                                                                                                   |==============================================================================================                             |  76%  |                                                                                                                                   |==============================================================================================                             |  77%  |                                                                                                                                   |===============================================================================================                            |  77%  |                                                                                                                                   |===============================================================================================                            |  78%  |                                                                                                                                   |================================================================================================                           |  78%  |                                                                                                                                   |=================================================================================================                          |  78%  |                                                                                                                                   |=================================================================================================                          |  79%  |                                                                                                                                   |==================================================================================================                         |  79%  |                                                                                                                                   |==================================================================================================                         |  80%  |                                                                                                                                   |===================================================================================================                        |  80%  |                                                                                                                                   |===================================================================================================                        |  81%  |                                                                                                                                   |====================================================================================================                       |  81%  |                                                                                                                                   |====================================================================================================                       |  82%  |                                                                                                                                   |=====================================================================================================                      |  82%  |                                                                                                                                   |=====================================================================================================                      |  83%  |                                                                                                                                   |======================================================================================================                     |  83%  |                                                                                                                                   |=======================================================================================================                    |  83%  |                                                                                                                                   |=======================================================================================================                    |  84%  |                                                                                                                                   |========================================================================================================                   |  84%  |                                                                                                                                   |========================================================================================================                   |  85%  |                                                                                                                                   |=========================================================================================================                  |  85%  |                                                                                                                                   |=========================================================================================================                  |  86%  |                                                                                                                                   |==========================================================================================================                 |  86%  |                                                                                                                                   |==========================================================================================================                 |  87%  |                                                                                                                                   |===========================================================================================================                |  87%  |                                                                                                                                   |============================================================================================================               |  87%  |                                                                                                                                   |============================================================================================================               |  88%  |                                                                                                                                   |=============================================================================================================              |  88%  |                                                                                                                                   |=============================================================================================================              |  89%  |                                                                                                                                   |==============================================================================================================             |  89%  |                                                                                                                                   |==============================================================================================================             |  90%  |                                                                                                                                   |===============================================================================================================            |  90%  |                                                                                                                                   |===============================================================================================================            |  91%  |                                                                                                                                   |================================================================================================================           |  91%  |                                                                                                                                   |=================================================================================================================          |  91%  |                                                                                                                                   |=================================================================================================================          |  92%  |                                                                                                                                   |==================================================================================================================         |  92%  |                                                                                                                                   |==================================================================================================================         |  93%  |                                                                                                                                   |===================================================================================================================        |  93%  |                                                                                                                                   |===================================================================================================================        |  94%  |                                                                                                                                   |====================================================================================================================       |  94%  |                                                                                                                                   |====================================================================================================================       |  95%  |                                                                                                                                   |=====================================================================================================================      |  95%  |                                                                                                                                   |=====================================================================================================================      |  96%  |                                                                                                                                   |======================================================================================================================     |  96%  |                                                                                                                                   |=======================================================================================================================    |  96%  |                                                                                                                                   |=======================================================================================================================    |  97%  |                                                                                                                                   |========================================================================================================================   |  97%  |                                                                                                                                   |========================================================================================================================   |  98%  |                                                                                                                                   |=========================================================================================================================  |  98%  |                                                                                                                                   |=========================================================================================================================  |  99%  |                                                                                                                                   |========================================================================================================================== |  99%  |                                                                                                                                   |========================================================================================================================== | 100%  |                                                                                                                                   |===========================================================================================================================| 100%
```

```
## Removing confirmed isotope entries:
```

```
## Recalculating peaklist for plotting:
```

```
##   |                                                                                                                                   |                                                                                                                           |   0%  |                                                                                                                                   |===============================                                                                                            |  25%  |                                                                                                                                   |==============================================================                                                             |  50%  |                                                                                                                                   |============================================================================================                               |  75%  |                                                                                                                                   |===========================================================================================================================| 100%
```

```
## Calculating adduct similarity:
```

```
##   |                                                                                                                                   |                                                                                                                           |   0%  |                                                                                                                                   |===============================                                                                                            |  25%  |                                                                                                                                   |==============================================================                                                             |  50%  |                                                                                                                                   |============================================================================================                               |  75%  |                                                                                                                                   |===========================================================================================================================| 100%
```

```
## 100 out of 3530
```

```
## 200 out of 3530
```

```
## 300 out of 3530
```

```
## 400 out of 3530
```

```
## 500 out of 3530
```

```
## 600 out of 3530
```

```
## 700 out of 3530
```

```
## 800 out of 3530
```

```
## 900 out of 3530
```

```
## 1000 out of 3530
```

```
## 1100 out of 3530
```

```
## 1200 out of 3530
```

```
## 1300 out of 3530
```

```
## 1400 out of 3530
```

```
## 1500 out of 3530
```

```
## 1600 out of 3530
```

```
## 1700 out of 3530
```

```
## 1800 out of 3530
```

```
## 1900 out of 3530
```

```
## 2000 out of 3530
```

```
## 2100 out of 3530
```

```
## 2200 out of 3530
```

```
## 2300 out of 3530
```

```
## 2400 out of 3530
```

```
## 2500 out of 3530
```

```
## 2600 out of 3530
```

```
## 2700 out of 3530
```

```
## 2800 out of 3530
```

```
## 2900 out of 3530
```

```
## 3000 out of 3530
```

```
## 3100 out of 3530
```

```
## 3200 out of 3530
```

```
## 3300 out of 3530
```

```
## 3400 out of 3530
```

```
## 3500 out of 3530
```

```
##   |                                                                                                                                   |                                                                                                                           |   0%  |                                                                                                                                   |===============================                                                                                            |  25%  |                                                                                                                                   |==============================================================                                                             |  50%  |                                                                                                                                   |============================================================================================                               |  75%  |                                                                                                                                   |===========================================================================================================================| 100%
```

```
## Computing isotopic elution similarity:
```

```
## Calculating isotope similarity:
```

```
##   |                                                                                                                                   |                                                                                                                           |   0%  |                                                                                                                                   |===============================                                                                                            |  25%  |                                                                                                                                   |==============================================================                                                             |  50%  |                                                                                                                                   |============================================================================================                               |  75%  |                                                                                                                                   |===========================================================================================================================| 100%
```

```
## Calculating isotopic fidelity metrics:
```

```
##   |                                                                                                                                   |                                                                                                                           |   0%  |                                                                                                                                   |=                                                                                                                          |   0%  |                                                                                                                                   |=                                                                                                                          |   1%  |                                                                                                                                   |==                                                                                                                         |   1%  |                                                                                                                                   |==                                                                                                                         |   2%  |                                                                                                                                   |===                                                                                                                        |   2%  |                                                                                                                                   |===                                                                                                                        |   3%  |                                                                                                                                   |====                                                                                                                       |   3%  |                                                                                                                                   |====                                                                                                                       |   4%  |                                                                                                                                   |=====                                                                                                                      |   4%  |                                                                                                                                   |======                                                                                                                     |   4%  |                                                                                                                                   |======                                                                                                                     |   5%  |                                                                                                                                   |=======                                                                                                                    |   5%  |                                                                                                                                   |=======                                                                                                                    |   6%  |                                                                                                                                   |========                                                                                                                   |   6%  |                                                                                                                                   |========                                                                                                                   |   7%  |                                                                                                                                   |=========                                                                                                                  |   7%  |                                                                                                                                   |=========                                                                                                                  |   8%  |                                                                                                                                   |==========                                                                                                                 |   8%  |                                                                                                                                   |==========                                                                                                                 |   9%  |                                                                                                                                   |===========                                                                                                                |   9%  |                                                                                                                                   |============                                                                                                               |   9%  |                                                                                                                                   |============                                                                                                               |  10%  |                                                                                                                                   |=============                                                                                                              |  10%  |                                                                                                                                   |=============                                                                                                              |  11%  |                                                                                                                                   |==============                                                                                                             |  11%  |                                                                                                                                   |==============                                                                                                             |  12%  |                                                                                                                                   |===============                                                                                                            |  12%  |                                                                                                                                   |===============                                                                                                            |  13%  |                                                                                                                                   |================                                                                                                           |  13%  |                                                                                                                                   |=================                                                                                                          |  13%  |                                                                                                                                   |=================                                                                                                          |  14%  |                                                                                                                                   |==================                                                                                                         |  14%  |                                                                                                                                   |==================                                                                                                         |  15%  |                                                                                                                                   |===================                                                                                                        |  15%  |                                                                                                                                   |===================                                                                                                        |  16%  |                                                                                                                                   |====================                                                                                                       |  16%  |                                                                                                                                   |====================                                                                                                       |  17%  |                                                                                                                                   |=====================                                                                                                      |  17%  |                                                                                                                                   |======================                                                                                                     |  17%  |                                                                                                                                   |======================                                                                                                     |  18%  |                                                                                                                                   |=======================                                                                                                    |  18%  |                                                                                                                                   |=======================                                                                                                    |  19%  |                                                                                                                                   |========================                                                                                                   |  19%  |                                                                                                                                   |========================                                                                                                   |  20%  |                                                                                                                                   |=========================                                                                                                  |  20%  |                                                                                                                                   |=========================                                                                                                  |  21%  |                                                                                                                                   |==========================                                                                                                 |  21%  |                                                                                                                                   |==========================                                                                                                 |  22%  |                                                                                                                                   |===========================                                                                                                |  22%  |                                                                                                                                   |============================                                                                                               |  22%  |                                                                                                                                   |============================                                                                                               |  23%  |                                                                                                                                   |=============================                                                                                              |  23%  |                                                                                                                                   |=============================                                                                                              |  24%  |                                                                                                                                   |==============================                                                                                             |  24%  |                                                                                                                                   |==============================                                                                                             |  25%  |                                                                                                                                   |===============================                                                                                            |  25%  |                                                                                                                                   |===============================                                                                                            |  26%  |                                                                                                                                   |================================                                                                                           |  26%  |                                                                                                                                   |=================================                                                                                          |  26%  |                                                                                                                                   |=================================                                                                                          |  27%  |                                                                                                                                   |==================================                                                                                         |  27%  |                                                                                                                                   |==================================                                                                                         |  28%  |                                                                                                                                   |===================================                                                                                        |  28%  |                                                                                                                                   |===================================                                                                                        |  29%  |                                                                                                                                   |====================================                                                                                       |  29%  |                                                                                                                                   |====================================                                                                                       |  30%  |                                                                                                                                   |=====================================                                                                                      |  30%  |                                                                                                                                   |======================================                                                                                     |  30%  |                                                                                                                                   |======================================                                                                                     |  31%  |                                                                                                                                   |=======================================                                                                                    |  31%  |                                                                                                                                   |=======================================                                                                                    |  32%  |                                                                                                                                   |========================================                                                                                   |  32%  |                                                                                                                                   |========================================                                                                                   |  33%  |                                                                                                                                   |=========================================                                                                                  |  33%  |                                                                                                                                   |=========================================                                                                                  |  34%  |                                                                                                                                   |==========================================                                                                                 |  34%  |                                                                                                                                   |==========================================                                                                                 |  35%  |                                                                                                                                   |===========================================                                                                                |  35%  |                                                                                                                                   |============================================                                                                               |  35%  |                                                                                                                                   |============================================                                                                               |  36%  |                                                                                                                                   |=============================================                                                                              |  36%  |                                                                                                                                   |=============================================                                                                              |  37%  |                                                                                                                                   |==============================================                                                                             |  37%  |                                                                                                                                   |==============================================                                                                             |  38%  |                                                                                                                                   |===============================================                                                                            |  38%  |                                                                                                                                   |===============================================                                                                            |  39%  |                                                                                                                                   |================================================                                                                           |  39%  |                                                                                                                                   |=================================================                                                                          |  39%  |                                                                                                                                   |=================================================                                                                          |  40%  |                                                                                                                                   |==================================================                                                                         |  40%  |                                                                                                                                   |==================================================                                                                         |  41%  |                                                                                                                                   |===================================================                                                                        |  41%  |                                                                                                                                   |===================================================                                                                        |  42%  |                                                                                                                                   |====================================================                                                                       |  42%  |                                                                                                                                   |====================================================                                                                       |  43%  |                                                                                                                                   |=====================================================                                                                      |  43%  |                                                                                                                                   |======================================================                                                                     |  43%  |                                                                                                                                   |======================================================                                                                     |  44%  |                                                                                                                                   |=======================================================                                                                    |  44%  |                                                                                                                                   |=======================================================                                                                    |  45%  |                                                                                                                                   |========================================================                                                                   |  45%  |                                                                                                                                   |========================================================                                                                   |  46%  |                                                                                                                                   |=========================================================                                                                  |  46%  |                                                                                                                                   |=========================================================                                                                  |  47%  |                                                                                                                                   |==========================================================                                                                 |  47%  |                                                                                                                                   |==========================================================                                                                 |  48%  |                                                                                                                                   |===========================================================                                                                |  48%  |                                                                                                                                   |============================================================                                                               |  48%  |                                                                                                                                   |============================================================                                                               |  49%  |                                                                                                                                   |=============================================================                                                              |  49%  |                                                                                                                                   |=============================================================                                                              |  50%  |                                                                                                                                   |==============================================================                                                             |  50%  |                                                                                                                                   |==============================================================                                                             |  51%  |                                                                                                                                   |===============================================================                                                            |  51%  |                                                                                                                                   |===============================================================                                                            |  52%  |                                                                                                                                   |================================================================                                                           |  52%  |                                                                                                                                   |=================================================================                                                          |  52%  |                                                                                                                                   |=================================================================                                                          |  53%  |                                                                                                                                   |==================================================================                                                         |  53%  |                                                                                                                                   |==================================================================                                                         |  54%  |                                                                                                                                   |===================================================================                                                        |  54%  |                                                                                                                                   |===================================================================                                                        |  55%  |                                                                                                                                   |====================================================================                                                       |  55%  |                                                                                                                                   |====================================================================                                                       |  56%  |                                                                                                                                   |=====================================================================                                                      |  56%  |                                                                                                                                   |=====================================================================                                                      |  57%  |                                                                                                                                   |======================================================================                                                     |  57%  |                                                                                                                                   |=======================================================================                                                    |  57%  |                                                                                                                                   |=======================================================================                                                    |  58%  |                                                                                                                                   |========================================================================                                                   |  58%  |                                                                                                                                   |========================================================================                                                   |  59%  |                                                                                                                                   |=========================================================================                                                  |  59%  |                                                                                                                                   |=========================================================================                                                  |  60%  |                                                                                                                                   |==========================================================================                                                 |  60%  |                                                                                                                                   |==========================================================================                                                 |  61%  |                                                                                                                                   |===========================================================================                                                |  61%  |                                                                                                                                   |============================================================================                                               |  61%  |                                                                                                                                   |============================================================================                                               |  62%  |                                                                                                                                   |=============================================================================                                              |  62%  |                                                                                                                                   |=============================================================================                                              |  63%  |                                                                                                                                   |==============================================================================                                             |  63%  |                                                                                                                                   |==============================================================================                                             |  64%  |                                                                                                                                   |===============================================================================                                            |  64%  |                                                                                                                                   |===============================================================================                                            |  65%  |                                                                                                                                   |================================================================================                                           |  65%  |                                                                                                                                   |=================================================================================                                          |  65%  |                                                                                                                                   |=================================================================================                                          |  66%  |                                                                                                                                   |==================================================================================                                         |  66%  |                                                                                                                                   |==================================================================================                                         |  67%  |                                                                                                                                   |===================================================================================                                        |  67%  |                                                                                                                                   |===================================================================================                                        |  68%  |                                                                                                                                   |====================================================================================                                       |  68%  |                                                                                                                                   |====================================================================================                                       |  69%  |                                                                                                                                   |=====================================================================================                                      |  69%  |                                                                                                                                   |=====================================================================================                                      |  70%  |                                                                                                                                   |======================================================================================                                     |  70%  |                                                                                                                                   |=======================================================================================                                    |  70%  |                                                                                                                                   |=======================================================================================                                    |  71%  |                                                                                                                                   |========================================================================================                                   |  71%  |                                                                                                                                   |========================================================================================                                   |  72%  |                                                                                                                                   |=========================================================================================                                  |  72%  |                                                                                                                                   |=========================================================================================                                  |  73%  |                                                                                                                                   |==========================================================================================                                 |  73%  |                                                                                                                                   |==========================================================================================                                 |  74%  |                                                                                                                                   |===========================================================================================                                |  74%  |                                                                                                                                   |============================================================================================                               |  74%  |                                                                                                                                   |============================================================================================                               |  75%  |                                                                                                                                   |=============================================================================================                              |  75%  |                                                                                                                                   |=============================================================================================                              |  76%  |                                                                                                                                   |==============================================================================================                             |  76%  |                                                                                                                                   |==============================================================================================                             |  77%  |                                                                                                                                   |===============================================================================================                            |  77%  |                                                                                                                                   |===============================================================================================                            |  78%  |                                                                                                                                   |================================================================================================                           |  78%  |                                                                                                                                   |=================================================================================================                          |  78%  |                                                                                                                                   |=================================================================================================                          |  79%  |                                                                                                                                   |==================================================================================================                         |  79%  |                                                                                                                                   |==================================================================================================                         |  80%  |                                                                                                                                   |===================================================================================================                        |  80%  |                                                                                                                                   |===================================================================================================                        |  81%  |                                                                                                                                   |====================================================================================================                       |  81%  |                                                                                                                                   |====================================================================================================                       |  82%  |                                                                                                                                   |=====================================================================================================                      |  82%  |                                                                                                                                   |=====================================================================================================                      |  83%  |                                                                                                                                   |======================================================================================================                     |  83%  |                                                                                                                                   |=======================================================================================================                    |  83%  |                                                                                                                                   |=======================================================================================================                    |  84%  |                                                                                                                                   |========================================================================================================                   |  84%  |                                                                                                                                   |========================================================================================================                   |  85%  |                                                                                                                                   |=========================================================================================================                  |  85%  |                                                                                                                                   |=========================================================================================================                  |  86%  |                                                                                                                                   |==========================================================================================================                 |  86%  |                                                                                                                                   |==========================================================================================================                 |  87%  |                                                                                                                                   |===========================================================================================================                |  87%  |                                                                                                                                   |============================================================================================================               |  87%  |                                                                                                                                   |============================================================================================================               |  88%  |                                                                                                                                   |=============================================================================================================              |  88%  |                                                                                                                                   |=============================================================================================================              |  89%  |                                                                                                                                   |==============================================================================================================             |  89%  |                                                                                                                                   |==============================================================================================================             |  90%  |                                                                                                                                   |===============================================================================================================            |  90%  |                                                                                                                                   |===============================================================================================================            |  91%  |                                                                                                                                   |================================================================================================================           |  91%  |                                                                                                                                   |=================================================================================================================          |  91%  |                                                                                                                                   |=================================================================================================================          |  92%  |                                                                                                                                   |==================================================================================================================         |  92%  |                                                                                                                                   |==================================================================================================================         |  93%  |                                                                                                                                   |===================================================================================================================        |  93%  |                                                                                                                                   |===================================================================================================================        |  94%  |                                                                                                                                   |====================================================================================================================       |  94%  |                                                                                                                                   |====================================================================================================================       |  95%  |                                                                                                                                   |=====================================================================================================================      |  95%  |                                                                                                                                   |=====================================================================================================================      |  96%  |                                                                                                                                   |======================================================================================================================     |  96%  |                                                                                                                                   |=======================================================================================================================    |  96%  |                                                                                                                                   |=======================================================================================================================    |  97%  |                                                                                                                                   |========================================================================================================================   |  97%  |                                                                                                                                   |========================================================================================================================   |  98%  |                                                                                                                                   |=========================================================================================================================  |  98%  |                                                                                                                                   |=========================================================================================================================  |  99%  |                                                                                                                                   |========================================================================================================================== |  99%  |                                                                                                                                   |========================================================================================================================== | 100%  |                                                                                                                                   |===========================================================================================================================| 100%
```

```
## Removing confirmed isotope entries:
```

```
## Recalculating peaklist for plotting:
```

```
##   |                                                                                                                                   |                                                                                                                           |   0%  |                                                                                                                                   |===============================                                                                                            |  25%  |                                                                                                                                   |==============================================================                                                             |  50%  |                                                                                                                                   |============================================================================================                               |  75%  |                                                                                                                                   |===========================================================================================================================| 100%
```

```
## Calculating adduct similarity:
```

```
##   |                                                                                                                                   |                                                                                                                           |   0%  |                                                                                                                                   |===============================                                                                                            |  25%  |                                                                                                                                   |==============================================================                                                             |  50%  |                                                                                                                                   |============================================================================================                               |  75%  |                                                                                                                                   |===========================================================================================================================| 100%
```

```
## 100 out of 14666
```

```
## 200 out of 14666
```

```
## 300 out of 14666
```

```
## 400 out of 14666
```

```
## 500 out of 14666
```

```
## 600 out of 14666
```

```
## 700 out of 14666
```

```
## 800 out of 14666
```

```
## 900 out of 14666
```

```
## 1000 out of 14666
```

```
## 1100 out of 14666
```

```
## 1200 out of 14666
```

```
## 1300 out of 14666
```

```
## 1400 out of 14666
```

```
## 1500 out of 14666
```

```
## 1600 out of 14666
```

```
## 1700 out of 14666
```

```
## 1800 out of 14666
```

```
## 1900 out of 14666
```

```
## 2000 out of 14666
```

```
## 2100 out of 14666
```

```
## 2200 out of 14666
```

```
## 2300 out of 14666
```

```
## 2400 out of 14666
```

```
## 2500 out of 14666
```

```
## 2600 out of 14666
```

```
## 2700 out of 14666
```

```
## 2800 out of 14666
```

```
## 2900 out of 14666
```

```
## 3000 out of 14666
```

```
## 3100 out of 14666
```

```
## 3200 out of 14666
```

```
## 3300 out of 14666
```

```
## 3400 out of 14666
```

```
## 3500 out of 14666
```

```
## 3600 out of 14666
```

```
## 3700 out of 14666
```

```
## 3800 out of 14666
```

```
## 3900 out of 14666
```

```
## 4000 out of 14666
```

```
## 4100 out of 14666
```

```
## 4200 out of 14666
```

```
## 4300 out of 14666
```

```
## 4400 out of 14666
```

```
## 4500 out of 14666
```

```
## 4600 out of 14666
```

```
## 4700 out of 14666
```

```
## 4800 out of 14666
```

```
## 4900 out of 14666
```

```
## 5000 out of 14666
```

```
## 5100 out of 14666
```

```
## 5200 out of 14666
```

```
## 5300 out of 14666
```

```
## 5400 out of 14666
```

```
## 5500 out of 14666
```

```
## 5600 out of 14666
```

```
## 5700 out of 14666
```

```
## 5800 out of 14666
```

```
## 5900 out of 14666
```

```
## 6000 out of 14666
```

```
## 6100 out of 14666
```

```
## 6200 out of 14666
```

```
## 6300 out of 14666
```

```
## 6400 out of 14666
```

```
## 6500 out of 14666
```

```
## 6600 out of 14666
```

```
## 6700 out of 14666
```

```
## 6800 out of 14666
```

```
## 6900 out of 14666
```

```
## 7000 out of 14666
```

```
## 7100 out of 14666
```

```
## 7200 out of 14666
```

```
## 7300 out of 14666
```

```
## 7400 out of 14666
```

```
## 7500 out of 14666
```

```
## 7600 out of 14666
```

```
## 7700 out of 14666
```

```
## 7800 out of 14666
```

```
## 7900 out of 14666
```

```
## 8000 out of 14666
```

```
## 8100 out of 14666
```

```
## 8200 out of 14666
```

```
## 8300 out of 14666
```

```
## 8400 out of 14666
```

```
## 8500 out of 14666
```

```
## 8600 out of 14666
```

```
## 8700 out of 14666
```

```
## 8800 out of 14666
```

```
## 8900 out of 14666
```

```
## 9000 out of 14666
```

```
## 9100 out of 14666
```

```
## 9200 out of 14666
```

```
## 9300 out of 14666
```

```
## 9400 out of 14666
```

```
## 9500 out of 14666
```

```
## 9600 out of 14666
```

```
## 9700 out of 14666
```

```
## 9800 out of 14666
```

```
## 9900 out of 14666
```

```
## 10000 out of 14666
```

```
## 10100 out of 14666
```

```
## 10200 out of 14666
```

```
## 10300 out of 14666
```

```
## 10400 out of 14666
```

```
## 10500 out of 14666
```

```
## 10600 out of 14666
```

```
## 10700 out of 14666
```

```
## 10800 out of 14666
```

```
## 10900 out of 14666
```

```
## 11000 out of 14666
```

```
## 11100 out of 14666
```

```
## 11200 out of 14666
```

```
## 11300 out of 14666
```

```
## 11400 out of 14666
```

```
## 11500 out of 14666
```

```
## 11600 out of 14666
```

```
## 11700 out of 14666
```

```
## 11800 out of 14666
```

```
## 11900 out of 14666
```

```
## 12000 out of 14666
```

```
## 12100 out of 14666
```

```
## 12200 out of 14666
```

```
## 12300 out of 14666
```

```
## 12400 out of 14666
```

```
## 12500 out of 14666
```

```
## 12600 out of 14666
```

```
## 12700 out of 14666
```

```
## 12800 out of 14666
```

```
## 12900 out of 14666
```

```
## 13000 out of 14666
```

```
## 13100 out of 14666
```

```
## 13200 out of 14666
```

```
## 13300 out of 14666
```

```
## 13400 out of 14666
```

```
## 13500 out of 14666
```

```
## 13600 out of 14666
```

```
## 13700 out of 14666
```

```
## 13800 out of 14666
```

```
## 13900 out of 14666
```

```
## 14000 out of 14666
```

```
## 14100 out of 14666
```

```
## 14200 out of 14666
```

```
## 14300 out of 14666
```

```
## 14400 out of 14666
```

```
## 14500 out of 14666
```

```
## 14600 out of 14666
```

```
##   |                                                                                                                                   |                                                                                                                           |   0%  |                                                                                                                                   |===============================                                                                                            |  25%  |                                                                                                                                   |==============================================================                                                             |  50%  |                                                                                                                                   |============================================================================================                               |  75%  |                                                                                                                                   |===========================================================================================================================| 100%
```

```
## Computing isotopic elution similarity:
```

```
## Calculating isotope similarity:
```

```
##   |                                                                                                                                   |                                                                                                                           |   0%  |                                                                                                                                   |===============================                                                                                            |  25%  |                                                                                                                                   |==============================================================                                                             |  50%  |                                                                                                                                   |============================================================================================                               |  75%  |                                                                                                                                   |===========================================================================================================================| 100%
```

```
## Calculating isotopic fidelity metrics:
```

```
##   |                                                                                                                                   |                                                                                                                           |   0%  |                                                                                                                                   |=                                                                                                                          |   0%  |                                                                                                                                   |=                                                                                                                          |   1%  |                                                                                                                                   |==                                                                                                                         |   1%  |                                                                                                                                   |==                                                                                                                         |   2%  |                                                                                                                                   |===                                                                                                                        |   2%  |                                                                                                                                   |===                                                                                                                        |   3%  |                                                                                                                                   |====                                                                                                                       |   3%  |                                                                                                                                   |====                                                                                                                       |   4%  |                                                                                                                                   |=====                                                                                                                      |   4%  |                                                                                                                                   |======                                                                                                                     |   4%  |                                                                                                                                   |======                                                                                                                     |   5%  |                                                                                                                                   |=======                                                                                                                    |   5%  |                                                                                                                                   |=======                                                                                                                    |   6%  |                                                                                                                                   |========                                                                                                                   |   6%  |                                                                                                                                   |========                                                                                                                   |   7%  |                                                                                                                                   |=========                                                                                                                  |   7%  |                                                                                                                                   |=========                                                                                                                  |   8%  |                                                                                                                                   |==========                                                                                                                 |   8%  |                                                                                                                                   |==========                                                                                                                 |   9%  |                                                                                                                                   |===========                                                                                                                |   9%  |                                                                                                                                   |============                                                                                                               |   9%  |                                                                                                                                   |============                                                                                                               |  10%  |                                                                                                                                   |=============                                                                                                              |  10%  |                                                                                                                                   |=============                                                                                                              |  11%  |                                                                                                                                   |==============                                                                                                             |  11%  |                                                                                                                                   |==============                                                                                                             |  12%  |                                                                                                                                   |===============                                                                                                            |  12%  |                                                                                                                                   |===============                                                                                                            |  13%  |                                                                                                                                   |================                                                                                                           |  13%  |                                                                                                                                   |=================                                                                                                          |  13%  |                                                                                                                                   |=================                                                                                                          |  14%  |                                                                                                                                   |==================                                                                                                         |  14%  |                                                                                                                                   |==================                                                                                                         |  15%  |                                                                                                                                   |===================                                                                                                        |  15%  |                                                                                                                                   |===================                                                                                                        |  16%  |                                                                                                                                   |====================                                                                                                       |  16%  |                                                                                                                                   |====================                                                                                                       |  17%  |                                                                                                                                   |=====================                                                                                                      |  17%  |                                                                                                                                   |======================                                                                                                     |  17%  |                                                                                                                                   |======================                                                                                                     |  18%  |                                                                                                                                   |=======================                                                                                                    |  18%  |                                                                                                                                   |=======================                                                                                                    |  19%  |                                                                                                                                   |========================                                                                                                   |  19%  |                                                                                                                                   |========================                                                                                                   |  20%  |                                                                                                                                   |=========================                                                                                                  |  20%  |                                                                                                                                   |=========================                                                                                                  |  21%  |                                                                                                                                   |==========================                                                                                                 |  21%  |                                                                                                                                   |==========================                                                                                                 |  22%  |                                                                                                                                   |===========================                                                                                                |  22%  |                                                                                                                                   |============================                                                                                               |  22%  |                                                                                                                                   |============================                                                                                               |  23%  |                                                                                                                                   |=============================                                                                                              |  23%  |                                                                                                                                   |=============================                                                                                              |  24%  |                                                                                                                                   |==============================                                                                                             |  24%  |                                                                                                                                   |==============================                                                                                             |  25%  |                                                                                                                                   |===============================                                                                                            |  25%  |                                                                                                                                   |===============================                                                                                            |  26%  |                                                                                                                                   |================================                                                                                           |  26%  |                                                                                                                                   |=================================                                                                                          |  26%  |                                                                                                                                   |=================================                                                                                          |  27%  |                                                                                                                                   |==================================                                                                                         |  27%  |                                                                                                                                   |==================================                                                                                         |  28%  |                                                                                                                                   |===================================                                                                                        |  28%  |                                                                                                                                   |===================================                                                                                        |  29%  |                                                                                                                                   |====================================                                                                                       |  29%  |                                                                                                                                   |====================================                                                                                       |  30%  |                                                                                                                                   |=====================================                                                                                      |  30%  |                                                                                                                                   |======================================                                                                                     |  30%  |                                                                                                                                   |======================================                                                                                     |  31%  |                                                                                                                                   |=======================================                                                                                    |  31%  |                                                                                                                                   |=======================================                                                                                    |  32%  |                                                                                                                                   |========================================                                                                                   |  32%  |                                                                                                                                   |========================================                                                                                   |  33%  |                                                                                                                                   |=========================================                                                                                  |  33%  |                                                                                                                                   |=========================================                                                                                  |  34%  |                                                                                                                                   |==========================================                                                                                 |  34%  |                                                                                                                                   |==========================================                                                                                 |  35%  |                                                                                                                                   |===========================================                                                                                |  35%  |                                                                                                                                   |============================================                                                                               |  35%  |                                                                                                                                   |============================================                                                                               |  36%  |                                                                                                                                   |=============================================                                                                              |  36%  |                                                                                                                                   |=============================================                                                                              |  37%  |                                                                                                                                   |==============================================                                                                             |  37%  |                                                                                                                                   |==============================================                                                                             |  38%  |                                                                                                                                   |===============================================                                                                            |  38%  |                                                                                                                                   |===============================================                                                                            |  39%  |                                                                                                                                   |================================================                                                                           |  39%  |                                                                                                                                   |=================================================                                                                          |  39%  |                                                                                                                                   |=================================================                                                                          |  40%  |                                                                                                                                   |==================================================                                                                         |  40%  |                                                                                                                                   |==================================================                                                                         |  41%  |                                                                                                                                   |===================================================                                                                        |  41%  |                                                                                                                                   |===================================================                                                                        |  42%  |                                                                                                                                   |====================================================                                                                       |  42%  |                                                                                                                                   |====================================================                                                                       |  43%  |                                                                                                                                   |=====================================================                                                                      |  43%  |                                                                                                                                   |======================================================                                                                     |  44%  |                                                                                                                                   |=======================================================                                                                    |  44%  |                                                                                                                                   |=======================================================                                                                    |  45%  |                                                                                                                                   |========================================================                                                                   |  45%  |                                                                                                                                   |========================================================                                                                   |  46%  |                                                                                                                                   |=========================================================                                                                  |  46%  |                                                                                                                                   |=========================================================                                                                  |  47%  |                                                                                                                                   |==========================================================                                                                 |  47%  |                                                                                                                                   |==========================================================                                                                 |  48%  |                                                                                                                                   |===========================================================                                                                |  48%  |                                                                                                                                   |============================================================                                                               |  48%  |                                                                                                                                   |============================================================                                                               |  49%  |                                                                                                                                   |=============================================================                                                              |  49%  |                                                                                                                                   |=============================================================                                                              |  50%  |                                                                                                                                   |==============================================================                                                             |  50%  |                                                                                                                                   |==============================================================                                                             |  51%  |                                                                                                                                   |===============================================================                                                            |  51%  |                                                                                                                                   |===============================================================                                                            |  52%  |                                                                                                                                   |================================================================                                                           |  52%  |                                                                                                                                   |=================================================================                                                          |  52%  |                                                                                                                                   |=================================================================                                                          |  53%  |                                                                                                                                   |==================================================================                                                         |  53%  |                                                                                                                                   |==================================================================                                                         |  54%  |                                                                                                                                   |===================================================================                                                        |  54%  |                                                                                                                                   |===================================================================                                                        |  55%  |                                                                                                                                   |====================================================================                                                       |  55%  |                                                                                                                                   |====================================================================                                                       |  56%  |                                                                                                                                   |=====================================================================                                                      |  56%  |                                                                                                                                   |======================================================================                                                     |  57%  |                                                                                                                                   |=======================================================================                                                    |  57%  |                                                                                                                                   |=======================================================================                                                    |  58%  |                                                                                                                                   |========================================================================                                                   |  58%  |                                                                                                                                   |========================================================================                                                   |  59%  |                                                                                                                                   |=========================================================================                                                  |  59%  |                                                                                                                                   |=========================================================================                                                  |  60%  |                                                                                                                                   |==========================================================================                                                 |  60%  |                                                                                                                                   |==========================================================================                                                 |  61%  |                                                                                                                                   |===========================================================================                                                |  61%  |                                                                                                                                   |============================================================================                                               |  61%  |                                                                                                                                   |============================================================================                                               |  62%  |                                                                                                                                   |=============================================================================                                              |  62%  |                                                                                                                                   |=============================================================================                                              |  63%  |                                                                                                                                   |==============================================================================                                             |  63%  |                                                                                                                                   |==============================================================================                                             |  64%  |                                                                                                                                   |===============================================================================                                            |  64%  |                                                                                                                                   |===============================================================================                                            |  65%  |                                                                                                                                   |================================================================================                                           |  65%  |                                                                                                                                   |=================================================================================                                          |  65%  |                                                                                                                                   |=================================================================================                                          |  66%  |                                                                                                                                   |==================================================================================                                         |  66%  |                                                                                                                                   |==================================================================================                                         |  67%  |                                                                                                                                   |===================================================================================                                        |  67%  |                                                                                                                                   |===================================================================================                                        |  68%  |                                                                                                                                   |====================================================================================                                       |  68%  |                                                                                                                                   |====================================================================================                                       |  69%  |                                                                                                                                   |=====================================================================================                                      |  69%  |                                                                                                                                   |=====================================================================================                                      |  70%  |                                                                                                                                   |======================================================================================                                     |  70%  |                                                                                                                                   |=======================================================================================                                    |  70%  |                                                                                                                                   |=======================================================================================                                    |  71%  |                                                                                                                                   |========================================================================================                                   |  71%  |                                                                                                                                   |========================================================================================                                   |  72%  |                                                                                                                                   |=========================================================================================                                  |  72%  |                                                                                                                                   |=========================================================================================                                  |  73%  |                                                                                                                                   |==========================================================================================                                 |  73%  |                                                                                                                                   |==========================================================================================                                 |  74%  |                                                                                                                                   |===========================================================================================                                |  74%  |                                                                                                                                   |============================================================================================                               |  74%  |                                                                                                                                   |============================================================================================                               |  75%  |                                                                                                                                   |=============================================================================================                              |  75%  |                                                                                                                                   |=============================================================================================                              |  76%  |                                                                                                                                   |==============================================================================================                             |  76%  |                                                                                                                                   |==============================================================================================                             |  77%  |                                                                                                                                   |===============================================================================================                            |  77%  |                                                                                                                                   |===============================================================================================                            |  78%  |                                                                                                                                   |================================================================================================                           |  78%  |                                                                                                                                   |=================================================================================================                          |  78%  |                                                                                                                                   |=================================================================================================                          |  79%  |                                                                                                                                   |==================================================================================================                         |  79%  |                                                                                                                                   |==================================================================================================                         |  80%  |                                                                                                                                   |===================================================================================================                        |  80%  |                                                                                                                                   |===================================================================================================                        |  81%  |                                                                                                                                   |====================================================================================================                       |  81%  |                                                                                                                                   |====================================================================================================                       |  82%  |                                                                                                                                   |=====================================================================================================                      |  82%  |                                                                                                                                   |=====================================================================================================                      |  83%  |                                                                                                                                   |======================================================================================================                     |  83%  |                                                                                                                                   |=======================================================================================================                    |  83%  |                                                                                                                                   |=======================================================================================================                    |  84%  |                                                                                                                                   |========================================================================================================                   |  84%  |                                                                                                                                   |========================================================================================================                   |  85%  |                                                                                                                                   |=========================================================================================================                  |  85%  |                                                                                                                                   |=========================================================================================================                  |  86%  |                                                                                                                                   |==========================================================================================================                 |  86%  |                                                                                                                                   |==========================================================================================================                 |  87%  |                                                                                                                                   |===========================================================================================================                |  87%  |                                                                                                                                   |============================================================================================================               |  87%  |                                                                                                                                   |============================================================================================================               |  88%  |                                                                                                                                   |=============================================================================================================              |  88%  |                                                                                                                                   |=============================================================================================================              |  89%  |                                                                                                                                   |==============================================================================================================             |  89%  |                                                                                                                                   |==============================================================================================================             |  90%  |                                                                                                                                   |===============================================================================================================            |  90%  |                                                                                                                                   |===============================================================================================================            |  91%  |                                                                                                                                   |================================================================================================================           |  91%  |                                                                                                                                   |=================================================================================================================          |  91%  |                                                                                                                                   |=================================================================================================================          |  92%  |                                                                                                                                   |==================================================================================================================         |  92%  |                                                                                                                                   |==================================================================================================================         |  93%  |                                                                                                                                   |===================================================================================================================        |  93%  |                                                                                                                                   |===================================================================================================================        |  94%  |                                                                                                                                   |====================================================================================================================       |  94%  |                                                                                                                                   |====================================================================================================================       |  95%  |                                                                                                                                   |=====================================================================================================================      |  95%  |                                                                                                                                   |=====================================================================================================================      |  96%  |                                                                                                                                   |======================================================================================================================     |  96%  |                                                                                                                                   |=======================================================================================================================    |  96%  |                                                                                                                                   |=======================================================================================================================    |  97%  |                                                                                                                                   |========================================================================================================================   |  97%  |                                                                                                                                   |========================================================================================================================   |  98%  |                                                                                                                                   |=========================================================================================================================  |  98%  |                                                                                                                                   |=========================================================================================================================  |  99%  |                                                                                                                                   |========================================================================================================================== |  99%  |                                                                                                                                   |========================================================================================================================== | 100%  |                                                                                                                                   |===========================================================================================================================| 100%
```

```
## Removing confirmed isotope entries:
```

```
## Recalculating peaklist for plotting:
```

```
##   |                                                                                                                                   |                                                                                                                           |   0%  |                                                                                                                                   |===============================                                                                                            |  25%  |                                                                                                                                   |==============================================================                                                             |  50%  |                                                                                                                                   |============================================================================================                               |  75%  |                                                                                                                                   |===========================================================================================================================| 100%
```

```
## Calculating adduct similarity:
```

```
##   |                                                                                                                                   |                                                                                                                           |   0%  |                                                                                                                                   |===============================                                                                                            |  25%  |                                                                                                                                   |==============================================================                                                             |  50%  |                                                                                                                                   |============================================================================================                               |  75%  |                                                                                                                                   |===========================================================================================================================| 100%
```

```
## 100 out of 4435
```

```
## 200 out of 4435
```

```
## 300 out of 4435
```

```
## 400 out of 4435
```

```
## 500 out of 4435
```

```
## 600 out of 4435
```

```
## 700 out of 4435
```

```
## 800 out of 4435
```

```
## 900 out of 4435
```

```
## 1000 out of 4435
```

```
## 1100 out of 4435
```

```
## 1200 out of 4435
```

```
## 1300 out of 4435
```

```
## 1400 out of 4435
```

```
## 1500 out of 4435
```

```
## 1600 out of 4435
```

```
## 1700 out of 4435
```

```
## 1800 out of 4435
```

```
## 1900 out of 4435
```

```
## 2000 out of 4435
```

```
## 2100 out of 4435
```

```
## 2200 out of 4435
```

```
## 2300 out of 4435
```

```
## 2400 out of 4435
```

```
## 2500 out of 4435
```

```
## 2600 out of 4435
```

```
## 2700 out of 4435
```

```
## 2800 out of 4435
```

```
## 2900 out of 4435
```

```
## 3000 out of 4435
```

```
## 3100 out of 4435
```

```
## 3200 out of 4435
```

```
## 3300 out of 4435
```

```
## 3400 out of 4435
```

```
## 3500 out of 4435
```

```
## 3600 out of 4435
```

```
## 3700 out of 4435
```

```
## 3800 out of 4435
```

```
## 3900 out of 4435
```

```
## 4000 out of 4435
```

```
## 4100 out of 4435
```

```
## 4200 out of 4435
```

```
## 4300 out of 4435
```

```
## 4400 out of 4435
```

```
##   |                                                                                                                                   |                                                                                                                           |   0%  |                                                                                                                                   |===============================                                                                                            |  25%  |                                                                                                                                   |==============================================================                                                             |  50%  |                                                                                                                                   |============================================================================================                               |  75%  |                                                                                                                                   |===========================================================================================================================| 100%
```

```
## Computing isotopic elution similarity:
```

```
## Calculating isotope similarity:
```

```
##   |                                                                                                                                   |                                                                                                                           |   0%  |                                                                                                                                   |===============================                                                                                            |  25%  |                                                                                                                                   |==============================================================                                                             |  50%  |                                                                                                                                   |============================================================================================                               |  75%  |                                                                                                                                   |===========================================================================================================================| 100%
```

```
## Calculating isotopic fidelity metrics:
```

```
##   |                                                                                                                                   |                                                                                                                           |   0%  |                                                                                                                                   |=                                                                                                                          |   0%  |                                                                                                                                   |=                                                                                                                          |   1%  |                                                                                                                                   |==                                                                                                                         |   1%  |                                                                                                                                   |==                                                                                                                         |   2%  |                                                                                                                                   |===                                                                                                                        |   2%  |                                                                                                                                   |===                                                                                                                        |   3%  |                                                                                                                                   |====                                                                                                                       |   3%  |                                                                                                                                   |====                                                                                                                       |   4%  |                                                                                                                                   |=====                                                                                                                      |   4%  |                                                                                                                                   |======                                                                                                                     |   5%  |                                                                                                                                   |=======                                                                                                                    |   5%  |                                                                                                                                   |=======                                                                                                                    |   6%  |                                                                                                                                   |========                                                                                                                   |   6%  |                                                                                                                                   |========                                                                                                                   |   7%  |                                                                                                                                   |=========                                                                                                                  |   7%  |                                                                                                                                   |=========                                                                                                                  |   8%  |                                                                                                                                   |==========                                                                                                                 |   8%  |                                                                                                                                   |==========                                                                                                                 |   9%  |                                                                                                                                   |===========                                                                                                                |   9%  |                                                                                                                                   |============                                                                                                               |   9%  |                                                                                                                                   |============                                                                                                               |  10%  |                                                                                                                                   |=============                                                                                                              |  10%  |                                                                                                                                   |=============                                                                                                              |  11%  |                                                                                                                                   |==============                                                                                                             |  11%  |                                                                                                                                   |==============                                                                                                             |  12%  |                                                                                                                                   |===============                                                                                                            |  12%  |                                                                                                                                   |===============                                                                                                            |  13%  |                                                                                                                                   |================                                                                                                           |  13%  |                                                                                                                                   |=================                                                                                                          |  13%  |                                                                                                                                   |=================                                                                                                          |  14%  |                                                                                                                                   |==================                                                                                                         |  14%  |                                                                                                                                   |==================                                                                                                         |  15%  |                                                                                                                                   |===================                                                                                                        |  15%  |                                                                                                                                   |===================                                                                                                        |  16%  |                                                                                                                                   |====================                                                                                                       |  16%  |                                                                                                                                   |====================                                                                                                       |  17%  |                                                                                                                                   |=====================                                                                                                      |  17%  |                                                                                                                                   |======================                                                                                                     |  17%  |                                                                                                                                   |======================                                                                                                     |  18%  |                                                                                                                                   |=======================                                                                                                    |  18%  |                                                                                                                                   |=======================                                                                                                    |  19%  |                                                                                                                                   |========================                                                                                                   |  19%  |                                                                                                                                   |========================                                                                                                   |  20%  |                                                                                                                                   |=========================                                                                                                  |  20%  |                                                                                                                                   |=========================                                                                                                  |  21%  |                                                                                                                                   |==========================                                                                                                 |  21%  |                                                                                                                                   |===========================                                                                                                |  22%  |                                                                                                                                   |============================                                                                                               |  22%  |                                                                                                                                   |============================                                                                                               |  23%  |                                                                                                                                   |=============================                                                                                              |  23%  |                                                                                                                                   |=============================                                                                                              |  24%  |                                                                                                                                   |==============================                                                                                             |  24%  |                                                                                                                                   |==============================                                                                                             |  25%  |                                                                                                                                   |===============================                                                                                            |  25%  |                                                                                                                                   |===============================                                                                                            |  26%  |                                                                                                                                   |================================                                                                                           |  26%  |                                                                                                                                   |=================================                                                                                          |  26%  |                                                                                                                                   |=================================                                                                                          |  27%  |                                                                                                                                   |==================================                                                                                         |  27%  |                                                                                                                                   |==================================                                                                                         |  28%  |                                                                                                                                   |===================================                                                                                        |  28%  |                                                                                                                                   |===================================                                                                                        |  29%  |                                                                                                                                   |====================================                                                                                       |  29%  |                                                                                                                                   |====================================                                                                                       |  30%  |                                                                                                                                   |=====================================                                                                                      |  30%  |                                                                                                                                   |======================================                                                                                     |  31%  |                                                                                                                                   |=======================================                                                                                    |  31%  |                                                                                                                                   |=======================================                                                                                    |  32%  |                                                                                                                                   |========================================                                                                                   |  32%  |                                                                                                                                   |========================================                                                                                   |  33%  |                                                                                                                                   |=========================================                                                                                  |  33%  |                                                                                                                                   |=========================================                                                                                  |  34%  |                                                                                                                                   |==========================================                                                                                 |  34%  |                                                                                                                                   |==========================================                                                                                 |  35%  |                                                                                                                                   |===========================================                                                                                |  35%  |                                                                                                                                   |============================================                                                                               |  35%  |                                                                                                                                   |============================================                                                                               |  36%  |                                                                                                                                   |=============================================                                                                              |  36%  |                                                                                                                                   |=============================================                                                                              |  37%  |                                                                                                                                   |==============================================                                                                             |  37%  |                                                                                                                                   |==============================================                                                                             |  38%  |                                                                                                                                   |===============================================                                                                            |  38%  |                                                                                                                                   |===============================================                                                                            |  39%  |                                                                                                                                   |================================================                                                                           |  39%  |                                                                                                                                   |=================================================                                                                          |  39%  |                                                                                                                                   |=================================================                                                                          |  40%  |                                                                                                                                   |==================================================                                                                         |  40%  |                                                                                                                                   |==================================================                                                                         |  41%  |                                                                                                                                   |===================================================                                                                        |  41%  |                                                                                                                                   |===================================================                                                                        |  42%  |                                                                                                                                   |====================================================                                                                       |  42%  |                                                                                                                                   |====================================================                                                                       |  43%  |                                                                                                                                   |=====================================================                                                                      |  43%  |                                                                                                                                   |======================================================                                                                     |  44%  |                                                                                                                                   |=======================================================                                                                    |  44%  |                                                                                                                                   |=======================================================                                                                    |  45%  |                                                                                                                                   |========================================================                                                                   |  45%  |                                                                                                                                   |========================================================                                                                   |  46%  |                                                                                                                                   |=========================================================                                                                  |  46%  |                                                                                                                                   |=========================================================                                                                  |  47%  |                                                                                                                                   |==========================================================                                                                 |  47%  |                                                                                                                                   |==========================================================                                                                 |  48%  |                                                                                                                                   |===========================================================                                                                |  48%  |                                                                                                                                   |============================================================                                                               |  48%  |                                                                                                                                   |============================================================                                                               |  49%  |                                                                                                                                   |=============================================================                                                              |  49%  |                                                                                                                                   |=============================================================                                                              |  50%  |                                                                                                                                   |==============================================================                                                             |  50%  |                                                                                                                                   |==============================================================                                                             |  51%  |                                                                                                                                   |===============================================================                                                            |  51%  |                                                                                                                                   |===============================================================                                                            |  52%  |                                                                                                                                   |================================================================                                                           |  52%  |                                                                                                                                   |=================================================================                                                          |  52%  |                                                                                                                                   |=================================================================                                                          |  53%  |                                                                                                                                   |==================================================================                                                         |  53%  |                                                                                                                                   |==================================================================                                                         |  54%  |                                                                                                                                   |===================================================================                                                        |  54%  |                                                                                                                                   |===================================================================                                                        |  55%  |                                                                                                                                   |====================================================================                                                       |  55%  |                                                                                                                                   |====================================================================                                                       |  56%  |                                                                                                                                   |=====================================================================                                                      |  56%  |                                                                                                                                   |======================================================================                                                     |  57%  |                                                                                                                                   |=======================================================================                                                    |  57%  |                                                                                                                                   |=======================================================================                                                    |  58%  |                                                                                                                                   |========================================================================                                                   |  58%  |                                                                                                                                   |========================================================================                                                   |  59%  |                                                                                                                                   |=========================================================================                                                  |  59%  |                                                                                                                                   |=========================================================================                                                  |  60%  |                                                                                                                                   |==========================================================================                                                 |  60%  |                                                                                                                                   |==========================================================================                                                 |  61%  |                                                                                                                                   |===========================================================================                                                |  61%  |                                                                                                                                   |============================================================================                                               |  61%  |                                                                                                                                   |============================================================================                                               |  62%  |                                                                                                                                   |=============================================================================                                              |  62%  |                                                                                                                                   |=============================================================================                                              |  63%  |                                                                                                                                   |==============================================================================                                             |  63%  |                                                                                                                                   |==============================================================================                                             |  64%  |                                                                                                                                   |===============================================================================                                            |  64%  |                                                                                                                                   |===============================================================================                                            |  65%  |                                                                                                                                   |================================================================================                                           |  65%  |                                                                                                                                   |=================================================================================                                          |  65%  |                                                                                                                                   |=================================================================================                                          |  66%  |                                                                                                                                   |==================================================================================                                         |  66%  |                                                                                                                                   |==================================================================================                                         |  67%  |                                                                                                                                   |===================================================================================                                        |  67%  |                                                                                                                                   |===================================================================================                                        |  68%  |                                                                                                                                   |====================================================================================                                       |  68%  |                                                                                                                                   |====================================================================================                                       |  69%  |                                                                                                                                   |=====================================================================================                                      |  69%  |                                                                                                                                   |======================================================================================                                     |  70%  |                                                                                                                                   |=======================================================================================                                    |  70%  |                                                                                                                                   |=======================================================================================                                    |  71%  |                                                                                                                                   |========================================================================================                                   |  71%  |                                                                                                                                   |========================================================================================                                   |  72%  |                                                                                                                                   |=========================================================================================                                  |  72%  |                                                                                                                                   |=========================================================================================                                  |  73%  |                                                                                                                                   |==========================================================================================                                 |  73%  |                                                                                                                                   |==========================================================================================                                 |  74%  |                                                                                                                                   |===========================================================================================                                |  74%  |                                                                                                                                   |============================================================================================                               |  74%  |                                                                                                                                   |============================================================================================                               |  75%  |                                                                                                                                   |=============================================================================================                              |  75%  |                                                                                                                                   |=============================================================================================                              |  76%  |                                                                                                                                   |==============================================================================================                             |  76%  |                                                                                                                                   |==============================================================================================                             |  77%  |                                                                                                                                   |===============================================================================================                            |  77%  |                                                                                                                                   |===============================================================================================                            |  78%  |                                                                                                                                   |================================================================================================                           |  78%  |                                                                                                                                   |=================================================================================================                          |  79%  |                                                                                                                                   |==================================================================================================                         |  79%  |                                                                                                                                   |==================================================================================================                         |  80%  |                                                                                                                                   |===================================================================================================                        |  80%  |                                                                                                                                   |===================================================================================================                        |  81%  |                                                                                                                                   |====================================================================================================                       |  81%  |                                                                                                                                   |====================================================================================================                       |  82%  |                                                                                                                                   |=====================================================================================================                      |  82%  |                                                                                                                                   |=====================================================================================================                      |  83%  |                                                                                                                                   |======================================================================================================                     |  83%  |                                                                                                                                   |=======================================================================================================                    |  83%  |                                                                                                                                   |=======================================================================================================                    |  84%  |                                                                                                                                   |========================================================================================================                   |  84%  |                                                                                                                                   |========================================================================================================                   |  85%  |                                                                                                                                   |=========================================================================================================                  |  85%  |                                                                                                                                   |=========================================================================================================                  |  86%  |                                                                                                                                   |==========================================================================================================                 |  86%  |                                                                                                                                   |==========================================================================================================                 |  87%  |                                                                                                                                   |===========================================================================================================                |  87%  |                                                                                                                                   |============================================================================================================               |  87%  |                                                                                                                                   |============================================================================================================               |  88%  |                                                                                                                                   |=============================================================================================================              |  88%  |                                                                                                                                   |=============================================================================================================              |  89%  |                                                                                                                                   |==============================================================================================================             |  89%  |                                                                                                                                   |==============================================================================================================             |  90%  |                                                                                                                                   |===============================================================================================================            |  90%  |                                                                                                                                   |===============================================================================================================            |  91%  |                                                                                                                                   |================================================================================================================           |  91%  |                                                                                                                                   |=================================================================================================================          |  91%  |                                                                                                                                   |=================================================================================================================          |  92%  |                                                                                                                                   |==================================================================================================================         |  92%  |                                                                                                                                   |==================================================================================================================         |  93%  |                                                                                                                                   |===================================================================================================================        |  93%  |                                                                                                                                   |===================================================================================================================        |  94%  |                                                                                                                                   |====================================================================================================================       |  94%  |                                                                                                                                   |====================================================================================================================       |  95%  |                                                                                                                                   |=====================================================================================================================      |  95%  |                                                                                                                                   |======================================================================================================================     |  96%  |                                                                                                                                   |=======================================================================================================================    |  96%  |                                                                                                                                   |=======================================================================================================================    |  97%  |                                                                                                                                   |========================================================================================================================   |  97%  |                                                                                                                                   |========================================================================================================================   |  98%  |                                                                                                                                   |=========================================================================================================================  |  98%  |                                                                                                                                   |=========================================================================================================================  |  99%  |                                                                                                                                   |========================================================================================================================== |  99%  |                                                                                                                                   |========================================================================================================================== | 100%  |                                                                                                                                   |===========================================================================================================================| 100%
```

```
## Removing confirmed isotope entries:
```

```
## Recalculating peaklist for plotting:
```

```
##   |                                                                                                                                   |                                                                                                                           |   0%  |                                                                                                                                   |===============================                                                                                            |  25%  |                                                                                                                                   |==============================================================                                                             |  50%  |                                                                                                                                   |============================================================================================                               |  75%  |                                                                                                                                   |===========================================================================================================================| 100%
```

```
## Calculating adduct similarity:
```

```
##   |                                                                                                                                   |                                                                                                                           |   0%  |                                                                                                                                   |===============================                                                                                            |  25%  |                                                                                                                                   |==============================================================                                                             |  50%  |                                                                                                                                   |============================================================================================                               |  75%  |                                                                                                                                   |===========================================================================================================================| 100%
```

```
## 100 out of 1157
```

```
## 200 out of 1157
```

```
## 300 out of 1157
```

```
## 400 out of 1157
```

```
## 500 out of 1157
```

```
## 600 out of 1157
```

```
## 700 out of 1157
```

```
## 800 out of 1157
```

```
## 900 out of 1157
```

```
## 1000 out of 1157
```

```
## 1100 out of 1157
```

```
##   |                                                                                                                                   |                                                                                                                           |   0%  |                                                                                                                                   |===============================                                                                                            |  25%  |                                                                                                                                   |==============================================================                                                             |  50%  |                                                                                                                                   |============================================================================================                               |  75%  |                                                                                                                                   |===========================================================================================================================| 100%
```

```r
saveRDS(myHermes, file = "D:/HermesResults/SurfaceWater/PosSOIFilter.rds")
```


```r
s <- getSOIpar("double")
negHermes <- SOIfinder(negHermes, params = s, fileID = c(1,2,2,3,3,3),
                      against = c(0,0,1,0,1,2))
```

```
## Starting density filtering:
```

```
## Running Density Filter
```

```
##   |                                                                                                                                   |                                                                                                                           |   0%  |                                                                                                                                   |===============================                                                                                            |  25%  |                                                                                                                                   |==============================================================                                                             |  50%  |                                                                                                                                   |============================================================================================                               |  75%  |                                                                                                                                   |===========================================================================================================================| 100%
```

```
## Running Density Interpreter
```

```
##   |                                                                                                                                   |                                                                                                                           |   0%  |                                                                                                                                   |===============================                                                                                            |  25%  |                                                                                                                                   |==============================================================                                                             |  50%  |                                                                                                                                   |============================================================================================                               |  75%  |                                                                                                                                   |===========================================================================================================================| 100%
```

```
## Running Density Filter
```

```
##   |                                                                                                                                   |                                                                                                                           |   0%  |                                                                                                                                   |===============================                                                                                            |  25%  |                                                                                                                                   |==============================================================                                                             |  50%  |                                                                                                                                   |============================================================================================                               |  75%  |                                                                                                                                   |===========================================================================================================================| 100%
```

```
## Running Density Interpreter
```

```
##   |                                                                                                                                   |                                                                                                                           |   0%  |                                                                                                                                   |===============================                                                                                            |  25%  |                                                                                                                                   |==============================================================                                                             |  50%  |                                                                                                                                   |============================================================================================                               |  75%  |                                                                                                                                   |===========================================================================================================================| 100%
```

```
## Now Grouping:
```

```
## Merging filter 2 out of 2 with previous SOI list
```

```
##   |                                                                                                                                   |                                                                                                                           |   0%  |                                                                                                                                   |===============================                                                                                            |  25%  |                                                                                                                                   |==============================================================                                                             |  50%  |                                                                                                                                   |============================================================================================                               |  75%  |                                                                                                                                   |===========================================================================================================================| 100%
```

```
## Adding entries unique from the first DF
```

```
## Adding entries unique from the second DF
```

```
## Solving inner conflicts, Round: 1
```

```
## Solving inner conflicts, Round: 2
```

```
## Solving inner conflicts, Round: 3
```

```
## Initial peak retrieval:
```

```
##   |                                                                                                                                   |                                                                                                                           |   0%  |                                                                                                                                   |===============================                                                                                            |  25%  |                                                                                                                                   |==============================================================                                                             |  50%  |                                                                                                                                   |============================================================================================                               |  75%  |                                                                                                                                   |===========================================================================================================================| 100%
```

```
## 
```

```
## Starting group characterization:
```

```
## Mass calculation:
```

```
## Shortening and selecting long groups:
```

```
##   |                                                                                                                                   |                                                                                                                           |   0%  |                                                                                                                                   |===============================                                                                                            |  25%  |                                                                                                                                   |==============================================================                                                             |  50%  |                                                                                                                                   |============================================================================================                               |  75%  |                                                                                                                                   |===========================================================================================================================| 100%
```

```
## 
```

```
## Width calculation:
```

```
## Maximum intensity calculation:
```

```
## Number of scans:
```

```
## Converting from ionic formula to F/A combinations:
```

```
## Generating peaklist for plotting:
```

```
##   |                                                                                                                                   |                                                                                                                           |   0%  |                                                                                                                                   |===============================                                                                                            |  25%  |                                                                                                                                   |==============================================================                                                             |  50%  |                                                                                                                                   |============================================================================================                               |  75%  |                                                                                                                                   |===========================================================================================================================| 100%
```

```
## Starting density filtering:
```

```
## Running Density Filter
```

```
##   |                                                                                                                                   |                                                                                                                           |   0%  |                                                                                                                                   |===============================                                                                                            |  25%  |                                                                                                                                   |==============================================================                                                             |  50%  |                                                                                                                                   |============================================================================================                               |  75%  |                                                                                                                                   |===========================================================================================================================| 100%
```

```
## Running Density Interpreter
```

```
##   |                                                                                                                                   |                                                                                                                           |   0%  |                                                                                                                                   |===============================                                                                                            |  25%  |                                                                                                                                   |==============================================================                                                             |  50%  |                                                                                                                                   |============================================================================================                               |  75%  |                                                                                                                                   |===========================================================================================================================| 100%
```

```
## Running Density Filter
```

```
##   |                                                                                                                                   |                                                                                                                           |   0%  |                                                                                                                                   |===============================                                                                                            |  25%  |                                                                                                                                   |==============================================================                                                             |  50%  |                                                                                                                                   |============================================================================================                               |  75%  |                                                                                                                                   |===========================================================================================================================| 100%
```

```
## Running Density Interpreter
```

```
##   |                                                                                                                                   |                                                                                                                           |   0%  |                                                                                                                                   |===============================                                                                                            |  25%  |                                                                                                                                   |==============================================================                                                             |  50%  |                                                                                                                                   |============================================================================================                               |  75%  |                                                                                                                                   |===========================================================================================================================| 100%
```

```
## Now Grouping:
```

```
## Merging filter 2 out of 2 with previous SOI list
```

```
##   |                                                                                                                                   |                                                                                                                           |   0%  |                                                                                                                                   |===============================                                                                                            |  25%  |                                                                                                                                   |==============================================================                                                             |  50%  |                                                                                                                                   |============================================================================================                               |  75%  |                                                                                                                                   |===========================================================================================================================| 100%
```

```
## Adding entries unique from the first DF
```

```
## Adding entries unique from the second DF
```

```
## Solving inner conflicts, Round: 1
```

```
## Solving inner conflicts, Round: 2
```

```
## Solving inner conflicts, Round: 3
```

```
## Initial peak retrieval:
```

```
##   |                                                                                                                                   |                                                                                                                           |   0%  |                                                                                                                                   |===============================                                                                                            |  25%  |                                                                                                                                   |==============================================================                                                             |  50%  |                                                                                                                                   |============================================================================================                               |  75%  |                                                                                                                                   |===========================================================================================================================| 100%
```

```
## 
```

```
## Starting group characterization:
```

```
## Mass calculation:
```

```
## Shortening and selecting long groups:
```

```
##   |                                                                                                                                   |                                                                                                                           |   0%  |                                                                                                                                   |===============================                                                                                            |  25%  |                                                                                                                                   |==============================================================                                                             |  50%  |                                                                                                                                   |============================================================================================                               |  75%  |                                                                                                                                   |===========================================================================================================================| 100%
```

```
## 
```

```
## Width calculation:
```

```
## Maximum intensity calculation:
```

```
## Number of scans:
```

```
## Converting from ionic formula to F/A combinations:
```

```
## Generating peaklist for plotting:
```

```
##   |                                                                                                                                   |                                                                                                                           |   0%  |                                                                                                                                   |===============================                                                                                            |  25%  |                                                                                                                                   |==============================================================                                                             |  50%  |                                                                                                                                   |============================================================================================                               |  75%  |                                                                                                                                   |===========================================================================================================================| 100%
```

```
## Starting density filtering:
```

```
## Running Density Filter
```

```
##   |                                                                                                                                   |                                                                                                                           |   0%  |                                                                                                                                   |===============================                                                                                            |  25%  |                                                                                                                                   |==============================================================                                                             |  50%  |                                                                                                                                   |============================================================================================                               |  75%  |                                                                                                                                   |===========================================================================================================================| 100%
```

```
## Running Density Interpreter
```

```
##   |                                                                                                                                   |                                                                                                                           |   0%  |                                                                                                                                   |===============================                                                                                            |  25%  |                                                                                                                                   |==============================================================                                                             |  50%  |                                                                                                                                   |============================================================================================                               |  75%  |                                                                                                                                   |===========================================================================================================================| 100%
```

```
## Running Density Filter
```

```
##   |                                                                                                                                   |                                                                                                                           |   0%  |                                                                                                                                   |===============================                                                                                            |  25%  |                                                                                                                                   |==============================================================                                                             |  50%  |                                                                                                                                   |============================================================================================                               |  75%  |                                                                                                                                   |===========================================================================================================================| 100%
```

```
## Running Density Interpreter
```

```
##   |                                                                                                                                   |                                                                                                                           |   0%  |                                                                                                                                   |===============================                                                                                            |  25%  |                                                                                                                                   |==============================================================                                                             |  50%  |                                                                                                                                   |============================================================================================                               |  75%  |                                                                                                                                   |===========================================================================================================================| 100%
```

```
## Now Grouping:
```

```
## Merging filter 2 out of 2 with previous SOI list
```

```
##   |                                                                                                                                   |                                                                                                                           |   0%  |                                                                                                                                   |===============================                                                                                            |  25%  |                                                                                                                                   |==============================================================                                                             |  50%  |                                                                                                                                   |============================================================================================                               |  75%  |                                                                                                                                   |===========================================================================================================================| 100%
```

```
## Adding entries unique from the first DF
```

```
## Adding entries unique from the second DF
```

```
## Solving inner conflicts, Round: 1
```

```
## Solving inner conflicts, Round: 2
```

```
## Solving inner conflicts, Round: 3
```

```
## Initial peak retrieval:
```

```
##   |                                                                                                                                   |                                                                                                                           |   0%  |                                                                                                                                   |===============================                                                                                            |  25%  |                                                                                                                                   |==============================================================                                                             |  50%  |                                                                                                                                   |============================================================================================                               |  75%  |                                                                                                                                   |===========================================================================================================================| 100%
```

```
## 
```

```
## Starting group characterization:
```

```
## Mass calculation:
```

```
## Shortening and selecting long groups:
```

```
##   |                                                                                                                                   |                                                                                                                           |   0%  |                                                                                                                                   |===============================                                                                                            |  25%  |                                                                                                                                   |==============================================================                                                             |  50%  |                                                                                                                                   |============================================================================================                               |  75%  |                                                                                                                                   |===========================================================================================================================| 100%
```

```
## Blank substraction:
```

```
##   |                                                                                                                                   |                                                                                                                           |   0%  |                                                                                                                                   |===============================                                                                                            |  25%  |                                                                                                                                   |==============================================================                                                             |  50%  |                                                                                                                                   |============================================================================================                               |  75%  |                                                                                                                                   |===========================================================================================================================| 100%
## 
##   |                                                                                                                                   |                                                                                                                           |   0%  |                                                                                                                                   |===============================                                                                                            |  25%  |                                                                                                                                   |==============================================================                                                             |  50%  |                                                                                                                                   |============================================================================================                               |  75%  |                                                                                                                                   |===========================================================================================================================| 100%
```

```
## 
```

```
## Width calculation:
```

```
## Maximum intensity calculation:
```

```
## Number of scans:
```

```
## Converting from ionic formula to F/A combinations:
```

```
## Generating peaklist for plotting:
```

```
##   |                                                                                                                                   |                                                                                                                           |   0%  |                                                                                                                                   |===============================                                                                                            |  25%  |                                                                                                                                   |==============================================================                                                             |  50%  |                                                                                                                                   |============================================================================================                               |  75%  |                                                                                                                                   |===========================================================================================================================| 100%
```

```
## Starting density filtering:
```

```
## Running Density Filter
```

```
##   |                                                                                                                                   |                                                                                                                           |   0%  |                                                                                                                                   |===============================                                                                                            |  25%  |                                                                                                                                   |==============================================================                                                             |  50%  |                                                                                                                                   |============================================================================================                               |  75%  |                                                                                                                                   |===========================================================================================================================| 100%
```

```
## Running Density Interpreter
```

```
##   |                                                                                                                                   |                                                                                                                           |   0%  |                                                                                                                                   |===============================                                                                                            |  25%  |                                                                                                                                   |==============================================================                                                             |  50%  |                                                                                                                                   |============================================================================================                               |  75%  |                                                                                                                                   |===========================================================================================================================| 100%
```

```
## Running Density Filter
```

```
##   |                                                                                                                                   |                                                                                                                           |   0%  |                                                                                                                                   |===============================                                                                                            |  25%  |                                                                                                                                   |==============================================================                                                             |  50%  |                                                                                                                                   |============================================================================================                               |  75%  |                                                                                                                                   |===========================================================================================================================| 100%
```

```
## Running Density Interpreter
```

```
##   |                                                                                                                                   |                                                                                                                           |   0%  |                                                                                                                                   |===============================                                                                                            |  25%  |                                                                                                                                   |==============================================================                                                             |  50%  |                                                                                                                                   |============================================================================================                               |  75%  |                                                                                                                                   |===========================================================================================================================| 100%
```

```
## Now Grouping:
```

```
## Merging filter 2 out of 2 with previous SOI list
```

```
##   |                                                                                                                                   |                                                                                                                           |   0%  |                                                                                                                                   |===============================                                                                                            |  25%  |                                                                                                                                   |==============================================================                                                             |  50%  |                                                                                                                                   |============================================================================================                               |  75%  |                                                                                                                                   |===========================================================================================================================| 100%
```

```
## Adding entries unique from the first DF
```

```
## Adding entries unique from the second DF
```

```
## Solving inner conflicts, Round: 1
```

```
## Solving inner conflicts, Round: 2
```

```
## Solving inner conflicts, Round: 3
```

```
## Initial peak retrieval:
```

```
##   |                                                                                                                                   |                                                                                                                           |   0%  |                                                                                                                                   |===============================                                                                                            |  25%  |                                                                                                                                   |==============================================================                                                             |  50%  |                                                                                                                                   |============================================================================================                               |  75%  |                                                                                                                                   |===========================================================================================================================| 100%
```

```
## 
```

```
## Starting group characterization:
```

```
## Mass calculation:
```

```
## Shortening and selecting long groups:
```

```
##   |                                                                                                                                   |                                                                                                                           |   0%  |                                                                                                                                   |===============================                                                                                            |  25%  |                                                                                                                                   |==============================================================                                                             |  50%  |                                                                                                                                   |============================================================================================                               |  75%  |                                                                                                                                   |===========================================================================================================================| 100%
```

```
## 
```

```
## Width calculation:
```

```
## Maximum intensity calculation:
```

```
## Number of scans:
```

```
## Converting from ionic formula to F/A combinations:
```

```
## Generating peaklist for plotting:
```

```
##   |                                                                                                                                   |                                                                                                                           |   0%  |                                                                                                                                   |===============================                                                                                            |  25%  |                                                                                                                                   |==============================================================                                                             |  50%  |                                                                                                                                   |============================================================================================                               |  75%  |                                                                                                                                   |===========================================================================================================================| 100%
```

```
## Starting density filtering:
```

```
## Running Density Filter
```

```
##   |                                                                                                                                   |                                                                                                                           |   0%  |                                                                                                                                   |===============================                                                                                            |  25%  |                                                                                                                                   |==============================================================                                                             |  50%  |                                                                                                                                   |============================================================================================                               |  75%  |                                                                                                                                   |===========================================================================================================================| 100%
```

```
## Running Density Interpreter
```

```
##   |                                                                                                                                   |                                                                                                                           |   0%  |                                                                                                                                   |===============================                                                                                            |  25%  |                                                                                                                                   |==============================================================                                                             |  50%  |                                                                                                                                   |============================================================================================                               |  75%  |                                                                                                                                   |===========================================================================================================================| 100%
```

```
## Running Density Filter
```

```
##   |                                                                                                                                   |                                                                                                                           |   0%  |                                                                                                                                   |===============================                                                                                            |  25%  |                                                                                                                                   |==============================================================                                                             |  50%  |                                                                                                                                   |============================================================================================                               |  75%  |                                                                                                                                   |===========================================================================================================================| 100%
```

```
## Running Density Interpreter
```

```
##   |                                                                                                                                   |                                                                                                                           |   0%  |                                                                                                                                   |===============================                                                                                            |  25%  |                                                                                                                                   |==============================================================                                                             |  50%  |                                                                                                                                   |============================================================================================                               |  75%  |                                                                                                                                   |===========================================================================================================================| 100%
```

```
## Now Grouping:
```

```
## Merging filter 2 out of 2 with previous SOI list
```

```
##   |                                                                                                                                   |                                                                                                                           |   0%  |                                                                                                                                   |===============================                                                                                            |  25%  |                                                                                                                                   |==============================================================                                                             |  50%  |                                                                                                                                   |============================================================================================                               |  75%  |                                                                                                                                   |===========================================================================================================================| 100%
```

```
## Adding entries unique from the first DF
```

```
## Adding entries unique from the second DF
```

```
## Solving inner conflicts, Round: 1
```

```
## Solving inner conflicts, Round: 2
```

```
## Solving inner conflicts, Round: 3
```

```
## Initial peak retrieval:
```

```
##   |                                                                                                                                   |                                                                                                                           |   0%  |                                                                                                                                   |===============================                                                                                            |  25%  |                                                                                                                                   |==============================================================                                                             |  50%  |                                                                                                                                   |============================================================================================                               |  75%  |                                                                                                                                   |===========================================================================================================================| 100%
```

```
## 
```

```
## Starting group characterization:
```

```
## Mass calculation:
```

```
## Shortening and selecting long groups:
```

```
##   |                                                                                                                                   |                                                                                                                           |   0%  |                                                                                                                                   |===============================                                                                                            |  25%  |                                                                                                                                   |==============================================================                                                             |  50%  |                                                                                                                                   |============================================================================================                               |  75%  |                                                                                                                                   |===========================================================================================================================| 100%
```

```
## Blank substraction:
```

```
##   |                                                                                                                                   |                                                                                                                           |   0%  |                                                                                                                                   |===============================                                                                                            |  25%  |                                                                                                                                   |==============================================================                                                             |  50%  |                                                                                                                                   |============================================================================================                               |  75%  |                                                                                                                                   |===========================================================================================================================| 100%
## 
##   |                                                                                                                                   |                                                                                                                           |   0%  |                                                                                                                                   |===============================                                                                                            |  25%  |                                                                                                                                   |==============================================================                                                             |  50%  |                                                                                                                                   |============================================================================================                               |  75%  |                                                                                                                                   |===========================================================================================================================| 100%
```

```
## 
```

```
## Width calculation:
```

```
## Maximum intensity calculation:
```

```
## Number of scans:
```

```
## Converting from ionic formula to F/A combinations:
```

```
## Generating peaklist for plotting:
```

```
##   |                                                                                                                                   |                                                                                                                           |   0%  |                                                                                                                                   |===============================                                                                                            |  25%  |                                                                                                                                   |==============================================================                                                             |  50%  |                                                                                                                                   |============================================================================================                               |  75%  |                                                                                                                                   |===========================================================================================================================| 100%
```

```
## Starting density filtering:
```

```
## Running Density Filter
```

```
##   |                                                                                                                                   |                                                                                                                           |   0%  |                                                                                                                                   |===============================                                                                                            |  25%  |                                                                                                                                   |==============================================================                                                             |  50%  |                                                                                                                                   |============================================================================================                               |  75%  |                                                                                                                                   |===========================================================================================================================| 100%
```

```
## Running Density Interpreter
```

```
##   |                                                                                                                                   |                                                                                                                           |   0%  |                                                                                                                                   |===============================                                                                                            |  25%  |                                                                                                                                   |==============================================================                                                             |  50%  |                                                                                                                                   |============================================================================================                               |  75%  |                                                                                                                                   |===========================================================================================================================| 100%
```

```
## Running Density Filter
```

```
##   |                                                                                                                                   |                                                                                                                           |   0%  |                                                                                                                                   |===============================                                                                                            |  25%  |                                                                                                                                   |==============================================================                                                             |  50%  |                                                                                                                                   |============================================================================================                               |  75%  |                                                                                                                                   |===========================================================================================================================| 100%
```

```
## Running Density Interpreter
```

```
##   |                                                                                                                                   |                                                                                                                           |   0%  |                                                                                                                                   |===============================                                                                                            |  25%  |                                                                                                                                   |==============================================================                                                             |  50%  |                                                                                                                                   |============================================================================================                               |  75%  |                                                                                                                                   |===========================================================================================================================| 100%
```

```
## Now Grouping:
```

```
## Merging filter 2 out of 2 with previous SOI list
```

```
##   |                                                                                                                                   |                                                                                                                           |   0%  |                                                                                                                                   |===============================                                                                                            |  25%  |                                                                                                                                   |==============================================================                                                             |  50%  |                                                                                                                                   |============================================================================================                               |  75%  |                                                                                                                                   |===========================================================================================================================| 100%
```

```
## Adding entries unique from the first DF
```

```
## Adding entries unique from the second DF
```

```
## Solving inner conflicts, Round: 1
```

```
## Solving inner conflicts, Round: 2
```

```
## Solving inner conflicts, Round: 3
```

```
## Initial peak retrieval:
```

```
##   |                                                                                                                                   |                                                                                                                           |   0%  |                                                                                                                                   |===============================                                                                                            |  25%  |                                                                                                                                   |==============================================================                                                             |  50%  |                                                                                                                                   |============================================================================================                               |  75%  |                                                                                                                                   |===========================================================================================================================| 100%
```

```
## 
```

```
## Starting group characterization:
```

```
## Mass calculation:
```

```
## Shortening and selecting long groups:
```

```
##   |                                                                                                                                   |                                                                                                                           |   0%  |                                                                                                                                   |===============================                                                                                            |  25%  |                                                                                                                                   |==============================================================                                                             |  50%  |                                                                                                                                   |============================================================================================                               |  75%  |                                                                                                                                   |===========================================================================================================================| 100%
```

```
## Blank substraction:
```

```
##   |                                                                                                                                   |                                                                                                                           |   0%  |                                                                                                                                   |===============================                                                                                            |  25%  |                                                                                                                                   |==============================================================                                                             |  50%  |                                                                                                                                   |============================================================================================                               |  75%  |                                                                                                                                   |===========================================================================================================================| 100%
## 
##   |                                                                                                                                   |                                                                                                                           |   0%  |                                                                                                                                   |===============================                                                                                            |  25%  |                                                                                                                                   |==============================================================                                                             |  50%  |                                                                                                                                   |============================================================================================                               |  75%  |                                                                                                                                   |===========================================================================================================================| 100%
```

```
## 
```

```
## Width calculation:
```

```
## Maximum intensity calculation:
```

```
## Number of scans:
```

```
## Converting from ionic formula to F/A combinations:
```

```
## Generating peaklist for plotting:
```

```
##   |                                                                                                                                   |                                                                                                                           |   0%  |                                                                                                                                   |===============================                                                                                            |  25%  |                                                                                                                                   |==============================================================                                                             |  50%  |                                                                                                                                   |============================================================================================                               |  75%  |                                                                                                                                   |===========================================================================================================================| 100%
```

```r
saveRDS(negHermes, file = "D:/HermesResults/SurfaceWater/NegSOI.rds")

for(i in 1:6){
  negHermes <- SOIcleaner(negHermes, soiid = i, minint = 10000, TRUE)
  negHermes <- ISFproc(negHermes, id = i)
}
```

```
## Computing isotopic elution similarity:
```

```
## Calculating isotope similarity:
```

```
##   |                                                                                                                                   |                                                                                                                           |   0%  |                                                                                                                                   |===============================                                                                                            |  25%  |                                                                                                                                   |==============================================================                                                             |  50%  |                                                                                                                                   |============================================================================================                               |  75%  |                                                                                                                                   |===========================================================================================================================| 100%
```

```
## Calculating isotopic fidelity metrics:
```

```
##   |                                                                                                                                   |                                                                                                                           |   0%  |                                                                                                                                   |=                                                                                                                          |   0%  |                                                                                                                                   |=                                                                                                                          |   1%  |                                                                                                                                   |==                                                                                                                         |   1%  |                                                                                                                                   |==                                                                                                                         |   2%  |                                                                                                                                   |===                                                                                                                        |   2%  |                                                                                                                                   |===                                                                                                                        |   3%  |                                                                                                                                   |====                                                                                                                       |   3%  |                                                                                                                                   |====                                                                                                                       |   4%  |                                                                                                                                   |=====                                                                                                                      |   4%  |                                                                                                                                   |======                                                                                                                     |   4%  |                                                                                                                                   |======                                                                                                                     |   5%  |                                                                                                                                   |=======                                                                                                                    |   5%  |                                                                                                                                   |=======                                                                                                                    |   6%  |                                                                                                                                   |========                                                                                                                   |   6%  |                                                                                                                                   |========                                                                                                                   |   7%  |                                                                                                                                   |=========                                                                                                                  |   7%  |                                                                                                                                   |=========                                                                                                                  |   8%  |                                                                                                                                   |==========                                                                                                                 |   8%  |                                                                                                                                   |===========                                                                                                                |   9%  |                                                                                                                                   |============                                                                                                               |   9%  |                                                                                                                                   |============                                                                                                               |  10%  |                                                                                                                                   |=============                                                                                                              |  10%  |                                                                                                                                   |=============                                                                                                              |  11%  |                                                                                                                                   |==============                                                                                                             |  11%  |                                                                                                                                   |==============                                                                                                             |  12%  |                                                                                                                                   |===============                                                                                                            |  12%  |                                                                                                                                   |===============                                                                                                            |  13%  |                                                                                                                                   |================                                                                                                           |  13%  |                                                                                                                                   |=================                                                                                                          |  13%  |                                                                                                                                   |=================                                                                                                          |  14%  |                                                                                                                                   |==================                                                                                                         |  14%  |                                                                                                                                   |==================                                                                                                         |  15%  |                                                                                                                                   |===================                                                                                                        |  15%  |                                                                                                                                   |===================                                                                                                        |  16%  |                                                                                                                                   |====================                                                                                                       |  16%  |                                                                                                                                   |====================                                                                                                       |  17%  |                                                                                                                                   |=====================                                                                                                      |  17%  |                                                                                                                                   |======================                                                                                                     |  18%  |                                                                                                                                   |=======================                                                                                                    |  18%  |                                                                                                                                   |=======================                                                                                                    |  19%  |                                                                                                                                   |========================                                                                                                   |  19%  |                                                                                                                                   |========================                                                                                                   |  20%  |                                                                                                                                   |=========================                                                                                                  |  20%  |                                                                                                                                   |=========================                                                                                                  |  21%  |                                                                                                                                   |==========================                                                                                                 |  21%  |                                                                                                                                   |==========================                                                                                                 |  22%  |                                                                                                                                   |===========================                                                                                                |  22%  |                                                                                                                                   |============================                                                                                               |  22%  |                                                                                                                                   |============================                                                                                               |  23%  |                                                                                                                                   |=============================                                                                                              |  23%  |                                                                                                                                   |=============================                                                                                              |  24%  |                                                                                                                                   |==============================                                                                                             |  24%  |                                                                                                                                   |==============================                                                                                             |  25%  |                                                                                                                                   |===============================                                                                                            |  25%  |                                                                                                                                   |===============================                                                                                            |  26%  |                                                                                                                                   |================================                                                                                           |  26%  |                                                                                                                                   |=================================                                                                                          |  26%  |                                                                                                                                   |=================================                                                                                          |  27%  |                                                                                                                                   |==================================                                                                                         |  27%  |                                                                                                                                   |==================================                                                                                         |  28%  |                                                                                                                                   |===================================                                                                                        |  28%  |                                                                                                                                   |===================================                                                                                        |  29%  |                                                                                                                                   |====================================                                                                                       |  29%  |                                                                                                                                   |====================================                                                                                       |  30%  |                                                                                                                                   |=====================================                                                                                      |  30%  |                                                                                                                                   |======================================                                                                                     |  31%  |                                                                                                                                   |=======================================                                                                                    |  31%  |                                                                                                                                   |=======================================                                                                                    |  32%  |                                                                                                                                   |========================================                                                                                   |  32%  |                                                                                                                                   |========================================                                                                                   |  33%  |                                                                                                                                   |=========================================                                                                                  |  33%  |                                                                                                                                   |=========================================                                                                                  |  34%  |                                                                                                                                   |==========================================                                                                                 |  34%  |                                                                                                                                   |==========================================                                                                                 |  35%  |                                                                                                                                   |===========================================                                                                                |  35%  |                                                                                                                                   |============================================                                                                               |  35%  |                                                                                                                                   |============================================                                                                               |  36%  |                                                                                                                                   |=============================================                                                                              |  36%  |                                                                                                                                   |=============================================                                                                              |  37%  |                                                                                                                                   |==============================================                                                                             |  37%  |                                                                                                                                   |==============================================                                                                             |  38%  |                                                                                                                                   |===============================================                                                                            |  38%  |                                                                                                                                   |===============================================                                                                            |  39%  |                                                                                                                                   |================================================                                                                           |  39%  |                                                                                                                                   |=================================================                                                                          |  39%  |                                                                                                                                   |=================================================                                                                          |  40%  |                                                                                                                                   |==================================================                                                                         |  40%  |                                                                                                                                   |==================================================                                                                         |  41%  |                                                                                                                                   |===================================================                                                                        |  41%  |                                                                                                                                   |===================================================                                                                        |  42%  |                                                                                                                                   |====================================================                                                                       |  42%  |                                                                                                                                   |====================================================                                                                       |  43%  |                                                                                                                                   |=====================================================                                                                      |  43%  |                                                                                                                                   |======================================================                                                                     |  44%  |                                                                                                                                   |=======================================================                                                                    |  44%  |                                                                                                                                   |=======================================================                                                                    |  45%  |                                                                                                                                   |========================================================                                                                   |  45%  |                                                                                                                                   |========================================================                                                                   |  46%  |                                                                                                                                   |=========================================================                                                                  |  46%  |                                                                                                                                   |=========================================================                                                                  |  47%  |                                                                                                                                   |==========================================================                                                                 |  47%  |                                                                                                                                   |==========================================================                                                                 |  48%  |                                                                                                                                   |===========================================================                                                                |  48%  |                                                                                                                                   |============================================================                                                               |  48%  |                                                                                                                                   |============================================================                                                               |  49%  |                                                                                                                                   |=============================================================                                                              |  49%  |                                                                                                                                   |=============================================================                                                              |  50%  |                                                                                                                                   |==============================================================                                                             |  50%  |                                                                                                                                   |==============================================================                                                             |  51%  |                                                                                                                                   |===============================================================                                                            |  51%  |                                                                                                                                   |===============================================================                                                            |  52%  |                                                                                                                                   |================================================================                                                           |  52%  |                                                                                                                                   |=================================================================                                                          |  52%  |                                                                                                                                   |=================================================================                                                          |  53%  |                                                                                                                                   |==================================================================                                                         |  53%  |                                                                                                                                   |==================================================================                                                         |  54%  |                                                                                                                                   |===================================================================                                                        |  54%  |                                                                                                                                   |===================================================================                                                        |  55%  |                                                                                                                                   |====================================================================                                                       |  55%  |                                                                                                                                   |====================================================================                                                       |  56%  |                                                                                                                                   |=====================================================================                                                      |  56%  |                                                                                                                                   |======================================================================                                                     |  57%  |                                                                                                                                   |=======================================================================                                                    |  57%  |                                                                                                                                   |=======================================================================                                                    |  58%  |                                                                                                                                   |========================================================================                                                   |  58%  |                                                                                                                                   |========================================================================                                                   |  59%  |                                                                                                                                   |=========================================================================                                                  |  59%  |                                                                                                                                   |=========================================================================                                                  |  60%  |                                                                                                                                   |==========================================================================                                                 |  60%  |                                                                                                                                   |==========================================================================                                                 |  61%  |                                                                                                                                   |===========================================================================                                                |  61%  |                                                                                                                                   |============================================================================                                               |  61%  |                                                                                                                                   |============================================================================                                               |  62%  |                                                                                                                                   |=============================================================================                                              |  62%  |                                                                                                                                   |=============================================================================                                              |  63%  |                                                                                                                                   |==============================================================================                                             |  63%  |                                                                                                                                   |==============================================================================                                             |  64%  |                                                                                                                                   |===============================================================================                                            |  64%  |                                                                                                                                   |===============================================================================                                            |  65%  |                                                                                                                                   |================================================================================                                           |  65%  |                                                                                                                                   |=================================================================================                                          |  65%  |                                                                                                                                   |=================================================================================                                          |  66%  |                                                                                                                                   |==================================================================================                                         |  66%  |                                                                                                                                   |==================================================================================                                         |  67%  |                                                                                                                                   |===================================================================================                                        |  67%  |                                                                                                                                   |===================================================================================                                        |  68%  |                                                                                                                                   |====================================================================================                                       |  68%  |                                                                                                                                   |====================================================================================                                       |  69%  |                                                                                                                                   |=====================================================================================                                      |  69%  |                                                                                                                                   |======================================================================================                                     |  70%  |                                                                                                                                   |=======================================================================================                                    |  70%  |                                                                                                                                   |=======================================================================================                                    |  71%  |                                                                                                                                   |========================================================================================                                   |  71%  |                                                                                                                                   |========================================================================================                                   |  72%  |                                                                                                                                   |=========================================================================================                                  |  72%  |                                                                                                                                   |=========================================================================================                                  |  73%  |                                                                                                                                   |==========================================================================================                                 |  73%  |                                                                                                                                   |==========================================================================================                                 |  74%  |                                                                                                                                   |===========================================================================================                                |  74%  |                                                                                                                                   |============================================================================================                               |  74%  |                                                                                                                                   |============================================================================================                               |  75%  |                                                                                                                                   |=============================================================================================                              |  75%  |                                                                                                                                   |=============================================================================================                              |  76%  |                                                                                                                                   |==============================================================================================                             |  76%  |                                                                                                                                   |==============================================================================================                             |  77%  |                                                                                                                                   |===============================================================================================                            |  77%  |                                                                                                                                   |===============================================================================================                            |  78%  |                                                                                                                                   |================================================================================================                           |  78%  |                                                                                                                                   |=================================================================================================                          |  78%  |                                                                                                                                   |=================================================================================================                          |  79%  |                                                                                                                                   |==================================================================================================                         |  79%  |                                                                                                                                   |==================================================================================================                         |  80%  |                                                                                                                                   |===================================================================================================                        |  80%  |                                                                                                                                   |===================================================================================================                        |  81%  |                                                                                                                                   |====================================================================================================                       |  81%  |                                                                                                                                   |====================================================================================================                       |  82%  |                                                                                                                                   |=====================================================================================================                      |  82%  |                                                                                                                                   |======================================================================================================                     |  83%  |                                                                                                                                   |=======================================================================================================                    |  83%  |                                                                                                                                   |=======================================================================================================                    |  84%  |                                                                                                                                   |========================================================================================================                   |  84%  |                                                                                                                                   |========================================================================================================                   |  85%  |                                                                                                                                   |=========================================================================================================                  |  85%  |                                                                                                                                   |=========================================================================================================                  |  86%  |                                                                                                                                   |==========================================================================================================                 |  86%  |                                                                                                                                   |==========================================================================================================                 |  87%  |                                                                                                                                   |===========================================================================================================                |  87%  |                                                                                                                                   |============================================================================================================               |  87%  |                                                                                                                                   |============================================================================================================               |  88%  |                                                                                                                                   |=============================================================================================================              |  88%  |                                                                                                                                   |=============================================================================================================              |  89%  |                                                                                                                                   |==============================================================================================================             |  89%  |                                                                                                                                   |==============================================================================================================             |  90%  |                                                                                                                                   |===============================================================================================================            |  90%  |                                                                                                                                   |===============================================================================================================            |  91%  |                                                                                                                                   |================================================================================================================           |  91%  |                                                                                                                                   |=================================================================================================================          |  92%  |                                                                                                                                   |==================================================================================================================         |  92%  |                                                                                                                                   |==================================================================================================================         |  93%  |                                                                                                                                   |===================================================================================================================        |  93%  |                                                                                                                                   |===================================================================================================================        |  94%  |                                                                                                                                   |====================================================================================================================       |  94%  |                                                                                                                                   |====================================================================================================================       |  95%  |                                                                                                                                   |=====================================================================================================================      |  95%  |                                                                                                                                   |=====================================================================================================================      |  96%  |                                                                                                                                   |======================================================================================================================     |  96%  |                                                                                                                                   |=======================================================================================================================    |  96%  |                                                                                                                                   |=======================================================================================================================    |  97%  |                                                                                                                                   |========================================================================================================================   |  97%  |                                                                                                                                   |========================================================================================================================   |  98%  |                                                                                                                                   |=========================================================================================================================  |  98%  |                                                                                                                                   |=========================================================================================================================  |  99%  |                                                                                                                                   |========================================================================================================================== |  99%  |                                                                                                                                   |========================================================================================================================== | 100%  |                                                                                                                                   |===========================================================================================================================| 100%
```

```
## Removing confirmed isotope entries:
```

```
## Recalculating peaklist for plotting:
```

```
##   |                                                                                                                                   |                                                                                                                           |   0%  |                                                                                                                                   |===============================                                                                                            |  25%  |                                                                                                                                   |==============================================================                                                             |  50%  |                                                                                                                                   |============================================================================================                               |  75%  |                                                                                                                                   |===========================================================================================================================| 100%
```

```
## Calculating adduct similarity:
```

```
##   |                                                                                                                                   |                                                                                                                           |   0%  |                                                                                                                                   |===============================                                                                                            |  25%  |                                                                                                                                   |==============================================================                                                             |  50%  |                                                                                                                                   |============================================================================================                               |  75%  |                                                                                                                                   |===========================================================================================================================| 100%
```

```
## 100 out of 2087
```

```
## 200 out of 2087
```

```
## 300 out of 2087
```

```
## 400 out of 2087
```

```
## 500 out of 2087
```

```
## 600 out of 2087
```

```
## 700 out of 2087
```

```
## 800 out of 2087
```

```
## 900 out of 2087
```

```
## 1000 out of 2087
```

```
## 1100 out of 2087
```

```
## 1200 out of 2087
```

```
## 1300 out of 2087
```

```
## 1400 out of 2087
```

```
## 1500 out of 2087
```

```
## 1600 out of 2087
```

```
## 1700 out of 2087
```

```
## 1800 out of 2087
```

```
## 1900 out of 2087
```

```
## 2000 out of 2087
```

```
##   |                                                                                                                                   |                                                                                                                           |   0%  |                                                                                                                                   |===============================                                                                                            |  25%  |                                                                                                                                   |==============================================================                                                             |  50%  |                                                                                                                                   |============================================================================================                               |  75%  |                                                                                                                                   |===========================================================================================================================| 100%
```

```
## Computing isotopic elution similarity:
```

```
## Calculating isotope similarity:
```

```
##   |                                                                                                                                   |                                                                                                                           |   0%  |                                                                                                                                   |===============================                                                                                            |  25%  |                                                                                                                                   |==============================================================                                                             |  50%  |                                                                                                                                   |============================================================================================                               |  75%  |                                                                                                                                   |===========================================================================================================================| 100%
```

```
## Calculating isotopic fidelity metrics:
```

```
##   |                                                                                                                                   |                                                                                                                           |   0%  |                                                                                                                                   |=                                                                                                                          |   0%  |                                                                                                                                   |=                                                                                                                          |   1%  |                                                                                                                                   |==                                                                                                                         |   1%  |                                                                                                                                   |==                                                                                                                         |   2%  |                                                                                                                                   |===                                                                                                                        |   2%  |                                                                                                                                   |===                                                                                                                        |   3%  |                                                                                                                                   |====                                                                                                                       |   3%  |                                                                                                                                   |====                                                                                                                       |   4%  |                                                                                                                                   |=====                                                                                                                      |   4%  |                                                                                                                                   |======                                                                                                                     |   4%  |                                                                                                                                   |======                                                                                                                     |   5%  |                                                                                                                                   |=======                                                                                                                    |   5%  |                                                                                                                                   |=======                                                                                                                    |   6%  |                                                                                                                                   |========                                                                                                                   |   6%  |                                                                                                                                   |========                                                                                                                   |   7%  |                                                                                                                                   |=========                                                                                                                  |   7%  |                                                                                                                                   |=========                                                                                                                  |   8%  |                                                                                                                                   |==========                                                                                                                 |   8%  |                                                                                                                                   |==========                                                                                                                 |   9%  |                                                                                                                                   |===========                                                                                                                |   9%  |                                                                                                                                   |============                                                                                                               |   9%  |                                                                                                                                   |============                                                                                                               |  10%  |                                                                                                                                   |=============                                                                                                              |  10%  |                                                                                                                                   |=============                                                                                                              |  11%  |                                                                                                                                   |==============                                                                                                             |  11%  |                                                                                                                                   |==============                                                                                                             |  12%  |                                                                                                                                   |===============                                                                                                            |  12%  |                                                                                                                                   |===============                                                                                                            |  13%  |                                                                                                                                   |================                                                                                                           |  13%  |                                                                                                                                   |=================                                                                                                          |  13%  |                                                                                                                                   |=================                                                                                                          |  14%  |                                                                                                                                   |==================                                                                                                         |  14%  |                                                                                                                                   |==================                                                                                                         |  15%  |                                                                                                                                   |===================                                                                                                        |  15%  |                                                                                                                                   |===================                                                                                                        |  16%  |                                                                                                                                   |====================                                                                                                       |  16%  |                                                                                                                                   |====================                                                                                                       |  17%  |                                                                                                                                   |=====================                                                                                                      |  17%  |                                                                                                                                   |======================                                                                                                     |  17%  |                                                                                                                                   |======================                                                                                                     |  18%  |                                                                                                                                   |=======================                                                                                                    |  18%  |                                                                                                                                   |=======================                                                                                                    |  19%  |                                                                                                                                   |========================                                                                                                   |  19%  |                                                                                                                                   |========================                                                                                                   |  20%  |                                                                                                                                   |=========================                                                                                                  |  20%  |                                                                                                                                   |=========================                                                                                                  |  21%  |                                                                                                                                   |==========================                                                                                                 |  21%  |                                                                                                                                   |==========================                                                                                                 |  22%  |                                                                                                                                   |===========================                                                                                                |  22%  |                                                                                                                                   |============================                                                                                               |  22%  |                                                                                                                                   |============================                                                                                               |  23%  |                                                                                                                                   |=============================                                                                                              |  23%  |                                                                                                                                   |=============================                                                                                              |  24%  |                                                                                                                                   |==============================                                                                                             |  24%  |                                                                                                                                   |==============================                                                                                             |  25%  |                                                                                                                                   |===============================                                                                                            |  25%  |                                                                                                                                   |===============================                                                                                            |  26%  |                                                                                                                                   |================================                                                                                           |  26%  |                                                                                                                                   |=================================                                                                                          |  26%  |                                                                                                                                   |=================================                                                                                          |  27%  |                                                                                                                                   |==================================                                                                                         |  27%  |                                                                                                                                   |==================================                                                                                         |  28%  |                                                                                                                                   |===================================                                                                                        |  28%  |                                                                                                                                   |===================================                                                                                        |  29%  |                                                                                                                                   |====================================                                                                                       |  29%  |                                                                                                                                   |====================================                                                                                       |  30%  |                                                                                                                                   |=====================================                                                                                      |  30%  |                                                                                                                                   |======================================                                                                                     |  31%  |                                                                                                                                   |=======================================                                                                                    |  31%  |                                                                                                                                   |=======================================                                                                                    |  32%  |                                                                                                                                   |========================================                                                                                   |  32%  |                                                                                                                                   |========================================                                                                                   |  33%  |                                                                                                                                   |=========================================                                                                                  |  33%  |                                                                                                                                   |=========================================                                                                                  |  34%  |                                                                                                                                   |==========================================                                                                                 |  34%  |                                                                                                                                   |==========================================                                                                                 |  35%  |                                                                                                                                   |===========================================                                                                                |  35%  |                                                                                                                                   |============================================                                                                               |  35%  |                                                                                                                                   |============================================                                                                               |  36%  |                                                                                                                                   |=============================================                                                                              |  36%  |                                                                                                                                   |=============================================                                                                              |  37%  |                                                                                                                                   |==============================================                                                                             |  37%  |                                                                                                                                   |==============================================                                                                             |  38%  |                                                                                                                                   |===============================================                                                                            |  38%  |                                                                                                                                   |===============================================                                                                            |  39%  |                                                                                                                                   |================================================                                                                           |  39%  |                                                                                                                                   |=================================================                                                                          |  39%  |                                                                                                                                   |=================================================                                                                          |  40%  |                                                                                                                                   |==================================================                                                                         |  40%  |                                                                                                                                   |==================================================                                                                         |  41%  |                                                                                                                                   |===================================================                                                                        |  41%  |                                                                                                                                   |===================================================                                                                        |  42%  |                                                                                                                                   |====================================================                                                                       |  42%  |                                                                                                                                   |====================================================                                                                       |  43%  |                                                                                                                                   |=====================================================                                                                      |  43%  |                                                                                                                                   |======================================================                                                                     |  44%  |                                                                                                                                   |=======================================================                                                                    |  44%  |                                                                                                                                   |=======================================================                                                                    |  45%  |                                                                                                                                   |========================================================                                                                   |  45%  |                                                                                                                                   |========================================================                                                                   |  46%  |                                                                                                                                   |=========================================================                                                                  |  46%  |                                                                                                                                   |=========================================================                                                                  |  47%  |                                                                                                                                   |==========================================================                                                                 |  47%  |                                                                                                                                   |==========================================================                                                                 |  48%  |                                                                                                                                   |===========================================================                                                                |  48%  |                                                                                                                                   |============================================================                                                               |  48%  |                                                                                                                                   |============================================================                                                               |  49%  |                                                                                                                                   |=============================================================                                                              |  49%  |                                                                                                                                   |=============================================================                                                              |  50%  |                                                                                                                                   |==============================================================                                                             |  50%  |                                                                                                                                   |==============================================================                                                             |  51%  |                                                                                                                                   |===============================================================                                                            |  51%  |                                                                                                                                   |===============================================================                                                            |  52%  |                                                                                                                                   |================================================================                                                           |  52%  |                                                                                                                                   |=================================================================                                                          |  52%  |                                                                                                                                   |=================================================================                                                          |  53%  |                                                                                                                                   |==================================================================                                                         |  53%  |                                                                                                                                   |==================================================================                                                         |  54%  |                                                                                                                                   |===================================================================                                                        |  54%  |                                                                                                                                   |===================================================================                                                        |  55%  |                                                                                                                                   |====================================================================                                                       |  55%  |                                                                                                                                   |====================================================================                                                       |  56%  |                                                                                                                                   |=====================================================================                                                      |  56%  |                                                                                                                                   |======================================================================                                                     |  57%  |                                                                                                                                   |=======================================================================                                                    |  57%  |                                                                                                                                   |=======================================================================                                                    |  58%  |                                                                                                                                   |========================================================================                                                   |  58%  |                                                                                                                                   |========================================================================                                                   |  59%  |                                                                                                                                   |=========================================================================                                                  |  59%  |                                                                                                                                   |=========================================================================                                                  |  60%  |                                                                                                                                   |==========================================================================                                                 |  60%  |                                                                                                                                   |==========================================================================                                                 |  61%  |                                                                                                                                   |===========================================================================                                                |  61%  |                                                                                                                                   |============================================================================                                               |  61%  |                                                                                                                                   |============================================================================                                               |  62%  |                                                                                                                                   |=============================================================================                                              |  62%  |                                                                                                                                   |=============================================================================                                              |  63%  |                                                                                                                                   |==============================================================================                                             |  63%  |                                                                                                                                   |==============================================================================                                             |  64%  |                                                                                                                                   |===============================================================================                                            |  64%  |                                                                                                                                   |===============================================================================                                            |  65%  |                                                                                                                                   |================================================================================                                           |  65%  |                                                                                                                                   |=================================================================================                                          |  65%  |                                                                                                                                   |=================================================================================                                          |  66%  |                                                                                                                                   |==================================================================================                                         |  66%  |                                                                                                                                   |==================================================================================                                         |  67%  |                                                                                                                                   |===================================================================================                                        |  67%  |                                                                                                                                   |===================================================================================                                        |  68%  |                                                                                                                                   |====================================================================================                                       |  68%  |                                                                                                                                   |====================================================================================                                       |  69%  |                                                                                                                                   |=====================================================================================                                      |  69%  |                                                                                                                                   |======================================================================================                                     |  70%  |                                                                                                                                   |=======================================================================================                                    |  70%  |                                                                                                                                   |=======================================================================================                                    |  71%  |                                                                                                                                   |========================================================================================                                   |  71%  |                                                                                                                                   |========================================================================================                                   |  72%  |                                                                                                                                   |=========================================================================================                                  |  72%  |                                                                                                                                   |=========================================================================================                                  |  73%  |                                                                                                                                   |==========================================================================================                                 |  73%  |                                                                                                                                   |==========================================================================================                                 |  74%  |                                                                                                                                   |===========================================================================================                                |  74%  |                                                                                                                                   |============================================================================================                               |  74%  |                                                                                                                                   |============================================================================================                               |  75%  |                                                                                                                                   |=============================================================================================                              |  75%  |                                                                                                                                   |=============================================================================================                              |  76%  |                                                                                                                                   |==============================================================================================                             |  76%  |                                                                                                                                   |==============================================================================================                             |  77%  |                                                                                                                                   |===============================================================================================                            |  77%  |                                                                                                                                   |===============================================================================================                            |  78%  |                                                                                                                                   |================================================================================================                           |  78%  |                                                                                                                                   |=================================================================================================                          |  78%  |                                                                                                                                   |=================================================================================================                          |  79%  |                                                                                                                                   |==================================================================================================                         |  79%  |                                                                                                                                   |==================================================================================================                         |  80%  |                                                                                                                                   |===================================================================================================                        |  80%  |                                                                                                                                   |===================================================================================================                        |  81%  |                                                                                                                                   |====================================================================================================                       |  81%  |                                                                                                                                   |====================================================================================================                       |  82%  |                                                                                                                                   |=====================================================================================================                      |  82%  |                                                                                                                                   |=====================================================================================================                      |  83%  |                                                                                                                                   |======================================================================================================                     |  83%  |                                                                                                                                   |=======================================================================================================                    |  83%  |                                                                                                                                   |=======================================================================================================                    |  84%  |                                                                                                                                   |========================================================================================================                   |  84%  |                                                                                                                                   |========================================================================================================                   |  85%  |                                                                                                                                   |=========================================================================================================                  |  85%  |                                                                                                                                   |=========================================================================================================                  |  86%  |                                                                                                                                   |==========================================================================================================                 |  86%  |                                                                                                                                   |==========================================================================================================                 |  87%  |                                                                                                                                   |===========================================================================================================                |  87%  |                                                                                                                                   |============================================================================================================               |  87%  |                                                                                                                                   |============================================================================================================               |  88%  |                                                                                                                                   |=============================================================================================================              |  88%  |                                                                                                                                   |=============================================================================================================              |  89%  |                                                                                                                                   |==============================================================================================================             |  89%  |                                                                                                                                   |==============================================================================================================             |  90%  |                                                                                                                                   |===============================================================================================================            |  90%  |                                                                                                                                   |===============================================================================================================            |  91%  |                                                                                                                                   |================================================================================================================           |  91%  |                                                                                                                                   |=================================================================================================================          |  91%  |                                                                                                                                   |=================================================================================================================          |  92%  |                                                                                                                                   |==================================================================================================================         |  92%  |                                                                                                                                   |==================================================================================================================         |  93%  |                                                                                                                                   |===================================================================================================================        |  93%  |                                                                                                                                   |===================================================================================================================        |  94%  |                                                                                                                                   |====================================================================================================================       |  94%  |                                                                                                                                   |====================================================================================================================       |  95%  |                                                                                                                                   |=====================================================================================================================      |  95%  |                                                                                                                                   |=====================================================================================================================      |  96%  |                                                                                                                                   |======================================================================================================================     |  96%  |                                                                                                                                   |=======================================================================================================================    |  96%  |                                                                                                                                   |=======================================================================================================================    |  97%  |                                                                                                                                   |========================================================================================================================   |  97%  |                                                                                                                                   |========================================================================================================================   |  98%  |                                                                                                                                   |=========================================================================================================================  |  98%  |                                                                                                                                   |=========================================================================================================================  |  99%  |                                                                                                                                   |========================================================================================================================== |  99%  |                                                                                                                                   |========================================================================================================================== | 100%  |                                                                                                                                   |===========================================================================================================================| 100%
```

```
## Removing confirmed isotope entries:
```

```
## Recalculating peaklist for plotting:
```

```
##   |                                                                                                                                   |                                                                                                                           |   0%  |                                                                                                                                   |===============================                                                                                            |  25%  |                                                                                                                                   |==============================================================                                                             |  50%  |                                                                                                                                   |============================================================================================                               |  75%  |                                                                                                                                   |===========================================================================================================================| 100%
```

```
## Calculating adduct similarity:
```

```
##   |                                                                                                                                   |                                                                                                                           |   0%  |                                                                                                                                   |===============================                                                                                            |  25%  |                                                                                                                                   |==============================================================                                                             |  50%  |                                                                                                                                   |============================================================================================                               |  75%  |                                                                                                                                   |===========================================================================================================================| 100%
```

```
## 100 out of 3436
```

```
## 200 out of 3436
```

```
## 300 out of 3436
```

```
## 400 out of 3436
```

```
## 500 out of 3436
```

```
## 600 out of 3436
```

```
## 700 out of 3436
```

```
## 800 out of 3436
```

```
## 900 out of 3436
```

```
## 1000 out of 3436
```

```
## 1100 out of 3436
```

```
## 1200 out of 3436
```

```
## 1300 out of 3436
```

```
## 1400 out of 3436
```

```
## 1500 out of 3436
```

```
## 1600 out of 3436
```

```
## 1700 out of 3436
```

```
## 1800 out of 3436
```

```
## 1900 out of 3436
```

```
## 2000 out of 3436
```

```
## 2100 out of 3436
```

```
## 2200 out of 3436
```

```
## 2300 out of 3436
```

```
## 2400 out of 3436
```

```
## 2500 out of 3436
```

```
## 2600 out of 3436
```

```
## 2700 out of 3436
```

```
## 2800 out of 3436
```

```
## 2900 out of 3436
```

```
## 3000 out of 3436
```

```
## 3100 out of 3436
```

```
## 3200 out of 3436
```

```
## 3300 out of 3436
```

```
## 3400 out of 3436
```

```
##   |                                                                                                                                   |                                                                                                                           |   0%  |                                                                                                                                   |===============================                                                                                            |  25%  |                                                                                                                                   |==============================================================                                                             |  50%  |                                                                                                                                   |============================================================================================                               |  75%  |                                                                                                                                   |===========================================================================================================================| 100%
```

```
## Computing isotopic elution similarity:
```

```
## Calculating isotope similarity:
```

```
##   |                                                                                                                                   |                                                                                                                           |   0%  |                                                                                                                                   |===============================                                                                                            |  25%  |                                                                                                                                   |==============================================================                                                             |  50%  |                                                                                                                                   |============================================================================================                               |  75%  |                                                                                                                                   |===========================================================================================================================| 100%
```

```
## Calculating isotopic fidelity metrics:
```

```
##   |                                                                                                                                   |                                                                                                                           |   0%  |                                                                                                                                   |=                                                                                                                          |   0%  |                                                                                                                                   |=                                                                                                                          |   1%  |                                                                                                                                   |==                                                                                                                         |   1%  |                                                                                                                                   |==                                                                                                                         |   2%  |                                                                                                                                   |===                                                                                                                        |   2%  |                                                                                                                                   |===                                                                                                                        |   3%  |                                                                                                                                   |====                                                                                                                       |   3%  |                                                                                                                                   |====                                                                                                                       |   4%  |                                                                                                                                   |=====                                                                                                                      |   4%  |                                                                                                                                   |======                                                                                                                     |   4%  |                                                                                                                                   |======                                                                                                                     |   5%  |                                                                                                                                   |=======                                                                                                                    |   5%  |                                                                                                                                   |=======                                                                                                                    |   6%  |                                                                                                                                   |========                                                                                                                   |   6%  |                                                                                                                                   |========                                                                                                                   |   7%  |                                                                                                                                   |=========                                                                                                                  |   7%  |                                                                                                                                   |=========                                                                                                                  |   8%  |                                                                                                                                   |==========                                                                                                                 |   8%  |                                                                                                                                   |===========                                                                                                                |   9%  |                                                                                                                                   |============                                                                                                               |   9%  |                                                                                                                                   |============                                                                                                               |  10%  |                                                                                                                                   |=============                                                                                                              |  10%  |                                                                                                                                   |=============                                                                                                              |  11%  |                                                                                                                                   |==============                                                                                                             |  11%  |                                                                                                                                   |==============                                                                                                             |  12%  |                                                                                                                                   |===============                                                                                                            |  12%  |                                                                                                                                   |===============                                                                                                            |  13%  |                                                                                                                                   |================                                                                                                           |  13%  |                                                                                                                                   |=================                                                                                                          |  13%  |                                                                                                                                   |=================                                                                                                          |  14%  |                                                                                                                                   |==================                                                                                                         |  14%  |                                                                                                                                   |==================                                                                                                         |  15%  |                                                                                                                                   |===================                                                                                                        |  15%  |                                                                                                                                   |===================                                                                                                        |  16%  |                                                                                                                                   |====================                                                                                                       |  16%  |                                                                                                                                   |====================                                                                                                       |  17%  |                                                                                                                                   |=====================                                                                                                      |  17%  |                                                                                                                                   |======================                                                                                                     |  17%  |                                                                                                                                   |======================                                                                                                     |  18%  |                                                                                                                                   |=======================                                                                                                    |  18%  |                                                                                                                                   |=======================                                                                                                    |  19%  |                                                                                                                                   |========================                                                                                                   |  19%  |                                                                                                                                   |========================                                                                                                   |  20%  |                                                                                                                                   |=========================                                                                                                  |  20%  |                                                                                                                                   |=========================                                                                                                  |  21%  |                                                                                                                                   |==========================                                                                                                 |  21%  |                                                                                                                                   |==========================                                                                                                 |  22%  |                                                                                                                                   |===========================                                                                                                |  22%  |                                                                                                                                   |============================                                                                                               |  22%  |                                                                                                                                   |============================                                                                                               |  23%  |                                                                                                                                   |=============================                                                                                              |  23%  |                                                                                                                                   |=============================                                                                                              |  24%  |                                                                                                                                   |==============================                                                                                             |  24%  |                                                                                                                                   |==============================                                                                                             |  25%  |                                                                                                                                   |===============================                                                                                            |  25%  |                                                                                                                                   |===============================                                                                                            |  26%  |                                                                                                                                   |================================                                                                                           |  26%  |                                                                                                                                   |=================================                                                                                          |  26%  |                                                                                                                                   |=================================                                                                                          |  27%  |                                                                                                                                   |==================================                                                                                         |  27%  |                                                                                                                                   |==================================                                                                                         |  28%  |                                                                                                                                   |===================================                                                                                        |  28%  |                                                                                                                                   |===================================                                                                                        |  29%  |                                                                                                                                   |====================================                                                                                       |  29%  |                                                                                                                                   |====================================                                                                                       |  30%  |                                                                                                                                   |=====================================                                                                                      |  30%  |                                                                                                                                   |======================================                                                                                     |  31%  |                                                                                                                                   |=======================================                                                                                    |  31%  |                                                                                                                                   |=======================================                                                                                    |  32%  |                                                                                                                                   |========================================                                                                                   |  32%  |                                                                                                                                   |========================================                                                                                   |  33%  |                                                                                                                                   |=========================================                                                                                  |  33%  |                                                                                                                                   |=========================================                                                                                  |  34%  |                                                                                                                                   |==========================================                                                                                 |  34%  |                                                                                                                                   |==========================================                                                                                 |  35%  |                                                                                                                                   |===========================================                                                                                |  35%  |                                                                                                                                   |============================================                                                                               |  35%  |                                                                                                                                   |============================================                                                                               |  36%  |                                                                                                                                   |=============================================                                                                              |  36%  |                                                                                                                                   |=============================================                                                                              |  37%  |                                                                                                                                   |==============================================                                                                             |  37%  |                                                                                                                                   |==============================================                                                                             |  38%  |                                                                                                                                   |===============================================                                                                            |  38%  |                                                                                                                                   |===============================================                                                                            |  39%  |                                                                                                                                   |================================================                                                                           |  39%  |                                                                                                                                   |=================================================                                                                          |  39%  |                                                                                                                                   |=================================================                                                                          |  40%  |                                                                                                                                   |==================================================                                                                         |  40%  |                                                                                                                                   |==================================================                                                                         |  41%  |                                                                                                                                   |===================================================                                                                        |  41%  |                                                                                                                                   |===================================================                                                                        |  42%  |                                                                                                                                   |====================================================                                                                       |  42%  |                                                                                                                                   |====================================================                                                                       |  43%  |                                                                                                                                   |=====================================================                                                                      |  43%  |                                                                                                                                   |======================================================                                                                     |  44%  |                                                                                                                                   |=======================================================                                                                    |  44%  |                                                                                                                                   |=======================================================                                                                    |  45%  |                                                                                                                                   |========================================================                                                                   |  45%  |                                                                                                                                   |========================================================                                                                   |  46%  |                                                                                                                                   |=========================================================                                                                  |  46%  |                                                                                                                                   |=========================================================                                                                  |  47%  |                                                                                                                                   |==========================================================                                                                 |  47%  |                                                                                                                                   |==========================================================                                                                 |  48%  |                                                                                                                                   |===========================================================                                                                |  48%  |                                                                                                                                   |============================================================                                                               |  48%  |                                                                                                                                   |============================================================                                                               |  49%  |                                                                                                                                   |=============================================================                                                              |  49%  |                                                                                                                                   |=============================================================                                                              |  50%  |                                                                                                                                   |==============================================================                                                             |  50%  |                                                                                                                                   |==============================================================                                                             |  51%  |                                                                                                                                   |===============================================================                                                            |  51%  |                                                                                                                                   |===============================================================                                                            |  52%  |                                                                                                                                   |================================================================                                                           |  52%  |                                                                                                                                   |=================================================================                                                          |  52%  |                                                                                                                                   |=================================================================                                                          |  53%  |                                                                                                                                   |==================================================================                                                         |  53%  |                                                                                                                                   |==================================================================                                                         |  54%  |                                                                                                                                   |===================================================================                                                        |  54%  |                                                                                                                                   |===================================================================                                                        |  55%  |                                                                                                                                   |====================================================================                                                       |  55%  |                                                                                                                                   |====================================================================                                                       |  56%  |                                                                                                                                   |=====================================================================                                                      |  56%  |                                                                                                                                   |======================================================================                                                     |  57%  |                                                                                                                                   |=======================================================================                                                    |  57%  |                                                                                                                                   |=======================================================================                                                    |  58%  |                                                                                                                                   |========================================================================                                                   |  58%  |                                                                                                                                   |========================================================================                                                   |  59%  |                                                                                                                                   |=========================================================================                                                  |  59%  |                                                                                                                                   |=========================================================================                                                  |  60%  |                                                                                                                                   |==========================================================================                                                 |  60%  |                                                                                                                                   |==========================================================================                                                 |  61%  |                                                                                                                                   |===========================================================================                                                |  61%  |                                                                                                                                   |============================================================================                                               |  61%  |                                                                                                                                   |============================================================================                                               |  62%  |                                                                                                                                   |=============================================================================                                              |  62%  |                                                                                                                                   |=============================================================================                                              |  63%  |                                                                                                                                   |==============================================================================                                             |  63%  |                                                                                                                                   |==============================================================================                                             |  64%  |                                                                                                                                   |===============================================================================                                            |  64%  |                                                                                                                                   |===============================================================================                                            |  65%  |                                                                                                                                   |================================================================================                                           |  65%  |                                                                                                                                   |=================================================================================                                          |  65%  |                                                                                                                                   |=================================================================================                                          |  66%  |                                                                                                                                   |==================================================================================                                         |  66%  |                                                                                                                                   |==================================================================================                                         |  67%  |                                                                                                                                   |===================================================================================                                        |  67%  |                                                                                                                                   |===================================================================================                                        |  68%  |                                                                                                                                   |====================================================================================                                       |  68%  |                                                                                                                                   |====================================================================================                                       |  69%  |                                                                                                                                   |=====================================================================================                                      |  69%  |                                                                                                                                   |======================================================================================                                     |  70%  |                                                                                                                                   |=======================================================================================                                    |  70%  |                                                                                                                                   |=======================================================================================                                    |  71%  |                                                                                                                                   |========================================================================================                                   |  71%  |                                                                                                                                   |========================================================================================                                   |  72%  |                                                                                                                                   |=========================================================================================                                  |  72%  |                                                                                                                                   |=========================================================================================                                  |  73%  |                                                                                                                                   |==========================================================================================                                 |  73%  |                                                                                                                                   |==========================================================================================                                 |  74%  |                                                                                                                                   |===========================================================================================                                |  74%  |                                                                                                                                   |============================================================================================                               |  74%  |                                                                                                                                   |============================================================================================                               |  75%  |                                                                                                                                   |=============================================================================================                              |  75%  |                                                                                                                                   |=============================================================================================                              |  76%  |                                                                                                                                   |==============================================================================================                             |  76%  |                                                                                                                                   |==============================================================================================                             |  77%  |                                                                                                                                   |===============================================================================================                            |  77%  |                                                                                                                                   |===============================================================================================                            |  78%  |                                                                                                                                   |================================================================================================                           |  78%  |                                                                                                                                   |=================================================================================================                          |  78%  |                                                                                                                                   |=================================================================================================                          |  79%  |                                                                                                                                   |==================================================================================================                         |  79%  |                                                                                                                                   |==================================================================================================                         |  80%  |                                                                                                                                   |===================================================================================================                        |  80%  |                                                                                                                                   |===================================================================================================                        |  81%  |                                                                                                                                   |====================================================================================================                       |  81%  |                                                                                                                                   |====================================================================================================                       |  82%  |                                                                                                                                   |=====================================================================================================                      |  82%  |                                                                                                                                   |=====================================================================================================                      |  83%  |                                                                                                                                   |======================================================================================================                     |  83%  |                                                                                                                                   |=======================================================================================================                    |  83%  |                                                                                                                                   |=======================================================================================================                    |  84%  |                                                                                                                                   |========================================================================================================                   |  84%  |                                                                                                                                   |========================================================================================================                   |  85%  |                                                                                                                                   |=========================================================================================================                  |  85%  |                                                                                                                                   |=========================================================================================================                  |  86%  |                                                                                                                                   |==========================================================================================================                 |  86%  |                                                                                                                                   |==========================================================================================================                 |  87%  |                                                                                                                                   |===========================================================================================================                |  87%  |                                                                                                                                   |============================================================================================================               |  87%  |                                                                                                                                   |============================================================================================================               |  88%  |                                                                                                                                   |=============================================================================================================              |  88%  |                                                                                                                                   |=============================================================================================================              |  89%  |                                                                                                                                   |==============================================================================================================             |  89%  |                                                                                                                                   |==============================================================================================================             |  90%  |                                                                                                                                   |===============================================================================================================            |  90%  |                                                                                                                                   |===============================================================================================================            |  91%  |                                                                                                                                   |================================================================================================================           |  91%  |                                                                                                                                   |=================================================================================================================          |  92%  |                                                                                                                                   |==================================================================================================================         |  92%  |                                                                                                                                   |==================================================================================================================         |  93%  |                                                                                                                                   |===================================================================================================================        |  93%  |                                                                                                                                   |===================================================================================================================        |  94%  |                                                                                                                                   |====================================================================================================================       |  94%  |                                                                                                                                   |====================================================================================================================       |  95%  |                                                                                                                                   |=====================================================================================================================      |  95%  |                                                                                                                                   |=====================================================================================================================      |  96%  |                                                                                                                                   |======================================================================================================================     |  96%  |                                                                                                                                   |=======================================================================================================================    |  96%  |                                                                                                                                   |=======================================================================================================================    |  97%  |                                                                                                                                   |========================================================================================================================   |  97%  |                                                                                                                                   |========================================================================================================================   |  98%  |                                                                                                                                   |=========================================================================================================================  |  98%  |                                                                                                                                   |=========================================================================================================================  |  99%  |                                                                                                                                   |========================================================================================================================== |  99%  |                                                                                                                                   |========================================================================================================================== | 100%  |                                                                                                                                   |===========================================================================================================================| 100%
```

```
## Removing confirmed isotope entries:
```

```
## Recalculating peaklist for plotting:
```

```
##   |                                                                                                                                   |                                                                                                                           |   0%  |                                                                                                                                   |===============================                                                                                            |  25%  |                                                                                                                                   |==============================================================                                                             |  50%  |                                                                                                                                   |============================================================================================                               |  75%  |                                                                                                                                   |===========================================================================================================================| 100%
```

```
## Calculating adduct similarity:
```

```
##   |                                                                                                                                   |                                                                                                                           |   0%  |                                                                                                                                   |===============================                                                                                            |  25%  |                                                                                                                                   |==============================================================                                                             |  50%  |                                                                                                                                   |============================================================================================                               |  75%  |                                                                                                                                   |===========================================================================================================================| 100%
```

```
## 100 out of 1851
```

```
## 200 out of 1851
```

```
## 300 out of 1851
```

```
## 400 out of 1851
```

```
## 500 out of 1851
```

```
## 600 out of 1851
```

```
## 700 out of 1851
```

```
## 800 out of 1851
```

```
## 900 out of 1851
```

```
## 1000 out of 1851
```

```
## 1100 out of 1851
```

```
## 1200 out of 1851
```

```
## 1300 out of 1851
```

```
## 1400 out of 1851
```

```
## 1500 out of 1851
```

```
## 1600 out of 1851
```

```
## 1700 out of 1851
```

```
## 1800 out of 1851
```

```
##   |                                                                                                                                   |                                                                                                                           |   0%  |                                                                                                                                   |===============================                                                                                            |  25%  |                                                                                                                                   |==============================================================                                                             |  50%  |                                                                                                                                   |============================================================================================                               |  75%  |                                                                                                                                   |===========================================================================================================================| 100%
```

```
## Computing isotopic elution similarity:
```

```
## Calculating isotope similarity:
```

```
##   |                                                                                                                                   |                                                                                                                           |   0%  |                                                                                                                                   |===============================                                                                                            |  25%  |                                                                                                                                   |==============================================================                                                             |  50%  |                                                                                                                                   |============================================================================================                               |  75%  |                                                                                                                                   |===========================================================================================================================| 100%
```

```
## Calculating isotopic fidelity metrics:
```

```
##   |                                                                                                                                   |                                                                                                                           |   0%  |                                                                                                                                   |=                                                                                                                          |   0%  |                                                                                                                                   |=                                                                                                                          |   1%  |                                                                                                                                   |==                                                                                                                         |   1%  |                                                                                                                                   |==                                                                                                                         |   2%  |                                                                                                                                   |===                                                                                                                        |   2%  |                                                                                                                                   |===                                                                                                                        |   3%  |                                                                                                                                   |====                                                                                                                       |   3%  |                                                                                                                                   |====                                                                                                                       |   4%  |                                                                                                                                   |=====                                                                                                                      |   4%  |                                                                                                                                   |======                                                                                                                     |   4%  |                                                                                                                                   |======                                                                                                                     |   5%  |                                                                                                                                   |=======                                                                                                                    |   5%  |                                                                                                                                   |=======                                                                                                                    |   6%  |                                                                                                                                   |========                                                                                                                   |   6%  |                                                                                                                                   |========                                                                                                                   |   7%  |                                                                                                                                   |=========                                                                                                                  |   7%  |                                                                                                                                   |=========                                                                                                                  |   8%  |                                                                                                                                   |==========                                                                                                                 |   8%  |                                                                                                                                   |==========                                                                                                                 |   9%  |                                                                                                                                   |===========                                                                                                                |   9%  |                                                                                                                                   |============                                                                                                               |   9%  |                                                                                                                                   |============                                                                                                               |  10%  |                                                                                                                                   |=============                                                                                                              |  10%  |                                                                                                                                   |=============                                                                                                              |  11%  |                                                                                                                                   |==============                                                                                                             |  11%  |                                                                                                                                   |==============                                                                                                             |  12%  |                                                                                                                                   |===============                                                                                                            |  12%  |                                                                                                                                   |===============                                                                                                            |  13%  |                                                                                                                                   |================                                                                                                           |  13%  |                                                                                                                                   |=================                                                                                                          |  13%  |                                                                                                                                   |=================                                                                                                          |  14%  |                                                                                                                                   |==================                                                                                                         |  14%  |                                                                                                                                   |==================                                                                                                         |  15%  |                                                                                                                                   |===================                                                                                                        |  15%  |                                                                                                                                   |===================                                                                                                        |  16%  |                                                                                                                                   |====================                                                                                                       |  16%  |                                                                                                                                   |====================                                                                                                       |  17%  |                                                                                                                                   |=====================                                                                                                      |  17%  |                                                                                                                                   |======================                                                                                                     |  17%  |                                                                                                                                   |======================                                                                                                     |  18%  |                                                                                                                                   |=======================                                                                                                    |  18%  |                                                                                                                                   |=======================                                                                                                    |  19%  |                                                                                                                                   |========================                                                                                                   |  19%  |                                                                                                                                   |========================                                                                                                   |  20%  |                                                                                                                                   |=========================                                                                                                  |  20%  |                                                                                                                                   |=========================                                                                                                  |  21%  |                                                                                                                                   |==========================                                                                                                 |  21%  |                                                                                                                                   |==========================                                                                                                 |  22%  |                                                                                                                                   |===========================                                                                                                |  22%  |                                                                                                                                   |============================                                                                                               |  22%  |                                                                                                                                   |============================                                                                                               |  23%  |                                                                                                                                   |=============================                                                                                              |  23%  |                                                                                                                                   |=============================                                                                                              |  24%  |                                                                                                                                   |==============================                                                                                             |  24%  |                                                                                                                                   |==============================                                                                                             |  25%  |                                                                                                                                   |===============================                                                                                            |  25%  |                                                                                                                                   |===============================                                                                                            |  26%  |                                                                                                                                   |================================                                                                                           |  26%  |                                                                                                                                   |=================================                                                                                          |  26%  |                                                                                                                                   |=================================                                                                                          |  27%  |                                                                                                                                   |==================================                                                                                         |  27%  |                                                                                                                                   |==================================                                                                                         |  28%  |                                                                                                                                   |===================================                                                                                        |  28%  |                                                                                                                                   |===================================                                                                                        |  29%  |                                                                                                                                   |====================================                                                                                       |  29%  |                                                                                                                                   |====================================                                                                                       |  30%  |                                                                                                                                   |=====================================                                                                                      |  30%  |                                                                                                                                   |======================================                                                                                     |  31%  |                                                                                                                                   |=======================================                                                                                    |  31%  |                                                                                                                                   |=======================================                                                                                    |  32%  |                                                                                                                                   |========================================                                                                                   |  32%  |                                                                                                                                   |========================================                                                                                   |  33%  |                                                                                                                                   |=========================================                                                                                  |  33%  |                                                                                                                                   |=========================================                                                                                  |  34%  |                                                                                                                                   |==========================================                                                                                 |  34%  |                                                                                                                                   |==========================================                                                                                 |  35%  |                                                                                                                                   |===========================================                                                                                |  35%  |                                                                                                                                   |============================================                                                                               |  35%  |                                                                                                                                   |============================================                                                                               |  36%  |                                                                                                                                   |=============================================                                                                              |  36%  |                                                                                                                                   |=============================================                                                                              |  37%  |                                                                                                                                   |==============================================                                                                             |  37%  |                                                                                                                                   |==============================================                                                                             |  38%  |                                                                                                                                   |===============================================                                                                            |  38%  |                                                                                                                                   |===============================================                                                                            |  39%  |                                                                                                                                   |================================================                                                                           |  39%  |                                                                                                                                   |=================================================                                                                          |  39%  |                                                                                                                                   |=================================================                                                                          |  40%  |                                                                                                                                   |==================================================                                                                         |  40%  |                                                                                                                                   |==================================================                                                                         |  41%  |                                                                                                                                   |===================================================                                                                        |  41%  |                                                                                                                                   |===================================================                                                                        |  42%  |                                                                                                                                   |====================================================                                                                       |  42%  |                                                                                                                                   |====================================================                                                                       |  43%  |                                                                                                                                   |=====================================================                                                                      |  43%  |                                                                                                                                   |======================================================                                                                     |  44%  |                                                                                                                                   |=======================================================                                                                    |  44%  |                                                                                                                                   |=======================================================                                                                    |  45%  |                                                                                                                                   |========================================================                                                                   |  45%  |                                                                                                                                   |========================================================                                                                   |  46%  |                                                                                                                                   |=========================================================                                                                  |  46%  |                                                                                                                                   |=========================================================                                                                  |  47%  |                                                                                                                                   |==========================================================                                                                 |  47%  |                                                                                                                                   |==========================================================                                                                 |  48%  |                                                                                                                                   |===========================================================                                                                |  48%  |                                                                                                                                   |============================================================                                                               |  48%  |                                                                                                                                   |============================================================                                                               |  49%  |                                                                                                                                   |=============================================================                                                              |  49%  |                                                                                                                                   |=============================================================                                                              |  50%  |                                                                                                                                   |==============================================================                                                             |  50%  |                                                                                                                                   |==============================================================                                                             |  51%  |                                                                                                                                   |===============================================================                                                            |  51%  |                                                                                                                                   |===============================================================                                                            |  52%  |                                                                                                                                   |================================================================                                                           |  52%  |                                                                                                                                   |=================================================================                                                          |  52%  |                                                                                                                                   |=================================================================                                                          |  53%  |                                                                                                                                   |==================================================================                                                         |  53%  |                                                                                                                                   |==================================================================                                                         |  54%  |                                                                                                                                   |===================================================================                                                        |  54%  |                                                                                                                                   |===================================================================                                                        |  55%  |                                                                                                                                   |====================================================================                                                       |  55%  |                                                                                                                                   |====================================================================                                                       |  56%  |                                                                                                                                   |=====================================================================                                                      |  56%  |                                                                                                                                   |======================================================================                                                     |  57%  |                                                                                                                                   |=======================================================================                                                    |  57%  |                                                                                                                                   |=======================================================================                                                    |  58%  |                                                                                                                                   |========================================================================                                                   |  58%  |                                                                                                                                   |========================================================================                                                   |  59%  |                                                                                                                                   |=========================================================================                                                  |  59%  |                                                                                                                                   |=========================================================================                                                  |  60%  |                                                                                                                                   |==========================================================================                                                 |  60%  |                                                                                                                                   |==========================================================================                                                 |  61%  |                                                                                                                                   |===========================================================================                                                |  61%  |                                                                                                                                   |============================================================================                                               |  61%  |                                                                                                                                   |============================================================================                                               |  62%  |                                                                                                                                   |=============================================================================                                              |  62%  |                                                                                                                                   |=============================================================================                                              |  63%  |                                                                                                                                   |==============================================================================                                             |  63%  |                                                                                                                                   |==============================================================================                                             |  64%  |                                                                                                                                   |===============================================================================                                            |  64%  |                                                                                                                                   |===============================================================================                                            |  65%  |                                                                                                                                   |================================================================================                                           |  65%  |                                                                                                                                   |=================================================================================                                          |  65%  |                                                                                                                                   |=================================================================================                                          |  66%  |                                                                                                                                   |==================================================================================                                         |  66%  |                                                                                                                                   |==================================================================================                                         |  67%  |                                                                                                                                   |===================================================================================                                        |  67%  |                                                                                                                                   |===================================================================================                                        |  68%  |                                                                                                                                   |====================================================================================                                       |  68%  |                                                                                                                                   |====================================================================================                                       |  69%  |                                                                                                                                   |=====================================================================================                                      |  69%  |                                                                                                                                   |======================================================================================                                     |  70%  |                                                                                                                                   |=======================================================================================                                    |  70%  |                                                                                                                                   |=======================================================================================                                    |  71%  |                                                                                                                                   |========================================================================================                                   |  71%  |                                                                                                                                   |========================================================================================                                   |  72%  |                                                                                                                                   |=========================================================================================                                  |  72%  |                                                                                                                                   |=========================================================================================                                  |  73%  |                                                                                                                                   |==========================================================================================                                 |  73%  |                                                                                                                                   |==========================================================================================                                 |  74%  |                                                                                                                                   |===========================================================================================                                |  74%  |                                                                                                                                   |============================================================================================                               |  74%  |                                                                                                                                   |============================================================================================                               |  75%  |                                                                                                                                   |=============================================================================================                              |  75%  |                                                                                                                                   |=============================================================================================                              |  76%  |                                                                                                                                   |==============================================================================================                             |  76%  |                                                                                                                                   |==============================================================================================                             |  77%  |                                                                                                                                   |===============================================================================================                            |  77%  |                                                                                                                                   |===============================================================================================                            |  78%  |                                                                                                                                   |================================================================================================                           |  78%  |                                                                                                                                   |=================================================================================================                          |  78%  |                                                                                                                                   |=================================================================================================                          |  79%  |                                                                                                                                   |==================================================================================================                         |  79%  |                                                                                                                                   |==================================================================================================                         |  80%  |                                                                                                                                   |===================================================================================================                        |  80%  |                                                                                                                                   |===================================================================================================                        |  81%  |                                                                                                                                   |====================================================================================================                       |  81%  |                                                                                                                                   |====================================================================================================                       |  82%  |                                                                                                                                   |=====================================================================================================                      |  82%  |                                                                                                                                   |=====================================================================================================                      |  83%  |                                                                                                                                   |======================================================================================================                     |  83%  |                                                                                                                                   |=======================================================================================================                    |  83%  |                                                                                                                                   |=======================================================================================================                    |  84%  |                                                                                                                                   |========================================================================================================                   |  84%  |                                                                                                                                   |========================================================================================================                   |  85%  |                                                                                                                                   |=========================================================================================================                  |  85%  |                                                                                                                                   |=========================================================================================================                  |  86%  |                                                                                                                                   |==========================================================================================================                 |  86%  |                                                                                                                                   |==========================================================================================================                 |  87%  |                                                                                                                                   |===========================================================================================================                |  87%  |                                                                                                                                   |============================================================================================================               |  87%  |                                                                                                                                   |============================================================================================================               |  88%  |                                                                                                                                   |=============================================================================================================              |  88%  |                                                                                                                                   |=============================================================================================================              |  89%  |                                                                                                                                   |==============================================================================================================             |  89%  |                                                                                                                                   |==============================================================================================================             |  90%  |                                                                                                                                   |===============================================================================================================            |  90%  |                                                                                                                                   |===============================================================================================================            |  91%  |                                                                                                                                   |================================================================================================================           |  91%  |                                                                                                                                   |=================================================================================================================          |  91%  |                                                                                                                                   |=================================================================================================================          |  92%  |                                                                                                                                   |==================================================================================================================         |  92%  |                                                                                                                                   |==================================================================================================================         |  93%  |                                                                                                                                   |===================================================================================================================        |  93%  |                                                                                                                                   |===================================================================================================================        |  94%  |                                                                                                                                   |====================================================================================================================       |  94%  |                                                                                                                                   |====================================================================================================================       |  95%  |                                                                                                                                   |=====================================================================================================================      |  95%  |                                                                                                                                   |=====================================================================================================================      |  96%  |                                                                                                                                   |======================================================================================================================     |  96%  |                                                                                                                                   |=======================================================================================================================    |  96%  |                                                                                                                                   |=======================================================================================================================    |  97%  |                                                                                                                                   |========================================================================================================================   |  97%  |                                                                                                                                   |========================================================================================================================   |  98%  |                                                                                                                                   |=========================================================================================================================  |  98%  |                                                                                                                                   |=========================================================================================================================  |  99%  |                                                                                                                                   |========================================================================================================================== |  99%  |                                                                                                                                   |========================================================================================================================== | 100%  |                                                                                                                                   |===========================================================================================================================| 100%
```

```
## Removing confirmed isotope entries:
```

```
## Recalculating peaklist for plotting:
```

```
##   |                                                                                                                                   |                                                                                                                           |   0%  |                                                                                                                                   |===============================                                                                                            |  25%  |                                                                                                                                   |==============================================================                                                             |  50%  |                                                                                                                                   |============================================================================================                               |  75%  |                                                                                                                                   |===========================================================================================================================| 100%
```

```
## Calculating adduct similarity:
```

```
##   |                                                                                                                                   |                                                                                                                           |   0%  |                                                                                                                                   |===============================                                                                                            |  25%  |                                                                                                                                   |==============================================================                                                             |  50%  |                                                                                                                                   |============================================================================================                               |  75%  |                                                                                                                                   |===========================================================================================================================| 100%
```

```
## 100 out of 3290
```

```
## 200 out of 3290
```

```
## 300 out of 3290
```

```
## 400 out of 3290
```

```
## 500 out of 3290
```

```
## 600 out of 3290
```

```
## 700 out of 3290
```

```
## 800 out of 3290
```

```
## 900 out of 3290
```

```
## 1000 out of 3290
```

```
## 1100 out of 3290
```

```
## 1200 out of 3290
```

```
## 1300 out of 3290
```

```
## 1400 out of 3290
```

```
## 1500 out of 3290
```

```
## 1600 out of 3290
```

```
## 1700 out of 3290
```

```
## 1800 out of 3290
```

```
## 1900 out of 3290
```

```
## 2000 out of 3290
```

```
## 2100 out of 3290
```

```
## 2200 out of 3290
```

```
## 2300 out of 3290
```

```
## 2400 out of 3290
```

```
## 2500 out of 3290
```

```
## 2600 out of 3290
```

```
## 2700 out of 3290
```

```
## 2800 out of 3290
```

```
## 2900 out of 3290
```

```
## 3000 out of 3290
```

```
## 3100 out of 3290
```

```
## 3200 out of 3290
```

```
##   |                                                                                                                                   |                                                                                                                           |   0%  |                                                                                                                                   |===============================                                                                                            |  25%  |                                                                                                                                   |==============================================================                                                             |  50%  |                                                                                                                                   |============================================================================================                               |  75%  |                                                                                                                                   |===========================================================================================================================| 100%
```

```
## Computing isotopic elution similarity:
```

```
## Calculating isotope similarity:
```

```
##   |                                                                                                                                   |                                                                                                                           |   0%  |                                                                                                                                   |===============================                                                                                            |  25%  |                                                                                                                                   |==============================================================                                                             |  50%  |                                                                                                                                   |============================================================================================                               |  75%  |                                                                                                                                   |===========================================================================================================================| 100%
```

```
## Calculating isotopic fidelity metrics:
```

```
##   |                                                                                                                                   |                                                                                                                           |   0%  |                                                                                                                                   |=                                                                                                                          |   0%  |                                                                                                                                   |=                                                                                                                          |   1%  |                                                                                                                                   |==                                                                                                                         |   1%  |                                                                                                                                   |==                                                                                                                         |   2%  |                                                                                                                                   |===                                                                                                                        |   2%  |                                                                                                                                   |===                                                                                                                        |   3%  |                                                                                                                                   |====                                                                                                                       |   3%  |                                                                                                                                   |====                                                                                                                       |   4%  |                                                                                                                                   |=====                                                                                                                      |   4%  |                                                                                                                                   |======                                                                                                                     |   4%  |                                                                                                                                   |======                                                                                                                     |   5%  |                                                                                                                                   |=======                                                                                                                    |   5%  |                                                                                                                                   |=======                                                                                                                    |   6%  |                                                                                                                                   |========                                                                                                                   |   6%  |                                                                                                                                   |========                                                                                                                   |   7%  |                                                                                                                                   |=========                                                                                                                  |   7%  |                                                                                                                                   |=========                                                                                                                  |   8%  |                                                                                                                                   |==========                                                                                                                 |   8%  |                                                                                                                                   |==========                                                                                                                 |   9%  |                                                                                                                                   |===========                                                                                                                |   9%  |                                                                                                                                   |============                                                                                                               |   9%  |                                                                                                                                   |============                                                                                                               |  10%  |                                                                                                                                   |=============                                                                                                              |  10%  |                                                                                                                                   |=============                                                                                                              |  11%  |                                                                                                                                   |==============                                                                                                             |  11%  |                                                                                                                                   |==============                                                                                                             |  12%  |                                                                                                                                   |===============                                                                                                            |  12%  |                                                                                                                                   |===============                                                                                                            |  13%  |                                                                                                                                   |================                                                                                                           |  13%  |                                                                                                                                   |=================                                                                                                          |  13%  |                                                                                                                                   |=================                                                                                                          |  14%  |                                                                                                                                   |==================                                                                                                         |  14%  |                                                                                                                                   |==================                                                                                                         |  15%  |                                                                                                                                   |===================                                                                                                        |  15%  |                                                                                                                                   |===================                                                                                                        |  16%  |                                                                                                                                   |====================                                                                                                       |  16%  |                                                                                                                                   |====================                                                                                                       |  17%  |                                                                                                                                   |=====================                                                                                                      |  17%  |                                                                                                                                   |======================                                                                                                     |  17%  |                                                                                                                                   |======================                                                                                                     |  18%  |                                                                                                                                   |=======================                                                                                                    |  18%  |                                                                                                                                   |=======================                                                                                                    |  19%  |                                                                                                                                   |========================                                                                                                   |  19%  |                                                                                                                                   |========================                                                                                                   |  20%  |                                                                                                                                   |=========================                                                                                                  |  20%  |                                                                                                                                   |=========================                                                                                                  |  21%  |                                                                                                                                   |==========================                                                                                                 |  21%  |                                                                                                                                   |==========================                                                                                                 |  22%  |                                                                                                                                   |===========================                                                                                                |  22%  |                                                                                                                                   |============================                                                                                               |  22%  |                                                                                                                                   |============================                                                                                               |  23%  |                                                                                                                                   |=============================                                                                                              |  23%  |                                                                                                                                   |=============================                                                                                              |  24%  |                                                                                                                                   |==============================                                                                                             |  24%  |                                                                                                                                   |==============================                                                                                             |  25%  |                                                                                                                                   |===============================                                                                                            |  25%  |                                                                                                                                   |===============================                                                                                            |  26%  |                                                                                                                                   |================================                                                                                           |  26%  |                                                                                                                                   |=================================                                                                                          |  26%  |                                                                                                                                   |=================================                                                                                          |  27%  |                                                                                                                                   |==================================                                                                                         |  27%  |                                                                                                                                   |==================================                                                                                         |  28%  |                                                                                                                                   |===================================                                                                                        |  28%  |                                                                                                                                   |===================================                                                                                        |  29%  |                                                                                                                                   |====================================                                                                                       |  29%  |                                                                                                                                   |====================================                                                                                       |  30%  |                                                                                                                                   |=====================================                                                                                      |  30%  |                                                                                                                                   |======================================                                                                                     |  31%  |                                                                                                                                   |=======================================                                                                                    |  31%  |                                                                                                                                   |=======================================                                                                                    |  32%  |                                                                                                                                   |========================================                                                                                   |  32%  |                                                                                                                                   |========================================                                                                                   |  33%  |                                                                                                                                   |=========================================                                                                                  |  33%  |                                                                                                                                   |=========================================                                                                                  |  34%  |                                                                                                                                   |==========================================                                                                                 |  34%  |                                                                                                                                   |==========================================                                                                                 |  35%  |                                                                                                                                   |===========================================                                                                                |  35%  |                                                                                                                                   |============================================                                                                               |  35%  |                                                                                                                                   |============================================                                                                               |  36%  |                                                                                                                                   |=============================================                                                                              |  36%  |                                                                                                                                   |=============================================                                                                              |  37%  |                                                                                                                                   |==============================================                                                                             |  37%  |                                                                                                                                   |==============================================                                                                             |  38%  |                                                                                                                                   |===============================================                                                                            |  38%  |                                                                                                                                   |===============================================                                                                            |  39%  |                                                                                                                                   |================================================                                                                           |  39%  |                                                                                                                                   |=================================================                                                                          |  39%  |                                                                                                                                   |=================================================                                                                          |  40%  |                                                                                                                                   |==================================================                                                                         |  40%  |                                                                                                                                   |==================================================                                                                         |  41%  |                                                                                                                                   |===================================================                                                                        |  41%  |                                                                                                                                   |===================================================                                                                        |  42%  |                                                                                                                                   |====================================================                                                                       |  42%  |                                                                                                                                   |====================================================                                                                       |  43%  |                                                                                                                                   |=====================================================                                                                      |  43%  |                                                                                                                                   |======================================================                                                                     |  44%  |                                                                                                                                   |=======================================================                                                                    |  44%  |                                                                                                                                   |=======================================================                                                                    |  45%  |                                                                                                                                   |========================================================                                                                   |  45%  |                                                                                                                                   |========================================================                                                                   |  46%  |                                                                                                                                   |=========================================================                                                                  |  46%  |                                                                                                                                   |=========================================================                                                                  |  47%  |                                                                                                                                   |==========================================================                                                                 |  47%  |                                                                                                                                   |==========================================================                                                                 |  48%  |                                                                                                                                   |===========================================================                                                                |  48%  |                                                                                                                                   |============================================================                                                               |  48%  |                                                                                                                                   |============================================================                                                               |  49%  |                                                                                                                                   |=============================================================                                                              |  49%  |                                                                                                                                   |=============================================================                                                              |  50%  |                                                                                                                                   |==============================================================                                                             |  50%  |                                                                                                                                   |==============================================================                                                             |  51%  |                                                                                                                                   |===============================================================                                                            |  51%  |                                                                                                                                   |===============================================================                                                            |  52%  |                                                                                                                                   |================================================================                                                           |  52%  |                                                                                                                                   |=================================================================                                                          |  52%  |                                                                                                                                   |=================================================================                                                          |  53%  |                                                                                                                                   |==================================================================                                                         |  53%  |                                                                                                                                   |==================================================================                                                         |  54%  |                                                                                                                                   |===================================================================                                                        |  54%  |                                                                                                                                   |===================================================================                                                        |  55%  |                                                                                                                                   |====================================================================                                                       |  55%  |                                                                                                                                   |====================================================================                                                       |  56%  |                                                                                                                                   |=====================================================================                                                      |  56%  |                                                                                                                                   |======================================================================                                                     |  57%  |                                                                                                                                   |=======================================================================                                                    |  57%  |                                                                                                                                   |=======================================================================                                                    |  58%  |                                                                                                                                   |========================================================================                                                   |  58%  |                                                                                                                                   |========================================================================                                                   |  59%  |                                                                                                                                   |=========================================================================                                                  |  59%  |                                                                                                                                   |=========================================================================                                                  |  60%  |                                                                                                                                   |==========================================================================                                                 |  60%  |                                                                                                                                   |==========================================================================                                                 |  61%  |                                                                                                                                   |===========================================================================                                                |  61%  |                                                                                                                                   |============================================================================                                               |  61%  |                                                                                                                                   |============================================================================                                               |  62%  |                                                                                                                                   |=============================================================================                                              |  62%  |                                                                                                                                   |=============================================================================                                              |  63%  |                                                                                                                                   |==============================================================================                                             |  63%  |                                                                                                                                   |==============================================================================                                             |  64%  |                                                                                                                                   |===============================================================================                                            |  64%  |                                                                                                                                   |===============================================================================                                            |  65%  |                                                                                                                                   |================================================================================                                           |  65%  |                                                                                                                                   |=================================================================================                                          |  65%  |                                                                                                                                   |=================================================================================                                          |  66%  |                                                                                                                                   |==================================================================================                                         |  66%  |                                                                                                                                   |==================================================================================                                         |  67%  |                                                                                                                                   |===================================================================================                                        |  67%  |                                                                                                                                   |===================================================================================                                        |  68%  |                                                                                                                                   |====================================================================================                                       |  68%  |                                                                                                                                   |====================================================================================                                       |  69%  |                                                                                                                                   |=====================================================================================                                      |  69%  |                                                                                                                                   |======================================================================================                                     |  70%  |                                                                                                                                   |=======================================================================================                                    |  70%  |                                                                                                                                   |=======================================================================================                                    |  71%  |                                                                                                                                   |========================================================================================                                   |  71%  |                                                                                                                                   |========================================================================================                                   |  72%  |                                                                                                                                   |=========================================================================================                                  |  72%  |                                                                                                                                   |=========================================================================================                                  |  73%  |                                                                                                                                   |==========================================================================================                                 |  73%  |                                                                                                                                   |==========================================================================================                                 |  74%  |                                                                                                                                   |===========================================================================================                                |  74%  |                                                                                                                                   |============================================================================================                               |  74%  |                                                                                                                                   |============================================================================================                               |  75%  |                                                                                                                                   |=============================================================================================                              |  75%  |                                                                                                                                   |=============================================================================================                              |  76%  |                                                                                                                                   |==============================================================================================                             |  76%  |                                                                                                                                   |==============================================================================================                             |  77%  |                                                                                                                                   |===============================================================================================                            |  77%  |                                                                                                                                   |===============================================================================================                            |  78%  |                                                                                                                                   |================================================================================================                           |  78%  |                                                                                                                                   |=================================================================================================                          |  78%  |                                                                                                                                   |=================================================================================================                          |  79%  |                                                                                                                                   |==================================================================================================                         |  79%  |                                                                                                                                   |==================================================================================================                         |  80%  |                                                                                                                                   |===================================================================================================                        |  80%  |                                                                                                                                   |===================================================================================================                        |  81%  |                                                                                                                                   |====================================================================================================                       |  81%  |                                                                                                                                   |====================================================================================================                       |  82%  |                                                                                                                                   |=====================================================================================================                      |  82%  |                                                                                                                                   |=====================================================================================================                      |  83%  |                                                                                                                                   |======================================================================================================                     |  83%  |                                                                                                                                   |=======================================================================================================                    |  83%  |                                                                                                                                   |=======================================================================================================                    |  84%  |                                                                                                                                   |========================================================================================================                   |  84%  |                                                                                                                                   |========================================================================================================                   |  85%  |                                                                                                                                   |=========================================================================================================                  |  85%  |                                                                                                                                   |=========================================================================================================                  |  86%  |                                                                                                                                   |==========================================================================================================                 |  86%  |                                                                                                                                   |==========================================================================================================                 |  87%  |                                                                                                                                   |===========================================================================================================                |  87%  |                                                                                                                                   |============================================================================================================               |  87%  |                                                                                                                                   |============================================================================================================               |  88%  |                                                                                                                                   |=============================================================================================================              |  88%  |                                                                                                                                   |=============================================================================================================              |  89%  |                                                                                                                                   |==============================================================================================================             |  89%  |                                                                                                                                   |==============================================================================================================             |  90%  |                                                                                                                                   |===============================================================================================================            |  90%  |                                                                                                                                   |===============================================================================================================            |  91%  |                                                                                                                                   |================================================================================================================           |  91%  |                                                                                                                                   |=================================================================================================================          |  91%  |                                                                                                                                   |=================================================================================================================          |  92%  |                                                                                                                                   |==================================================================================================================         |  92%  |                                                                                                                                   |==================================================================================================================         |  93%  |                                                                                                                                   |===================================================================================================================        |  93%  |                                                                                                                                   |===================================================================================================================        |  94%  |                                                                                                                                   |====================================================================================================================       |  94%  |                                                                                                                                   |====================================================================================================================       |  95%  |                                                                                                                                   |=====================================================================================================================      |  95%  |                                                                                                                                   |=====================================================================================================================      |  96%  |                                                                                                                                   |======================================================================================================================     |  96%  |                                                                                                                                   |=======================================================================================================================    |  96%  |                                                                                                                                   |=======================================================================================================================    |  97%  |                                                                                                                                   |========================================================================================================================   |  97%  |                                                                                                                                   |========================================================================================================================   |  98%  |                                                                                                                                   |=========================================================================================================================  |  98%  |                                                                                                                                   |=========================================================================================================================  |  99%  |                                                                                                                                   |========================================================================================================================== |  99%  |                                                                                                                                   |========================================================================================================================== | 100%  |                                                                                                                                   |===========================================================================================================================| 100%
```

```
## Removing confirmed isotope entries:
```

```
## Recalculating peaklist for plotting:
```

```
##   |                                                                                                                                   |                                                                                                                           |   0%  |                                                                                                                                   |===============================                                                                                            |  25%  |                                                                                                                                   |==============================================================                                                             |  50%  |                                                                                                                                   |============================================================================================                               |  75%  |                                                                                                                                   |===========================================================================================================================| 100%
```

```
## Calculating adduct similarity:
```

```
##   |                                                                                                                                   |                                                                                                                           |   0%  |                                                                                                                                   |===============================                                                                                            |  25%  |                                                                                                                                   |==============================================================                                                             |  50%  |                                                                                                                                   |============================================================================================                               |  75%  |                                                                                                                                   |===========================================================================================================================| 100%
```

```
## 100 out of 1824
```

```
## 200 out of 1824
```

```
## 300 out of 1824
```

```
## 400 out of 1824
```

```
## 500 out of 1824
```

```
## 600 out of 1824
```

```
## 700 out of 1824
```

```
## 800 out of 1824
```

```
## 900 out of 1824
```

```
## 1000 out of 1824
```

```
## 1100 out of 1824
```

```
## 1200 out of 1824
```

```
## 1300 out of 1824
```

```
## 1400 out of 1824
```

```
## 1500 out of 1824
```

```
## 1600 out of 1824
```

```
## 1700 out of 1824
```

```
## 1800 out of 1824
```

```
##   |                                                                                                                                   |                                                                                                                           |   0%  |                                                                                                                                   |===============================                                                                                            |  25%  |                                                                                                                                   |==============================================================                                                             |  50%  |                                                                                                                                   |============================================================================================                               |  75%  |                                                                                                                                   |===========================================================================================================================| 100%
```

```
## Computing isotopic elution similarity:
```

```
## Calculating isotope similarity:
```

```
##   |                                                                                                                                   |                                                                                                                           |   0%  |                                                                                                                                   |===============================                                                                                            |  25%  |                                                                                                                                   |==============================================================                                                             |  50%  |                                                                                                                                   |============================================================================================                               |  75%  |                                                                                                                                   |===========================================================================================================================| 100%
```

```
## Calculating isotopic fidelity metrics:
```

```
##   |                                                                                                                                   |                                                                                                                           |   0%  |                                                                                                                                   |=                                                                                                                          |   1%  |                                                                                                                                   |==                                                                                                                         |   1%  |                                                                                                                                   |==                                                                                                                         |   2%  |                                                                                                                                   |===                                                                                                                        |   2%  |                                                                                                                                   |===                                                                                                                        |   3%  |                                                                                                                                   |====                                                                                                                       |   3%  |                                                                                                                                   |=====                                                                                                                      |   4%  |                                                                                                                                   |======                                                                                                                     |   5%  |                                                                                                                                   |=======                                                                                                                    |   5%  |                                                                                                                                   |=======                                                                                                                    |   6%  |                                                                                                                                   |========                                                                                                                   |   6%  |                                                                                                                                   |========                                                                                                                   |   7%  |                                                                                                                                   |=========                                                                                                                  |   7%  |                                                                                                                                   |==========                                                                                                                 |   8%  |                                                                                                                                   |===========                                                                                                                |   9%  |                                                                                                                                   |============                                                                                                               |   9%  |                                                                                                                                   |============                                                                                                               |  10%  |                                                                                                                                   |=============                                                                                                              |  10%  |                                                                                                                                   |=============                                                                                                              |  11%  |                                                                                                                                   |==============                                                                                                             |  11%  |                                                                                                                                   |===============                                                                                                            |  12%  |                                                                                                                                   |================                                                                                                           |  13%  |                                                                                                                                   |=================                                                                                                          |  14%  |                                                                                                                                   |==================                                                                                                         |  15%  |                                                                                                                                   |===================                                                                                                        |  15%  |                                                                                                                                   |===================                                                                                                        |  16%  |                                                                                                                                   |====================                                                                                                       |  16%  |                                                                                                                                   |====================                                                                                                       |  17%  |                                                                                                                                   |=====================                                                                                                      |  17%  |                                                                                                                                   |======================                                                                                                     |  18%  |                                                                                                                                   |=======================                                                                                                    |  19%  |                                                                                                                                   |========================                                                                                                   |  19%  |                                                                                                                                   |========================                                                                                                   |  20%  |                                                                                                                                   |=========================                                                                                                  |  20%  |                                                                                                                                   |=========================                                                                                                  |  21%  |                                                                                                                                   |==========================                                                                                                 |  21%  |                                                                                                                                   |===========================                                                                                                |  22%  |                                                                                                                                   |============================                                                                                               |  23%  |                                                                                                                                   |=============================                                                                                              |  23%  |                                                                                                                                   |=============================                                                                                              |  24%  |                                                                                                                                   |==============================                                                                                             |  24%  |                                                                                                                                   |==============================                                                                                             |  25%  |                                                                                                                                   |===============================                                                                                            |  25%  |                                                                                                                                   |================================                                                                                           |  26%  |                                                                                                                                   |=================================                                                                                          |  27%  |                                                                                                                                   |==================================                                                                                         |  27%  |                                                                                                                                   |==================================                                                                                         |  28%  |                                                                                                                                   |===================================                                                                                        |  28%  |                                                                                                                                   |===================================                                                                                        |  29%  |                                                                                                                                   |====================================                                                                                       |  29%  |                                                                                                                                   |=====================================                                                                                      |  30%  |                                                                                                                                   |======================================                                                                                     |  31%  |                                                                                                                                   |=======================================                                                                                    |  31%  |                                                                                                                                   |=======================================                                                                                    |  32%  |                                                                                                                                   |========================================                                                                                   |  32%  |                                                                                                                                   |========================================                                                                                   |  33%  |                                                                                                                                   |=========================================                                                                                  |  33%  |                                                                                                                                   |==========================================                                                                                 |  34%  |                                                                                                                                   |===========================================                                                                                |  35%  |                                                                                                                                   |============================================                                                                               |  35%  |                                                                                                                                   |============================================                                                                               |  36%  |                                                                                                                                   |=============================================                                                                              |  36%  |                                                                                                                                   |=============================================                                                                              |  37%  |                                                                                                                                   |==============================================                                                                             |  37%  |                                                                                                                                   |==============================================                                                                             |  38%  |                                                                                                                                   |===============================================                                                                            |  38%  |                                                                                                                                   |===============================================                                                                            |  39%  |                                                                                                                                   |================================================                                                                           |  39%  |                                                                                                                                   |=================================================                                                                          |  40%  |                                                                                                                                   |==================================================                                                                         |  41%  |                                                                                                                                   |===================================================                                                                        |  41%  |                                                                                                                                   |===================================================                                                                        |  42%  |                                                                                                                                   |====================================================                                                                       |  42%  |                                                                                                                                   |====================================================                                                                       |  43%  |                                                                                                                                   |=====================================================                                                                      |  43%  |                                                                                                                                   |======================================================                                                                     |  44%  |                                                                                                                                   |=======================================================                                                                    |  45%  |                                                                                                                                   |========================================================                                                                   |  45%  |                                                                                                                                   |========================================================                                                                   |  46%  |                                                                                                                                   |=========================================================                                                                  |  46%  |                                                                                                                                   |=========================================================                                                                  |  47%  |                                                                                                                                   |==========================================================                                                                 |  47%  |                                                                                                                                   |===========================================================                                                                |  48%  |                                                                                                                                   |============================================================                                                               |  49%  |                                                                                                                                   |=============================================================                                                              |  49%  |                                                                                                                                   |=============================================================                                                              |  50%  |                                                                                                                                   |==============================================================                                                             |  50%  |                                                                                                                                   |==============================================================                                                             |  51%  |                                                                                                                                   |===============================================================                                                            |  51%  |                                                                                                                                   |================================================================                                                           |  52%  |                                                                                                                                   |=================================================================                                                          |  53%  |                                                                                                                                   |==================================================================                                                         |  53%  |                                                                                                                                   |==================================================================                                                         |  54%  |                                                                                                                                   |===================================================================                                                        |  54%  |                                                                                                                                   |===================================================================                                                        |  55%  |                                                                                                                                   |====================================================================                                                       |  55%  |                                                                                                                                   |=====================================================================                                                      |  56%  |                                                                                                                                   |======================================================================                                                     |  57%  |                                                                                                                                   |=======================================================================                                                    |  57%  |                                                                                                                                   |=======================================================================                                                    |  58%  |                                                                                                                                   |========================================================================                                                   |  58%  |                                                                                                                                   |========================================================================                                                   |  59%  |                                                                                                                                   |=========================================================================                                                  |  59%  |                                                                                                                                   |==========================================================================                                                 |  60%  |                                                                                                                                   |===========================================================================                                                |  61%  |                                                                                                                                   |============================================================================                                               |  61%  |                                                                                                                                   |============================================================================                                               |  62%  |                                                                                                                                   |=============================================================================                                              |  62%  |                                                                                                                                   |=============================================================================                                              |  63%  |                                                                                                                                   |==============================================================================                                             |  63%  |                                                                                                                                   |==============================================================================                                             |  64%  |                                                                                                                                   |===============================================================================                                            |  64%  |                                                                                                                                   |===============================================================================                                            |  65%  |                                                                                                                                   |================================================================================                                           |  65%  |                                                                                                                                   |=================================================================================                                          |  66%  |                                                                                                                                   |==================================================================================                                         |  67%  |                                                                                                                                   |===================================================================================                                        |  67%  |                                                                                                                                   |===================================================================================                                        |  68%  |                                                                                                                                   |====================================================================================                                       |  68%  |                                                                                                                                   |====================================================================================                                       |  69%  |                                                                                                                                   |=====================================================================================                                      |  69%  |                                                                                                                                   |======================================================================================                                     |  70%  |                                                                                                                                   |=======================================================================================                                    |  71%  |                                                                                                                                   |========================================================================================                                   |  71%  |                                                                                                                                   |========================================================================================                                   |  72%  |                                                                                                                                   |=========================================================================================                                  |  72%  |                                                                                                                                   |=========================================================================================                                  |  73%  |                                                                                                                                   |==========================================================================================                                 |  73%  |                                                                                                                                   |===========================================================================================                                |  74%  |                                                                                                                                   |============================================================================================                               |  75%  |                                                                                                                                   |=============================================================================================                              |  75%  |                                                                                                                                   |=============================================================================================                              |  76%  |                                                                                                                                   |==============================================================================================                             |  76%  |                                                                                                                                   |==============================================================================================                             |  77%  |                                                                                                                                   |===============================================================================================                            |  77%  |                                                                                                                                   |================================================================================================                           |  78%  |                                                                                                                                   |=================================================================================================                          |  79%  |                                                                                                                                   |==================================================================================================                         |  79%  |                                                                                                                                   |==================================================================================================                         |  80%  |                                                                                                                                   |===================================================================================================                        |  80%  |                                                                                                                                   |===================================================================================================                        |  81%  |                                                                                                                                   |====================================================================================================                       |  81%  |                                                                                                                                   |=====================================================================================================                      |  82%  |                                                                                                                                   |======================================================================================================                     |  83%  |                                                                                                                                   |=======================================================================================================                    |  83%  |                                                                                                                                   |=======================================================================================================                    |  84%  |                                                                                                                                   |========================================================================================================                   |  84%  |                                                                                                                                   |========================================================================================================                   |  85%  |                                                                                                                                   |=========================================================================================================                  |  85%  |                                                                                                                                   |==========================================================================================================                 |  86%  |                                                                                                                                   |===========================================================================================================                |  87%  |                                                                                                                                   |============================================================================================================               |  88%  |                                                                                                                                   |=============================================================================================================              |  89%  |                                                                                                                                   |==============================================================================================================             |  89%  |                                                                                                                                   |==============================================================================================================             |  90%  |                                                                                                                                   |===============================================================================================================            |  90%  |                                                                                                                                   |===============================================================================================================            |  91%  |                                                                                                                                   |================================================================================================================           |  91%  |                                                                                                                                   |=================================================================================================================          |  92%  |                                                                                                                                   |==================================================================================================================         |  93%  |                                                                                                                                   |===================================================================================================================        |  93%  |                                                                                                                                   |===================================================================================================================        |  94%  |                                                                                                                                   |====================================================================================================================       |  94%  |                                                                                                                                   |====================================================================================================================       |  95%  |                                                                                                                                   |=====================================================================================================================      |  95%  |                                                                                                                                   |======================================================================================================================     |  96%  |                                                                                                                                   |=======================================================================================================================    |  97%  |                                                                                                                                   |========================================================================================================================   |  97%  |                                                                                                                                   |========================================================================================================================   |  98%  |                                                                                                                                   |=========================================================================================================================  |  98%  |                                                                                                                                   |=========================================================================================================================  |  99%  |                                                                                                                                   |========================================================================================================================== |  99%  |                                                                                                                                   |===========================================================================================================================| 100%
```

```
## Removing confirmed isotope entries:
```

```
## Recalculating peaklist for plotting:
```

```
##   |                                                                                                                                   |                                                                                                                           |   0%  |                                                                                                                                   |===============================                                                                                            |  25%  |                                                                                                                                   |==============================================================                                                             |  50%  |                                                                                                                                   |============================================================================================                               |  75%  |                                                                                                                                   |===========================================================================================================================| 100%
```

```
## Calculating adduct similarity:
```

```
##   |                                                                                                                                   |                                                                                                                           |   0%  |                                                                                                                                   |===============================                                                                                            |  25%  |                                                                                                                                   |==============================================================                                                             |  50%  |                                                                                                                                   |============================================================================================                               |  75%  |                                                                                                                                   |===========================================================================================================================| 100%
```

```
## 100 out of 186
```

```
##   |                                                                                                                                   |                                                                                                                           |   0%  |                                                                                                                                   |===============================                                                                                            |  25%  |                                                                                                                                   |==============================================================                                                             |  50%  |                                                                                                                                   |============================================================================================                               |  75%  |                                                                                                                                   |===========================================================================================================================| 100%
```

```r
saveRDS(negHermes, file = "D:/HermesResults/SurfaceWater/NegSOIFilter.rds")
```

## Quality control - Plotting



## Generating an Inclusion List
Once you are satisfied with a SOI list it is time to generate an inclusion list
to perform MSMS fragmentation of the SOI and proceed with the workflow.

Select the desired SOI list and the IL parameters:

 * **mode**: We offer three different inclusion list modes which offer some 
flexibility when designing the number of MSMS injections. "Full" will use all
SOIs, "only" will just pick the SOIs that are associated with a list of adducts
ad and "priorizate" will, if possible, only pick the first adduct in ad for a
certain formula, then the next, and so on. This last mode avoids fragmenting 
redundant adduct signals.
 * **ad**: The list of adducts to pick/prioritize. Check your 
object@metadata@ExpParam@adlist to ensure proper spelling.
 * **filtermz**: **Really** important, is the mz tolerance of your quadrupole filter.
 It has to be properly set so that SOIs overlapping in a mz-rt neighborhood are 
 grouped together into a single IL entry. Feel free to experiment, but larger mz
 values can yield convoluted the MSMS spectra. 


```r
myHermes <- genIL(myHermes, id = 3, par = ILParam(priorization = "full"))
```

```
## Now processing the IL:
```

```
## Grouping the entries
```

```
## Tidying the entries
```

```
## Done!
```

```r
myHermes <- genIL(myHermes, id = 3, par = ILParam(priorization = "yes", ad = c("M+H", "M+NH4", "M+CH3OH+H")))
```

```
## Now processing the IL:
```

```
## Prioritizing M+H
```

```
## Prioritizing M+NH4
```

```
## Prioritizing M+CH3OH+H
```

```
## Grouping the entries
```

```
## Tidying the entries
```

```
## Done!
```

```r
myHermes <- genIL(myHermes, id = 6, par = ILParam(priorization = "full"))
```

```
## Now processing the IL:
```

```
## Grouping the entries
```

```
## Tidying the entries
```

```
## Done!
```

```r
myHermes <- genIL(myHermes, id = 6, par = ILParam(priorization = "yes", ad = c("M+H", "M+NH4", "M+CH3OH+H")))
```

```
## Now processing the IL:
```

```
## Prioritizing M+H
```

```
## Prioritizing M+NH4
```

```
## Prioritizing M+CH3OH+H
```

```
## Grouping the entries
```

```
## Tidying the entries
```

```
## Done!
```

```r
myHermes <- genIL(myHermes, id = 4, par = ILParam(priorization = "full"))
```

```
## Now processing the IL:
```

```
## Grouping the entries
```

```
## Tidying the entries
```

```
## Done!
```

```r
myHermes <- genIL(myHermes, id = 4, par = ILParam(priorization = "yes", ad = c("M+H", "M+NH4", "M+CH3OH+H")))
```

```
## Now processing the IL:
```

```
## Prioritizing M+H
```

```
## Prioritizing M+NH4
```

```
## Prioritizing M+CH3OH+H
```

```
## Grouping the entries
```

```
## Tidying the entries
```

```
## Done!
```

```r
saveRDS(myHermes, file = "D:/HermesResults/SurfaceWater/PosSOIFilter.rds")


negHermes <- genIL(negHermes, id = 3, par = ILParam(priorization = "full"))
```

```
## Now processing the IL:
```

```
## Grouping the entries
```

```
## Tidying the entries
```

```
## Done!
```

```r
negHermes <- genIL(negHermes, id = 3, par = ILParam(priorization = "yes", ad = c("M-H")))
```

```
## Now processing the IL:
```

```
## Prioritizing M-H
```

```
## Grouping the entries
```

```
## Tidying the entries
```

```
## Done!
```

```r
negHermes <- genIL(negHermes, id = 6, par = ILParam(priorization = "full"))
```

```
## Now processing the IL:
```

```
## Grouping the entries
```

```
## Tidying the entries
```

```
## Done!
```

```r
negHermes <- genIL(negHermes, id = 6, par = ILParam(priorization = "yes", ad = c("M-H")))
```

```
## Now processing the IL:
```

```
## Prioritizing M-H
```

```
## Grouping the entries
```

```
## Tidying the entries
```

```
## Done!
```

```r
negHermes <- genIL(negHermes, id = 4, par = ILParam(priorization = "full"))
```

```
## Now processing the IL:
```

```
## Grouping the entries
```

```
## Tidying the entries
```

```
## Done!
```

```r
negHermes <- genIL(negHermes, id = 4, par = ILParam(priorization = "yes", ad = c("M-H")))
```

```
## Now processing the IL:
```

```
## Prioritizing M-H
```

```
## Grouping the entries
```

```
## Tidying the entries
```

```
## Done!
```

```r
saveRDS(negHermes, file = "D:/HermesResults/SurfaceWater/NegSOIFilter.rds")
```

Creating an inclusion list will generate an instance of an MS2Exp, a container
which will serve to store your IL, the MSMS data you will later acquire and the
identifications that arise from said data.

With the inclusion list ready you can just export it to a csv format that can be
used as input for the MSMS experiment. Just select the IL to export, where to 
save it, the maximum number of ions to monitor at the same time and whether
you want to get separate csv files or just one csv.


```r
# exportIL(myHermes, 1, folder = getwd(), maxOver = 10, sepFiles = TRUE)
# exportIL(negHermes, 1, folder = getwd(), maxOver = 10, sepFiles = TRUE)
```

# MS2 Processing

After performing the MSMS runs with the exported IL you can process the
generated data with RHermes and obtain compound identifications.

## Using MS2Proc


```r
ms2f <- list.files("D:/ABrunner/MSMSdata", pattern = ".*pos.*.mzML", full.names = TRUE)[seq(1,31)]
myHermes <- MS2Proc(myHermes, 1, MS2files = ms2f, sstype = "regular", referenceDB = "D:/MS2ID_B2R_20201113_083214.rds", useDB = TRUE, mincos = 0.5)
myHermes <- MS2Proc(myHermes, 2, MS2files = ms2f,sstype = "regular",  referenceDB = "D:/MS2ID_B2R_20201113_083214.rds", useDB = TRUE, mincos = 0.5)


ms2f <- list.files("D:/ABrunner/MSMSdata", pattern = ".*pos.*.mzML", full.names = TRUE)[seq(33,45)]
myHermes <- MS2Proc(myHermes, 3, MS2files = ms2f,sstype = "regular",  referenceDB = "D:/MS2ID_B2R_20201113_083214.rds", useDB = TRUE, mincos = 0.5)
myHermes <- MS2Proc(myHermes, 4, MS2files = ms2f,sstype = "regular",  referenceDB = "D:/MS2ID_B2R_20201113_083214.rds", useDB = TRUE, mincos = 0.5)

ms2f <- list.files("D:/ABrunner/MSMSdata", pattern = ".*pos.*.mzML", full.names = TRUE)[c(1:31, seq(33,45))]
myHermes <- MS2Proc(myHermes, 5, MS2files = ms2f,sstype = "regular",  referenceDB = "D:/MS2ID_B2R_20201113_083214.rds", useDB = TRUE, mincos = 0.5)
myHermes <- MS2Proc(myHermes, 6, MS2files = ms2f,sstype = "regular",  referenceDB = "D:/MS2ID_B2R_20201113_083214.rds", useDB = TRUE, mincos = 0.5)

saveRDS(myHermes, file = "D:/HermesResults/SurfaceWater/MS1_pos_Ident_Water.rds")
```




```r
ms2f <- list.files("D:/ABrunner/MSMSdata", pattern = ".*neg.*.mzML", full.names = TRUE)[seq(3,33)]
negHermes <- MS2Proc(negHermes, 1, MS2files = ms2f,sstype = "regular",  referenceDB = "D:/MS2ID_B2R_20201113_083214.rds", useDB = TRUE, mincos = 0.5)
negHermes <- MS2Proc(negHermes, 2, MS2files = ms2f,sstype = "regular",  referenceDB = "D:/MS2ID_B2R_20201113_083214.rds", useDB = TRUE, mincos = 0.5)


ms2f <- list.files("D:/ABrunner/MSMSdata", pattern = ".*neg.*.mzML", full.names = TRUE)[seq(35,45)]
negHermes <- MS2Proc(negHermes, 3, MS2files = ms2f,sstype = "regular",  referenceDB = "D:/MS2ID_B2R_20201113_083214.rds", useDB = TRUE, mincos = 0.5)
negHermes <- MS2Proc(negHermes, 4, MS2files = ms2f,sstype = "regular",  referenceDB = "D:/MS2ID_B2R_20201113_083214.rds", useDB = TRUE, mincos = 0.5)
saveRDS(negHermes, file = "D:/HermesResults/SurfaceWater/MS1_neg_Ident_Water.rds")
```


```r
sink("D:/HermesResults/SurfaceWater/sessionInfo.txt")
sessionInfo()
```

```
## R version 3.6.3 (2020-02-29)
## Platform: x86_64-w64-mingw32/x64 (64-bit)
## Running under: Windows 10 x64 (build 18363)
## 
## Matrix products: default
## 
## locale:
## [1] LC_COLLATE=Spanish_Spain.1252  LC_CTYPE=Spanish_Spain.1252    LC_MONETARY=Spanish_Spain.1252 LC_NUMERIC=C                  
## [5] LC_TIME=Spanish_Spain.1252    
## 
## attached base packages:
## [1] stats     graphics  grDevices utils     datasets  methods   base     
## 
## other attached packages:
##  [1] BiocParallel_1.20.1 keras_2.2.5.0       data.table_1.12.8   RCy3_2.6.3          gtools_3.8.1        forcats_0.4.0      
##  [7] stringr_1.4.0       dplyr_0.8.4         purrr_0.3.3         readr_1.3.1         tidyr_1.0.2         tibble_3.0.1       
## [13] tidyverse_1.3.0     viridis_0.5.1       viridisLite_0.3.0   ggplot2_3.3.0       CHNOSZ_1.3.4        magrittr_1.5       
## [19] mzR_2.20.0          Rcpp_1.0.4.6        igraph_1.2.4.2      enviPat_2.4         RHermes_0.99.9007  
## 
## loaded via a namespace (and not attached):
##  [1] colorspace_1.4-1     ellipsis_0.3.0       XVector_0.26.0       base64enc_0.1-3      fs_1.3.1             rstudioapi_0.11     
##  [7] DT_0.13              fansi_0.4.1          lubridate_1.7.4      xml2_1.2.2           codetools_0.2-16     R.methodsS3_1.7.1   
## [13] ncdf4_1.17           knitr_1.28           zeallot_0.1.0        jsonlite_1.6.1       packrat_0.5.0        broom_0.5.6         
## [19] dbplyr_1.4.2         png_0.1-7            R.oo_1.23.0          tfruns_1.4           shinydashboard_0.7.1 graph_1.64.0        
## [25] shiny_1.5.0          compiler_3.6.3       httr_1.4.1           backports_1.1.7      assertthat_0.2.1     fastmap_1.0.1       
## [31] lazyeval_0.2.2       cli_2.0.2            later_1.0.0          visNetwork_2.0.9     htmltools_0.5.0      tools_3.6.3         
## [37] gtable_0.3.0         glue_1.4.0           rappdirs_0.3.1       Biobase_2.46.0       cellranger_1.1.0     slickR_0.4.9        
## [43] vctrs_0.2.4          Biostrings_2.54.0    RJSONIO_1.3-1.4      nlme_3.1-144         iterators_1.0.12     xfun_0.12           
## [49] rvest_0.3.5          mime_0.9             lifecycle_0.2.0      XML_3.99-0.3         zlibbioc_1.32.0      scales_1.1.0        
## [55] hms_0.5.3            doSNOW_1.0.18        promises_1.1.0       ProtGenerics_1.18.0  parallel_3.6.3       yaml_2.2.1          
## [61] reticulate_1.14      gridExtra_2.3        stringi_1.4.5        S4Vectors_0.24.3     tensorflow_2.0.0     foreach_1.4.7       
## [67] BiocGenerics_0.32.0  rlang_0.4.6          pkgconfig_2.0.3      evaluate_0.14        lattice_0.20-38      htmlwidgets_1.5.1   
## [73] tidyselect_1.0.0     R6_2.4.1             IRanges_2.20.2       snow_0.4-3           generics_0.0.2       DBI_1.1.0           
## [79] pillar_1.4.3         haven_2.2.0          whisker_0.4          withr_2.2.0          KEGGREST_1.26.1      modelr_0.1.5        
## [85] crayon_1.3.4         shinyWidgets_0.5.3   plotly_4.9.2.1       rmarkdown_2.1        grid_3.6.3           readxl_1.3.1        
## [91] reprex_0.3.0         digest_0.6.25        xtable_1.8-4         httpuv_1.5.2         R.utils_2.9.2        stats4_3.6.3        
## [97] munsell_0.5.0
```

```r
paste("Date: ", date())
```

```
## [1] "Date:  Fri Jan 01 08:10:40 2021"
```

```r
sink()
rm(list = ls())
```




