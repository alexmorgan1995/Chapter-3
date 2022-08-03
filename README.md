# Chapter 3 Model - Understanding the role of global food trade on the transmission dynamics of antibiotic-resistant foodborne bacteria

## Foodborne AMR Model

This repository contains the files used in the third data chapter of the PhD (Chapter 4). This chapter looks at exploring the effect of food product importation on the transmission dynamics of AMR occuring within a domestic country of interest. This builds off the model developed in the second data chapter of the PhD (Chapter 3). This was conducted for one main case study the UK and *Salmonella* spp. in fattening pigs. 

This chapter was split into two parts, with both a "simple" and a "complex" model being created exploring less/more heterogeneity in food product importation.

**As of [03/21] all finalised models and data can be found in this repository**

## Folder Structure

All data used for both sections of the chapter and for the ABC-SMC model fitting procedure can be found in the ```/Model_Fit_Data``` folder. Part 1 analyses and model fitting can be found in the ```/Part_1_Basic``` folder, while Part 2 (complex heterogenous model) can be found in the ```/Part_2_Complex``` folder. 

### Section 1 Basic Model Folder

The main files for the analysis and the model fitting proceduce can be found in the main folder. This includes descriptive analyses, uncertainty analyses, trajectory plots, and model fitting files. There are also a number of supplementary analysis, found in the main folder and also in the ```/no_import_fit``` folder. This includes a model fit where import (psi) was set to zero and compared against the main text model. There is also a deprecated folder containing older versions of the R files, ```/Deprecated```.

### Section 2 Complex Model Folder

The main files for the analysis and the model fitting proceduce can be found in the main folder. This includes descriptive analyses, uncertainty analyses, trajectory plots, and model fitting files. There are also a number of supplementary analysis, found in the ```/Supplementary``` folder. This includes a model fit where import (psi) was set to zero and compared against the main text model. There is also a deprecated folder containing older versions of the R files, ```/Deprecated```.

## Model Details

All code was written and run in R-Studio. The packages required to run the code are as follows:

```library("deSolve"); library("fast"); library("sensitivity"); library("ggplot2"); library("plotly"); library("reshape2")```

deSolve is the primary package used to solve the ODEs used in the modelling approach.

