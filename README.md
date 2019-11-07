# Chapter-3
## Model code for Chapter 3 - Livestock Heterogeneity 

This chapter looks at exploring the effect of both livestock antibiotic usage heterogeneity and differential levels of food-product import at the country level . This is in relation to a specific country of interest - simplified in the following code using the UK as an example.

## Code Details

All scripts were written using R-Studio. The base packages used to run this code are:

`library("deSolve"); library("fast"); library("sensitivity"); library("ggplot2"); library("plotly"); library("reshape2")`

There are also multiple complexities of AMR models used in this chapter. These include:

- AMR model with Domestic and Imported Livestock Food Products
  - ```Chapter 3 Model v3 - very simplified - Imp vs Dom.R```
- AMR model with Domestic, EU and non-EU imported Food Products - "Simplified"
- AMR model with Domestic and 7 Importing Countries - "Expanded"
  - Denmark
  - Germany
  - Netherlands
  - Ireland
  - Poland
  - Spain
  - Belgium
