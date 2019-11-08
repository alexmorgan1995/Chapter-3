# Chapter-3
## Model code for Chapter 3 - Livestock Heterogeneity 

This chapter looks at exploring the effect of both livestock antibiotic usage heterogeneity and differential levels of food-product import at the country level . This is in relation to a specific country of interest - simplified in the following code using the UK as an example.

The most current iteration for the simplified model is: ```Chapter 3 Model v3 - very simplified -  Imp vs Dom.R```

This the model I have been most recently working on. This models a simple domestic/imported livestock population and food product fraction which causes exposure for the human population in the domestic country of interest. This model was created to create a very general overview on the model dynamics of spatially heterogenous livestock populations, livestock antibiotic usage heterogeneity and the impact of food usage from differing sources. 

**There are a number of issues with the current iteration, essentially that the sensitivity analysis is currently not working and needs to be fixed.**

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
