# SIGMA

# The Structure-Informed Genetic Missense mutation Assessor

This repository provides the code and data for the paper "SIGMA leverages protein structural information to predict the pathogenicity of missense variants".

### BACKGROUND

The task of assessing the pathogenicity of missense variants is key to interpreting genetic data. A promising strategy involves evaluating variant effects within the context of protein structure. However, the scarcity of known 3D protein structures has hindered the exploitation of structural information, as protein structure is critically important for ensuring proper molecular function. The advent of AlphaFold2 has begun to change this situation by enabling extensive, accurate predictions of 3D structures. Here, we introduce SIGMA, which utilizes AlphaFold2 predictions to evaluate the effects of missense variants in the context of predicted protein structures.

### What is SIGMA

SIGMA is an advanced tool designed to predict the pathogenicity of missense variants in proteins by utilizing structural information from AlphaFold2 predictions, and it demonstrates superior performance in comparison to other existing predictors.

<img src="https://github.com/zhq921/SIGMA/blob/main/imgs/workflow.png" width = "50%" />

### Online tool for the SIGMA model

To facilitate the application of SIGMA, we pre-computed SIGMA scores for over 48 million possible missense variants across 3,454 disease-associated genes and developed an interactive online platform (https://www.sigma-pred.org/).

<img src="https://github.com/zhq921/SIGMA/blob/main/imgs/SIGMA_prediction.png" width = "70%" />

### Who do I talk to?

- Hengqiang Zhao (zhaohq921@gmail.com)
