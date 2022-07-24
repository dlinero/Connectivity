## Mapping priority sites to improve protected areas connectivity and bird conservation in Colombia
##### Daniela Linero, MSc - [National Audubon Society](https://www.audubon.org/) 
##### Jorge Velásquez, Ph.D. - [National Audubon Society](https://www.audubon.org/) 
##### Camilo Correa-Ayram, Ph.D. - [Departamento de Ecología y Territorio, Universidad Javeriana](https://fear.javeriana.edu.co/departamento-de-ecologia-y-territorio) 
\

### Project abstract
***
Colombia's ecosystems are home to the highest richness of bird species worldwide, yet they are being lost rapidly due to the growing threats of climate and land-use change. The National Audubon Society seeks to counteract the effects of these threats by supporting the establishment of a representative and ecologically connected system of protected areas in Colombia. The purpose of this project is to identify, for the first time on a national scale, the most critical places for conserving and restoring the connectivity among Colombian protected areas. We will harness the power of Azure cloud computing to build high-resolution species distribution models for ten priority birds based on open-source occurrence data and environmental GIS layers. With the resulting maps, we will apply least-cost and circuit theory algorithms to model large-scale connectivity priorities among national and sub-national protected areas. 


### GitHub repository structure
***


```
¦
+---data
¦   +---   01_cleanOccurrences                     : Species occurrences obtained from eBird
¦   +---   02_SDMs                                 : Raw and processed environmental data
¦   +---   species_distributions                   : Shapefiles of the species distributions
¦   
+---outputs
¦    +---   01_cleanOccurrences                    : Clean occurrences
¦    +---   02_SDMs                                : Results of species distribution models under different frameworks
¦
+---scripts  
¦    +---   01_cleanOccurrences.R                  : R script to clean occurrences
¦    +---   01_cleanOccurrences_markdown.html      : Pdf file describing cleaning procedure
¦    +---   01_cleanOccurrences_markdown.Rmd       : Markdown file describing cleaning procedure
¦    +---   01_cleanOccurrences_markdown.tex       : Latex file associated with markdown
¦    +---   02_Build_sampling_probability_map.R    : R script to build sampling probability map based on birds occurrences
¦    +---   02_Prepare_ESA_data.R                  : R script to prepare land cover data for SDMs
¦    +---   02_SDMs_correcting_samplingBias.R      : SDM models under different frameworks to correct for sampling bias                        
¦    +---   02_SDMs_correcting_samplingBias_loop.R : Loop that accelerates SDM models to correct for sampling bias
¦    +---   02_SDMs_models.R                       : SDM models following the simplest steps of Wallace
¦    +---   02_SDMs_models_ESA.R                   : SDM models incorporating land cover data
¦    +---   02_SDMs_thinning.R                     : Script to do spatial thinning of occurrences
¦    +---   01_cleanOccurrences.R                  : Script to clean occurrences
¦    +---   03_composite_corridor_maps.R           : Reclassify corridors for each species to create composite maps
¦    +---   03_creation_resistance_surfaces.R      : Create movement resistance maps for each focal species
¦    +---   03_Selecting_cores_&_study_areas.R     : Selecting protected areas to connect for each species

¦
¦
¦   README.md                                      : Description of the repository
¦   Connectivity.Rproj                             : RStudio project file 


```


### Contact
***

Feel free to email me at daniela.linero@audubon.org 
