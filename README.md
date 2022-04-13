Code used to generate results, figures and tables for DOI: .

The data used to generate this publication were obtained from the Neotoma Paleoecology Database (NPD). These are not provided alongside this repository, but are freely available. We have included code to access NPD data as used in this project in the "data_acquisition_processing.R" file.

Three R scripts are provided in this repository, along with summarised data used for for the statistical models and analyses in the manuscript. These scripts makes extensive use of the code folding functionality in RStudio. Alt + O is the default shortcut on Windows and Linux versions of RStudio to collapse all folds.

data_acquisition_processing.R is a script to access NPD data via the API and clean pollen taxonomy using up-to-date data from The Plant List.

reproduce_analyses.R is a minimum working script that contains all of the data and analyses to produce the figures and tables found in the main text, and many of the supplementary analyses. This, this script makes use of files in respository sub-folders, so please point your working directory at the start of this script to its folder path, and include all other directories in the repository.

full_analyses.R shows the complete data analysis pathway (from the output of data_acquisition_processing.R), as well as main and supplementary analyses. This script requires respository sub-folders etc as per the analysis script above.
