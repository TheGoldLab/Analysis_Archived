README file for code associated with:

A bias-variance trade-off governs individual differences in on-line learning in an unpredictable environment

Christopher M. Glaze, Alexandre L.S. Filipowicz, Joseph W. Kable, Vijay Balasubramanian, and Joshua I. Gold

Nature Human Behavior
2018

DEPENDENCIES
————————————
1. Tested on Matlab 9.0.0.341360 (R2016a)
2. Uses some Gold Lab Utilities available here: https://github.com/TheGoldLab/Lab-Matlab-Utilities

CONTENTS
————————
1. Fig_*.m … code to generate figures in the paper, numbered accordingly
2. getDataInfo.m … Utility function to point to where you put the data (see: “SET DATA DIRECTORY HERE”) 
3. plotModelFreeAnalyses.m … utility used by some of the figure-generating functions
4. Models/ … c/mex files to implement the models from the paper. See the README in that directory for more details.

TO GENERATE FIGURES
———————————————————
1. Put everything on your path
2. Get the data by contacting Joshua Gold at jigold@pennmedicine.upenn.edu
3. Update the file ‘getDataInfo.m’ (see: “SET DATA DIRECTORY HERE”) with references to the location of the installed data directories (they will be provided as ‘Data/Analysis’ and ‘Data/Raw’)
4. Run the files: “Fig_*.m”

