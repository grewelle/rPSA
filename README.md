# rPSA
Software package for the Revised PSA as defined by the manuscript: "Redefining Risk in Data-poor Fisheries", Grewelle et al. 2020

Implementation of rPSA software requires a Python IDE such as Pycharm https://www.jetbrains.com/pycharm/.  
An additional recommended requirement is installation of Anaconda https://docs.anaconda.com/anaconda/install/.

Download the Python file PSA_software.py and two template spreadsheets: productivity.csv and susceptibility.csv.

Open PSA_software.py and edit lines 137 and 138 to include the paths of the two downloaded spreadsheets.

Open each of the spreadsheets and edit according to your desired specifications while following these rules:
1) 1st four rows must be kept and edited as:
    - row 1: second column is user input for type of model used for productivity/susceptibility.  Type either 'additive' or 'multiplicative'.
    - row 2: second and third columns are user input for the cut-offs used in scoring attributes (e.g. 1=0-25%, 2=25-75%, 3=75-100%).  In the example here, the two cut-offs are 25% and 75% and should be represented as '1.0/4.0' and '3.0/4.0', respectively.  Default is 33% and 67%.
    - row 3: second column is user input for the maximum number of attributes used in the analysis for productivity or susceptibility.  Some species/stock are not evaluated with the maximum number of attributes. This is taken into account automatically in the code.
    - row 4: second and third columns are user input for the thresholds used in categorizing vulnerability scores as low, medium, high. Default is equal (i.e. 1/3rd each for low, medium, high).  Formatting is same as row 2.
    
2) 5th row is blank
3) 6th row is a header labeling columns of attribute values and weights (label as you wish)
4) 7th row onwards are attribute values and weights.  If attribute is not used for the species, leave cell blank.  For each species, the number of attributes used is calculated by the program based on the number of filled columns.
5) For the 7th row onwards:
    - column 1: enables higher level labeling
    - column 2: species/stock labeling
6) Leave two blank columns between attribute values and associated weights
7) Weight columns must be in the same order as the attribute values (from left to right)
