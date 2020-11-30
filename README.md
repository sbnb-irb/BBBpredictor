In order to create the image and run the container:

### Go to the singularity folder
cd singularity

#### Build the image
sudo singularity build BBBpredictor.simg BBBpredictor.def

### You're ready to run the container as any executable on the example input file provided in the singularity/TEST folder

./BBBpredictor.simg TEST/all_smiles.txt

### Usage:

Blood-brain barrier (BBB) predictor:
------------------------------------
Classifier that will predict if a chemical coumpound will pass the Blood-brain barrier.
SBNB lab (IRB Barcelona) - Nov 2020

USAGE: BBBpredictor myfile.csv

Where my myfile.csv is an input text file that contains two columns separated by commas or spaces/tabs
COLUMN 1: compound id (ex:1)
COLUMN 2: compound SMILES string

NOTE: Lines starting with '#' will be ignored

OUTPUT: csv file containing the prediction for each compound

PREDICTION LEGEND: passes the BBB?
0: no
1: yes
-1: the molecular signature could not be calculated for this coumpound
