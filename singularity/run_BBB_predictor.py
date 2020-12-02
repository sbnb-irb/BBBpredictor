import sys, os
import joblib
import numpy as np

usage="""
Blood-brain barrier (BBB) predictor:
------------------------------------
Classifier that will predict if a chemical coumpound will pass the Blood-brain barrier.
SBNB lab (IRB Barcelona) - Nov 2020.

USAGE: BBBpredictor myfile.tsv [SMILES or INCHI]

Where my myfile.tsv is an input text file that contains two columns separated by TABS.
COLUMN 1: compound unique id (ex:1 or a molecule name without tab characters inside).
COLUMN 2: compound SMILES string.

NOTE: The default format is SMILES and can be turned into InChI if the second argument INCHI is provided.
NOTE: Lines starting with '#' as well as empty lines are ignored.

OUTPUT: tsv file containing the prediction for each compound

PREDICTION LEGEND: passes the BBB?
0: no
1: yes
-1: the molecular signature could not be calculated for this coumpound
"""

def formatArray(listofLists):
    # if np.nan is present, replace by -1
    listofLists = [l if not np.any(np.isnan(l)) else [-1 for x in l] for l in listofLists]
    return np.array(list([list(l) for l in listofLists]))


if len(sys.argv) <2:
    print(usage)
    sys.exit(1)


# Martino's signaturizer
from signaturizer import Signaturizer

inputFile=sys.argv[1]
entries=dict()
outputDic= dict()

formatMol= "SMILES"
if len(sys.argv)>2 and sys.argv[2].upper() == 'INCHI':
    formatMol='INCHI'

print("INPUT FORMAT: {}\n".format(formatMol))

rootModels="/opt/NNmodels"    # in simg
#rootModels="NNmodels"        # test

path_to_model="/opt/rf_from_signZ_paper_full.joblib"   # in simg
#path_to_model="../MLmodels/rf_from_signZ_paper_full.joblib"    # test

listModels=[os.path.join(rootModels,x+y) for x in "A B C D E".split() for y in "1 2 3 4 5".split()]


RFmodel= joblib.load(path_to_model)

# Parsing the input file
with open(inputFile) as f:
    try:
        lines=f.readlines()
        for l in lines:
            if len(l.strip())>0 and not l.startswith('#'):
                li= l.split('\t')

                entries[li[0].strip()] = li[1].strip()

    except Exception as e:
        print("ERROR")
        print(e)
        print("Quitting, please check your input file!\n")
        print(usage)
        sys.exit(1)

# Output dict
outputDic= dict()

# load the predictor for ALL spaces (the global signature)
print("Initializing the signaturizer, please wait..")
sign = Signaturizer(model_name=listModels, local=True)


# Submit all molecules to the signaturizer
smiles = list(entries.values())
clefs= list(entries.keys())

print("\n INPUT:")
print(smiles)

# Predict the signatures
print("\nCalculating the signatures for the {:d} input molecules:\n".format(len(entries)))
signZlist = sign.predict(smiles, keytype=formatMol).signature
print("SHAPE of the returned signatures",signZlist.shape,'\n')

assert (len(signZlist) == len(clefs)),"ERROR: The signature list and compound id list don't have the same length!"


# dict compoundid: signZ 
for i in range(len(clefs)):
    outputDic[clefs[i]] = signZlist[i]

X=formatArray(signZlist)

y_hat = RFmodel.predict(X)

y_final=[]

# Check whether the signature is all -1 (there was a problem)
problems =[]
for i in range(X.shape[0]):
    if np.all(X[i] == -1):
        problems.append(clefs[i])
        y_final.append(-1)
    else:
        y_final.append(y_hat[i])

if len(problems)>0:
    print("WARNING: no signature for the following {:d} coumpounds".format(len(problems)))
    print(problems)

outFile= os.path.basename(inputFile).split('.')[0]+"_BBBpredictions.tsv"
with open(outFile,'w') as f:
    f.write("#Brain blood barrier predictions\n")
    f.write("#1: pass, 0: does not pass, -1: error\n")
    f.write("#Compound\tBBBpred\n\n")
    for i in range(len(clefs)):
        f.write(str(clefs[i])+"\t"+str(y_final[i])+'\n')

print("Predictions written in", outFile)


# print the result to std output
# print("#Compound id\tprediction\n")
# for k,v in outputDic.items():
 #   print(k,'\t',v)

