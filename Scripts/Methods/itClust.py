
import ItClust as ic
import scanpy as sc
import os
import scanpy.external as sce
from numpy.random import seed
from tensorflow.compat.v1 import set_random_seed
import pandas as pd
import numpy as np
import warnings
import sys


##################################################

def get_data(folder, x):
	data = sc.read_csv((folder+"/data_"+ x +".csv"), delimiter=",", first_column_names=False)
	with open(folder + "/cells_" + x +".csv") as f:
		lines = f.read().splitlines()
	data.obs_names = lines
	return data 

print("Starting...")

os.environ["CUDA_VISIBLE_DEVICES"]="1"
warnings.filterwarnings("ignore")

args = sys.argv

trainFolder = args[1]
testFolder = args[2] 

print(args[3])
os.chdir(args[3])

#Set seeds
seed(20180806)
np.random.seed(10)
set_random_seed (20180806) # on GPU may be some other default

# Source data
print("Get source data...")

adata_train = get_data(trainFolder, "train")

print("Get source meta...")
celltypes = pd.read_csv((trainFolder + "/meta_train.csv"), sep=',')

adata_train.obs["celltype"] = celltypes.class_.values
adata_train.obs_names_make_unique(join="-")
adata_train.var_names_make_unique(join="-")
print(adata_train)
print(adata_train.obs.shape)


print("Get query data...")
adata_test =   get_data(testFolder, "test")

adata_test.var_names = adata_train.var_names
adata_test.obs_names_make_unique(join="-")
adata_test.var_names_make_unique(join="-")

print(adata_test)
print(adata_test.obs.shape)

print(" Fit ItClust model...")
clf= ic.transfer_learning_clf()

print("..........................................................................................")
print(adata_train.X)
print("---------")
print(adata_test.X)
print("..........................................................................................")
clf.fit(adata_train, adata_test)

print("Prediction...")
pred, prob, cell_type_pred=clf.predict()
pred.head()
print("DONE")