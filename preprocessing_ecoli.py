"""This script extract data from the original file and store the data in a
pickle file
"""
import csv
import numpy as np
import pickle
from python.toxicogenomic.utils.is_number import is_number

# Load data from the dataset file
file_dataset = open("/home/changyale/dataset/toxicogenomic/"+\
        "EcoliData_3Dataset_WithoutMovAvg.csv","rb")
reader = csv.reader(file_dataset)
lines = [line for line in reader]
file_dataset.close()

mtr_dataset = np.array(lines)

# Number of time points
n_times = 25
assert (mtr_dataset.shape[1]-1)%n_times == 0

# Number of genes used
n_genes = (mtr_dataset.shape[1]-1)/n_times

# Number of conditions
n_conditions = mtr_dataset.shape[0]

# Store the name of chemicals
conditions_name = []

# Construct a 3-d array to store the dataset
# dim 1--> conditions; dim 2--> genes; dim 3--> times;
data = np.zeros((n_conditions,n_genes,n_times))
for i in range(n_conditions):
    for j in range(n_genes):
        for k in range(n_times):
            if is_number(mtr_dataset[i,1+j*n_times+k]) == False:
                print([i,1+j*n_times+k,mtr_dataset[i,1+j*n_times+k]])
            #data[i,j,k] = float(mtr_dataset[i,1+j*n_times+k])
    conditions_name.append(mtr_dataset[i,0])


