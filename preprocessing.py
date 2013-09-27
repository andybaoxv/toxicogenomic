import csv
import numpy as np
import pickle

# Load data from the dataset file
file_dataset = open("12 DBPs and other chemicals.csv","rb")
reader = csv.reader(file_dataset)
lines = [line for line in reader]
file_dataset.close()

mtr_dataset = np.array(lines)

## Define the dimensions of the data matrix
# Number of Genes used
n_genes = 39

# Number of time points in collecting Gene expression data
n_times = 25
assert mtr_dataset.shape[0]%(n_genes+1) == 0
assert (mtr_dataset.shape[1]-6)%(n_times+1) == 0

# Number of toxicants used
n_toxicants = mtr_dataset.shape[0]/(n_genes+1)

# Number of dose concentration values used
n_concentrations = (mtr_dataset.shape[1]-1)/(n_times+1)

# Construct a 4-dim array to store the dataset
# dim 1--> toxicants; dim 2--> concentrations; dim 3--> genes; dim 4-->times;
data = np.zeros((n_toxicants,n_concentrations,n_genes,n_times))

# Extract data from the original dataset
# Specifically to this dataset "12 DBPs and other chemicals.csv"
# Toxicants 0-10 and 16-20 have the same delimiter(taking one cell)
# Toxicants 11-16 have the same delimiter(taking two cells)
for index_toxicants in list(set(range(n_toxicants))-set(range(11,16))):
    for index_concentrations in range(n_concentrations):
        for i in range(n_genes):
            for j in range(n_times):
                data[index_toxicants,index_concentrations,i,j] = float(\
                        mtr_dataset[1+index_toxicants*(n_genes+1)+i,\
                        2+index_concentrations*(n_times+1)+j])

for index_toxicants in range(11,16):
    for index_concentrations in range(n_concentrations):
        for i in range(n_genes):
            for j in range(n_times):
                data[index_toxicants,index_concentrations,i,j] = float(\
                        mtr_dataset[1+index_toxicants*(n_genes+1)+i,\
                        2+index_concentrations*(n_times+2)+j])

# Extract Gene and pathway information
# genes_name contains 39 names of different genes
# pathways_type contains corresponding pathway type for the genes
genes_name = []
pathways_type = []
for i in range(1,1+n_genes):
    genes_name.append(mtr_dataset[i,1])
    if mtr_dataset[i,0] != '':
        pathways_type.append(mtr_dataset[i,0])
    else:
        pathways_type.append(pathways_type[-1])

# Extract Toxicants information
toxicants_name = range(n_toxicants)
for i in range(n_toxicants):
    if i in range(0,7)+range(8,11):
        toxicants_name[i] = mtr_dataset[i*(n_genes+1),2]
    else:
        toxicants_name[i] = mtr_dataset[i*(n_genes+1),0]

file_dataset = open("dataset_12dbps.pkl","wb")
pickle.dump([data,genes_name,pathways_type,toxicants_name],file_dataset)
file_dataset.close()
