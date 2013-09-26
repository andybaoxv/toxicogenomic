import csv
import numpy as np
import matplotlib.pyplot as plt

# Load data from the dataset file
file_dataset = open("12 DBPs and other chemicals.csv","rb")
reader = csv.reader(file_dataset)
lines = [line for line in reader]

mtr_dataset = np.array(lines)

## Define the dimensions of the data matrix
# Number of Genes used
n_genes = 39

# Number of time points in collecting Gene expression data
n_times = 25
assert mtr_dataset.shape[0]%(n_genes+1) == 0
assert (mtr_dataset.shape[1]-1)%(n_times+1) == 0

# Number of toxicants used
n_toxicants = mtr_dataset.shape[0]/(n_genes+1)

# Number of dose concentration values used
n_concentrations = (mtr_dataset.shape[1]-1)/(n_times+1)

"""
# Choose dataset
for index_toxicants in range(n_toxicants):
    for index_concentrations in range(n_concentrations):
        tmp_1 = (n_genes+1)*index_toxicants
        tmp_2 = (n_times+2)*index_concentrations
        data_tmp = mtr_dataset[tmp_1+1:tmp_1+1+n_genes,tmp_2+2:tmp_2+2+n_times]
        
        data = np.zeros((n_genes,n_times))
        for i in range(n_genes):
            for j in range(n_times):
                data[i,j] = float(data_tmp[i,j])

        # Minus the DC level
        for i in range(n_genes):
            data[i,:] = data[i,:]-np.mean(data[i,:])
        
        # Draw the plot
        fig = plt.figure()

        for i in range(n_genes):
            plt.plot(range(n_times),data[i,:])

plt.show()
"""
