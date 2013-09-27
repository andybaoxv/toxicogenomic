import numpy as np
import pickle

from nitime import utils
from nitime import algorithms as alg
from nitime.timeseries import TimeSeries
from nitime.viz import plot_tseries

# Read dataset and related information
file_dataset = open("dataset_12dbps.pkl","rb")
data,genes_name,pathways_type = pickle.load(file_dataset)
file_dataset.close()

# Obtain size of the data
n_toxicants,n_concentrations,n_genes,n_times = data.shape

#########################################################################
# Specify the order of autoregressive model
n_ar_order = 5
#########################################################################

# Fit the curve of each gene under each condition(1 toxicant+ 1 concentration)
n_conditions = n_toxicants*n_concentrations
n_coefs = n_ar_order*n_genes
mtr_coef = np.zeros((n_conditions,n_coefs))

for i in range(n_toxicants):
    for j in range(n_concentrations):
        for k in range(n_genes):
            coefs_est, sigma_est = alg.AR_est_YW(data[i,j,k,:],n_ar_order)
            mtr_coef[i*n_concentrations+j,n_ar_order*k:n_ar_order*(k+1)] = \
                    coefs_est

# Apply PCA on the matrix of coeficient dataset

