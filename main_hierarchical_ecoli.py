import numpy as np
import pickle
from optparse import OptionParser

from nitime import utils
from nitime import algorithms as alg
from nitime.timeseries import TimeSeries
from nitime.viz import plot_tseries

import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import dendrogram
from my_hierarchical_clustering import my_hierarchical_clustering
from sklearn.metrics.pairwise import euclidean_distances
import copy

desc = """This script permits you to: 
1) Extract features from the time series curve of each gene under each 
condition(1 toxicant + 1 concentration). 
Note that the user need to specify the order of auto-regressive model. 

2) After feature extraction, Euclidean metric could be applied to generate
a distance matrix, then hierarchical clustering is applied.

3) The output is a file containing the dendrogram in png format.
"""

# User Interface with this program
parser = OptionParser(description=desc)
parser.add_option('--order',
        help='Order of auto-regressive model.',
        dest='n_ar_order',nargs=1,
        metavar='<n_ar_order>',type='int',
        default=3)
parser.add_option('--linkage',
        help='Linkage methods in hierarchical clustering:\
                (single,complete,average)',
        dest='method_linkage',nargs=1,
        metavar='<method_linkage>',default='complete')

(options, args) = parser.parse_args()


# Read dataset and related information
file_dataset = open("EcoliData_3Dataset_WithoutMovAvg.pkl","rb")
data,conditions_name = pickle.load(file_dataset)
file_dataset.close()

# Obtain size of the data
n_conditions,n_genes,n_times = data.shape

# Specify the order of autoregressive model
n_ar_order = options.n_ar_order

# Fit the curve of each gene under each condition
n_coefs = n_ar_order*n_genes
mtr_coef = np.zeros((n_conditions,n_coefs))
sigma_est_avg = 0

for i in range(n_conditions):
    for j in range(n_genes):
        coefs_est,sigma_est = alg.AR_est_YW(data[i,j,:],n_ar_order)
        mtr_coef[i,n_ar_order*j:n_ar_order*(j+1)] = coefs_est
        sigma_est_avg += sigma_est**2

# Obtain average noise level
sigma_est_avg = sigma_est_avg/n_conditions/n_genes
print("AR order: "+str(n_ar_order)+"\tAverage Noise Level: "+str(sigma_est_avg))

n_instances,n_features = mtr_coef.shape

# Obtain distance matrix
mtr_dis = euclidean_distances(mtr_coef,mtr_coef)
mtr_sim = np.max(mtr_dis)-mtr_dis
mtr_lin = my_hierarchical_clustering(mtr_sim,method=options.method_linkage)

# Draw dendrogram of hierarchical clustering
fig = plt.figure(figsize=(30,30))
labels = conditions_name
tmp = copy.deepcopy(mtr_lin)
for i in range(tmp.shape[0]):
    tmp[i,2] = np.max(mtr_dis)-tmp[i,2]
dend = dendrogram(tmp,labels=labels,orientation='right')
plt.savefig('dendrogram_hierarchical_ecoli.png')
print("The result is stored in dendrogram_hierarchical_ecoli.png")

