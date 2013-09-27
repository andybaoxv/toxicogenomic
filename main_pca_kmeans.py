import numpy as np
import pickle

from nitime import utils
from nitime import algorithms as alg
from nitime.timeseries import TimeSeries
from nitime.viz import plot_tseries

from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

# Read dataset and related information
file_dataset = open("dataset_12dbps.pkl","rb")
data,genes_name,pathways_type = pickle.load(file_dataset)
file_dataset.close()

# Obtain size of the data
n_toxicants,n_concentrations,n_genes,n_times = data.shape

#########################################################################
# Specify the order of autoregressive model
n_ar_order = 3
# Specify the number of clusters of conditions 
n_clusters = 4
# Apply PCA(True) or not(False)
flag_pca = False
# Visualization dimension(2 or 3)
visual_dim = 3
#########################################################################

# Fit the curve of each gene under each condition(1 toxicant+ 1 concentration)
n_conditions = n_toxicants*n_concentrations
n_coefs = n_ar_order*n_genes
mtr_coef = np.zeros((n_conditions,n_coefs))
sigma_est_avg = 0

for i in range(n_toxicants):
    for j in range(n_concentrations):
        for k in range(n_genes):
            coefs_est, sigma_est = alg.AR_est_YW(data[i,j,k,:],n_ar_order)
            mtr_coef[i*n_concentrations+j,n_ar_order*k:n_ar_order*(k+1)] = \
                    coefs_est
            sigma_est_avg += sigma_est**2

# Obtain average noise level
sigma_est_avg = sigma_est_avg/n_conditions
print("AR order: "+str(n_ar_order)+"\tNoise level: "+str(sigma_est_avg))

# n_instances = n_conditions, n_features = n_coefs
n_instances,n_features = mtr_coef.shape

# Apply PCA on the matrix of coeficient dataset
if flag_pca == True:
    clf_pca = PCA(n_components=3,copy=True,whiten=False)
    mtr_coef = clf_pca.fit_transform(mtr_coef)

# Apply K-means clustering on the dataset
clf_kmeans = KMeans(n_clusters=n_clusters,init='random',n_jobs=-1)
labels = clf_kmeans.fit_predict(mtr_coef)

# Visualize the clustering result in 3D
clf_pca_1 = PCA(n_components=visual_dim,copy=True,whiten=False)
mtr_coef_ld = clf_pca_1.fit_transform(mtr_coef)
fig = plt.figure()
color = ['g','y','r','k','y','k','c','w']

if visual_dim == 3:
    ax = fig.add_subplot(111,projection='3d')
    for i in range(n_instances):
        ax.scatter(mtr_coef_ld[i,0],mtr_coef_ld[i,1],mtr_coef_ld[i,2],\
                c=color[labels[i]])
elif visual_dim == 2:
    for i in range(n_instances):
        plt.scatter(mtr_coef_ld[i,0],mtr_coef_ld[i,1],c=color[labels[i]])
else:
    print "ERROR: visual_dim must be 2 or 3"
plt.show()
