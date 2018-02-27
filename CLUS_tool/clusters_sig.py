#*********************************
#         clusters_sig           *
#*********************************
'''
H_0: There are no regimes ---> multi-normal distribution PDF
Synthetic datasets modelled on the PCs of the original data are computed
(synthetic PCs have the same lag-1, mean and standard deviation of the original PCs)
SIGNIFICANCE = % of times that the optimal variance ratio found by clustering the real dataset
               exceed the optimal ratio found by clustering the syntetic dataset
'''

# Standard packages
#from netCDF4 import Dataset, num2date, datetime
import numpy as np
#import pandas as pd
import sys
import pickle
#import matplotlib.pyplot as plt
import os
import datetime


print('**************************************************************')
print('Running {0}'.format(sys.argv[0]))
print('**************************************************************')
dir_WRtool    = sys.argv[1]  # WRTOOL directory
dir_OUTPUT    = sys.argv[2]  # OUTPUT directory
name_outputs  = sys.argv[3]  # name of the outputs
numens        = int(sys.argv[4])  # total number of ensemble members
numpcs        = sys.argv[5]  # number of principal components
perc          = sys.argv[6]  # percentage of explained variance for EOF analysis
enstoselect   = sys.argv[7]  # ensemble member to analyse
season        = sys.argv[8]  # selected season
freq          = sys.argv[9]  # frequency of data

if numpcs!='no':
    numpcs=int(sys.argv[5])
if perc!='no':
    perc=int(sys.argv[6])
if (perc=='no' and numpcs=='no') or (perc!='no' and numpcs!='no'):
    raise ValueError('You have to specify either "perc" or "numpcs".')
if enstoselect!='no':
    enstoselect=int(sys.argv[7])

# User-defined libraries
sys.path.insert(0,dir_WRtool+'WRtool/')
import ctool
import ctp

# OUTPUT DIRECTORY
OUTPUTdir=dir_OUTPUT+'OUTPUT/'
OUTfig=OUTPUTdir+'OUTfig/'
OUTnc=OUTPUTdir+'OUTnc/'
OUTtxt=OUTPUTdir+'OUTtxt/'
OUTpy=OUTPUTdir+'OUTpy/'

name_outputs=name_outputs+'_{0}pcs'.format(numpcs)

if enstoselect!='no':
    print('Ensemble member number {0} is selected for the analysis'.format(enstoselect))
    name_outputs=name_outputs+'_'+str(enstoselect)
    print(name_outputs)
    numens=1
else:
    print('All the ensemble members are concatenated one after the other prior to the analysis')


# LOAD DATA PROCESSED BY clusters_comparison.py
#______________________________________

# load PCunscal
name_savedvar='{0}PCunscal_{1}.txt'.format(OUTtxt,name_outputs)
PCunscal=np.loadtxt(name_savedvar)
print('(times, numPCs) ={0}'.format(PCunscal.shape))

# load cluster optimal variance ratio
namef='{0}varopt_2to6clus_{1}.txt'.format(OUTtxt,name_outputs)
varopt=np.loadtxt(namef)
print('The optimal ratio for cluster partition is the ratio of the variance among cluster centroids and the intra-cluster variance')
print('Optimal ratio for cluster partition of {0} clusters = {1}'.format(PCunscal.shape[1],varopt))

npart=100      #100
nrsamp=500     #500
print('\nNumber of synthetic PCs used for significance computation = {0}'.format(nrsamp))
if freq=='day':
    ndis=int(PCunscal.shape[0]/(30*len(season)))-1        #28 (e.g. DJF data: ndis=number of winters-1)

elif freq=='mon':
    ndis=int(PCunscal.shape[0]/(12*len(season)))-1
else:
    print('input data frequency is neither day nor mon')
print('check: number of years={0}'.format(ndis+1))


# Compute significance

pc=np.transpose(PCunscal)

#=======serial=========
#significance=ctool.cluster_toolkit.clus_sig(nrsamp,npart,ndis,pc,varopt)

#=======parallel=======
start = datetime.datetime.now()
significance=ctp.cluster_toolkit_parallel.clus_sig_p(nrsamp,npart,ndis,pc,varopt)
end = datetime.datetime.now()
print('significance computation took me %s seconds' %(end-start))

print 'significance =', significance
print '{0:.16f}'.format(significance[2])
print 'significance for 4 clusters =', significance[2]

namef='{0}sig_2to6clus_{1}nrsamp_{2}.txt'.format(OUTtxt,nrsamp,name_outputs)
np.savetxt(namef,significance)


print('\n******************************************************************************')
print('END {0}'.format(sys.argv[0]))
print('*********************************************************************************')

