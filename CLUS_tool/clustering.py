#********************************
#          clustering           *
#********************************

# Standard packages
#from netCDF4 import Dataset, num2date, datetime
import numpy as np
import pandas as pd
import sys
#import matplotlib.pyplot as plt
import os

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
numclus       = int(sys.argv[8])  # number of clusters

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
from clus_manipulate import cluster_orderingFREQ, compute_clusterpatterns
from readsavencfield import read3Dncfield, save_N_2Dfields

# OUTPUT DIRECTORY
OUTPUTdir=dir_OUTPUT+'OUTPUT/'
OUTfig=OUTPUTdir+'OUTfig/'
OUTnc=OUTPUTdir+'OUTnc/'
OUTtxt=OUTPUTdir+'OUTtxt/'
OUTpy=OUTPUTdir+'OUTpy/'


# LOAD DATA PROCESSED BY precompute.py
#______________________________________
anomaly_filenames='area_anomaly_'+name_outputs
fn = [i for i in os.listdir(OUTnc) \
    if os.path.isfile(os.path.join(OUTnc,i)) and anomaly_filenames in i]
fn.sort()
print('\nInput {0} anomaly fields netCDF files (computed by precompute.py):'.format(numens))
print('\n'.join(fn))
var_ensList=[]
for ens in range(numens):
    #print(fn[ens])
    ifile=OUTnc+fn[ens]
    var_area, lat_area, lon_area, dates, time_units, var_units= read3Dncfield(ifile)
    var_ensList.append(var_area)

syr=pd.to_datetime(dates).year[0]
eyr=pd.to_datetime(dates).year[-1]

print('\nINPUT DATA: {0} anomaly fields netCDF files (computed by precompute.py)'.format(numens))
print('Data cover the time period from {0} to {1}'.format(syr,eyr))
print('Data dimensions for one ensemble member are {0}'.format(var_ensList[0].shape))
print('Number of ensemble members is {0}'.format(len(var_ensList)))

name_outputs=name_outputs+'_{0}pcs'.format(numpcs)

if enstoselect!='no':
    print('Ensemble member number {0} is selected for the analysis'.format(enstoselect))
    var_ens=[]
    var_ens.append(var_ensList[enstoselect])
    print(len(var_ens),var_ens[0].shape)
    del var_ensList
    var_ensList=var_ens
    name_outputs=name_outputs+'_'+str(enstoselect)
    print(name_outputs)
    numens=1
else:
    print('All the ensemble members are concatenated one after the other prior to the analysis')


# LOAD DATA PROCESSED BY compute.py
#______________________________________
name_savedvar='{0}varfrac_{1}.txt'.format(OUTtxt,name_outputs)
varfrac=np.loadtxt(name_savedvar)
acc=np.cumsum(varfrac*100)
if perc!='no':
    # Find how many PCs explain a certain percentage of variance
    # (find the mode relative to the percentage closest to perc, but bigger than perc)
    numpcs=min(enumerate(acc), key=lambda x: x[1]<=perc)[0]+1
    print('\nThe number of PCs that explain the percentage closest to {0}% of variance (but grater than {0}%) is {1}'.format(perc,numpcs))
    exctperc=min(enumerate(acc), key=lambda x: x[1]<=perc)[1]
if numpcs!='no':
    exctperc=acc[numpcs-1]
print('(the first {0} PCs explain exactly the {1}% of variance)'.format(numpcs,"%.2f" %exctperc))


# load PCunscal
name_savedvar='{0}PCunscal_{1}.txt'.format(OUTtxt,name_outputs)
PCunscal=np.loadtxt(name_savedvar)
print('(times, numPCs) ={0}'.format(PCunscal.shape))

pc=np.transpose(PCunscal)


# k-means analysis using the subset of PCs
#______________________________________
print('number of clusters: {0}'.format(numclus))

npart=100      #100
varopt=[]
indcl0=[]
centr=[]
print('npart ={0}'.format(npart))
for ncl in range(2,7):
    nfcl,indcl1,centr1,varopt1,iseed=ctool.cluster_toolkit.clus_opt(ncl,npart,pc)
    indcl0.append(np.subtract(indcl1,1))        #indcl1 starts from 1, indcl0 starts from 0
    centr.append(centr1)
    varopt.append(varopt1)
indcl=indcl0[numclus-2]
centr=centr[numclus-2]

# save cluster index
namef='{0}indcl_{1}clus_{2}.txt'.format(OUTtxt,numclus,name_outputs)
np.savetxt(namef,indcl,fmt='%d')

# save cluster centroids
namef='{0}centr_{1}clus_{2}.txt'.format(OUTtxt,numclus,name_outputs)
np.savetxt(namef,centr)

# save cluster optimal variance ratio (this is needed for significance computation: clusters_sig.py)
namef='{0}varopt_2to6clus_{2}.txt'.format(OUTtxt,numclus,name_outputs)
np.savetxt(namef,varopt, fmt='%1.10f')


# Cluster ordering in decreasing frequency
#_______________________
centrORD,indclORD=cluster_orderingFREQ(indcl,centr,numclus)

# save cluster index
namef='{0}indclORD_{1}clus_{2}.txt'.format(OUTtxt,numclus,name_outputs)
np.savetxt(namef,indclORD,fmt='%d')

# save cluster centroids
namef='{0}centrORD_{1}clus_{2}.txt'.format(OUTtxt,numclus,name_outputs)
np.savetxt(namef,centrORD)

# COMPUTE AND SAVE CLUSTERS PATTERNS
#_______________________
# compute cluster patterns
cluspattORD=compute_clusterpatterns(numclus,var_ensList,indclORD)
# save cluster patterns
varsave='cluspattern'
ofile='{0}cluspatternORD_{1}clus_{2}.nc'.format(OUTnc,numclus,name_outputs)
save_N_2Dfields(lat_area,lon_area,np.array(cluspattORD),varsave,var_units,ofile)

print('Cluster pattern netCDF variable (ordered by decreasing frequency of occurrence) is saved as\n{0}'.format(ofile))
print('____________________________________________________________________________________________________________________')


print('\n******************************************************************************')
print('END {0}'.format(sys.argv[0]))
print('*********************************************************************************')


## PLOT AND SAVE CLUSTERS
##_______________________
## plot the cluster patterns, all in one panel
#namef='{0}indclORD_{1}clus_{2}.txt'.format(OUTtxt,numclus,name_outputs)
#indcl=np.loadtxt(namef)
#if enstoselect!='no':
#    tit='{0} {1} {2}{3} {4} {5} {6} ({7}-{8})'.format(varname,model,kind,enstoselect,res,season,area,syr,eyr)
#else:
#    tit='{0} {1} {2} {3} {4} {5} ({6}-{7})'.format(varname,model,kind,res,season,area,syr,eyr)
#ax,cluspatt=computeplot_clusters(area,lon_area,lat_area,numclus,numpcs,var_ensList,indcl,tit)
#
#namef='{0}clus_patterns_{1}clus_{2}.eps'.format(OUTfig,numclus,name_outputs)
#ax.figure.savefig(namef)#bbox_inches='tight')
#print('Clusters eps figure for {0} weather regimes is saved as\n{1}'.format(area,namef))
#print('____________________________________________________________________________________________________________________')
#
## save cluster patterns
#varsave='cluspattern'
#ofile='{0}cluspattern_{1}clus_{2}.nc'.format(OUTnc,numclus,name_outputs)
#save_N_2Dfields(lat_area,lon_area,np.array(cluspatt),varsave,'m',ofile)
#print('Cluster pattern netCDF variable for {0} weather regimes is saved as\n{1}'.format(area,ofile))
#print('____________________________________________________________________________________________________________________')
# quick plot
#tit='cluster pattern 1'
#____________Plot the lon-lat map (Orthographic Projection)
#fig = plot_ortho(cluspatt[0], lat_area, lon_area, clat=50, clon=0, tit=tit)
#fig.show()
##plt.show(block=True)

