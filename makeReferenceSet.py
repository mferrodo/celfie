## Import libraries
import glob
import os
import re
import sys, getopt
import numpy as np
import pandas as pd
import seaborn as sns
import subprocess

# Parallelisation options
import multiprocessing
from multiprocessing import Process, Manager, Pool
cpuCount = (multiprocessing.cpu_count() - 2)

### BEDtools > 2.27.1 is needed!!
### Tested with py >3.4

os.chdir("/Users/rmvpaeme/Repos/2003_CelFiE/NBL_reference_set/")

# Name of output
trainFileName = "reference_set.bed"

# Folder to store intermediate files
tmp_folder = "/tmp/"

## train files for CancerLocator
# Infinium data
NBL_bismark_folder = "./NBL"
NBL_bismark_files = glob.glob(os.path.join(NBL_bismark_folder, "*.cov.gz"))

cfDNA_bismark_folder = "./cfDNA"
cfDNA_bismark_files = glob.glob(os.path.join(cfDNA_bismark_folder, "*.cov.gz"))

# The file containing the features (= the intersect between HM450K data and RRBS data, see GitHub README)
clusters = pd.read_csv("/Users/rmvpaeme/Repos/2003_CelFiE/NBL_reference_set/output/results_cfDNA_NBL_regs_nochr.csv", sep="\t",usecols=[0,1,2], skiprows=[0], header=None, index_col=None)
clusters[3] = clusters.index
clusterFile = "RRBSregions"
clusters.to_csv(tmp_folder + "%s.txt" % clusterFile, header=None, index=None, sep='\t', mode = 'w')
clusters = clusters.drop([3], axis = 1) # Use empty index to later extract all the clusters from, so that every sample has the same number of clusters


print("Generating %s" % trainFileName)


def generateTrain_NGS(inputfile, label, file_name):

    #inputfile = './NBL/DNA044134_S32.cov.gz'
    #file_name = 'a'
    df = pd.read_csv(inputfile, sep="\t",usecols=[0,1,2,3,4,5], header=None, compression = "gzip")
    df[3] = df[3]/100
    df.to_csv(tmp_folder + "%s.txt" % file_name , header=None, index=None, sep='\t', mode = 'w')

    outfile = open(tmp_folder + '%s_intersect.txt' % file_name, 'w')
    print("     Running bedtools intersect on %s.txt..." % file_name)
    arg = "bedtools intersect -wb -b %s%s.txt -a %s%s.txt" % (tmp_folder, clusterFile, tmp_folder, file_name)
    arg = arg.split()
    proc = subprocess.Popen(args=arg, stdout=outfile, universal_newlines=True).wait()


    if type == "individ_cpg":
        df = pd.read_csv(tmp_folder + '%s_intersect.txt' % file_name, sep="\t", usecols = [0,1,2,3,4,5,9], header=None )
        df = df[[0,1,2,3,4,5,9]]
        df.to_csv(tmp_folder + "%s_reordered.txt" % file_name, header=None, index=None, sep='\t', mode = 'w')

        df = pd.read_csv(tmp_folder + "%s_reordered.txt" % file_name, sep="\t", header=None)

        df.sort_values(by=[0,1,2], inplace=True)

        df[7] = df[4] + df[5]         # Get total depth

        df["index"] = df[0].astype(str)+'_'+df[1].astype(str)+'_'+df[2].astype(str)
        df.index = df["index"]
        df.index.name = None
        df = df.drop([0,1,2,3,5,6, "index"], axis = 1)

        os.remove(tmp_folder + '%s_intersect.txt' % file_name)
        os.remove(tmp_folder + "%s.txt" % file_name)

        return df
    else:
        df = pd.read_csv(tmp_folder + '%s_intersect.txt' % file_name, sep="\t", usecols=[6,7,8,3,4,5,9], header=None )
        df = df[[6,7,8,3,4,5,9]]
        df.to_csv(tmp_folder + "%s_reordered.txt" % file_name, header=None, index=None, sep='\t', mode = 'w')

        arg = "bedtools groupby -i %s%s_reordered.txt -g 1-3,7 -c 5,6 -o sum" % (tmp_folder, file_name)
        arg = arg.split()
        outfile = open(tmp_folder + '%s_clustered.txt' % file_name, 'w')
        print("     Running bedtools groupby on %s.txt..." %file_name)
        proc = subprocess.Popen(args=arg, stdout=outfile, universal_newlines=True).wait()

        df = pd.read_csv(tmp_folder + '%s_clustered.txt' % file_name, sep="\t", header=None, index_col = 3 )
        df.index.name = None

        df[6] = df[4]/(df[4] + df[5])  # Get beta value per cluster

        df = df[[0,1,2,6,4,5]]         # Reorder
        df.sort_values(by=[0,1,2], inplace=True)

        df[7] = df[4] + df[5]         # Get total depth

        df[[7,6,4,5]] = df[[7,6,4,5]].mask(df[7] < 30)         # Mark all clusters lower than 30 reads with NA

        df = df.drop([0,1,2,5,6], axis = 1)

        os.remove(tmp_folder + '%s_intersect.txt' % file_name)
        os.remove(tmp_folder + "%s_clustered.txt" % file_name)
        os.remove(tmp_folder + "%s.txt" % file_name)

        return df

# Similar to import_test_files
with Manager() as manager:
    #Define empty list
    trainFile_list = manager.list()
    tumorGroup_list =  manager.list()

    #trainFile_list = []
#    tumorGroup_list = []

    def import_NBL_train(x):
        file = x
    #    file = './NBL/DNA044133_S31.cov.gz'
    #    file = './NBL/DNA044135_S33.cov.gz'

        file_name = os.path.splitext(os.path.basename(file))[0]
        df = generateTrain_NGS(inputfile = file, label = "NBL", file_name = file_name)
        tumorGroup_list.append(df)

    pool = Pool(cpuCount)
    pool.map(import_NBL_train, NBL_bismark_files)

    tumorGroup = pd.concat(tumorGroup_list, axis = 1)
    tumorGroup[8] = tumorGroup[4].mean(axis= 1)
    tumorGroup[9] = tumorGroup[7].mean(axis= 1)
    tumorGroup = tumorGroup.drop([4,7], axis = 1)
    trainFile_list.append(tumorGroup)

    #tumorGroup_list = []
    tumorGroup_list =  manager.list()
    def import_cfDNA_train(x):
        file = x
        #file = './NBL/DNA044134_S32.cov.gz'
        #file = './NBL/DNA044133_S31.cov.gz'
        file_name = os.path.splitext(os.path.basename(file))[0]
        df = generateTrain_NGS(inputfile = file, label = "cfDNA", file_name = file_name)
        tumorGroup_list.append(df)

    pool = Pool(cpuCount)
    pool.map(import_cfDNA_train, cfDNA_bismark_files)

    tumorGroup = pd.concat(tumorGroup_list, axis = 1)
    tumorGroup[8] = tumorGroup[4].mean(axis= 1)
    tumorGroup[9] = tumorGroup[7].mean(axis= 1)
    tumorGroup = tumorGroup.drop([4,7], axis = 1)
    trainFile_list.append(tumorGroup)


    # Generate full matrix from list
    trainFile = pd.concat(trainFile_list, axis = 1)
    if type == "individ_cpg":
        a = trainFile
        a["index"] = a.index
        a[['chr','start', 'stop']] = a["index"].str.split('_',expand=True)
        a = a.drop(["index"], axis = 1 )
        def move_column_inplace(df, col, pos):
            col = df.pop(col)
            df.insert(pos, col.name, col)
        move_column_inplace(a, "chr", 0)
        move_column_inplace(a, "start", 1)
        move_column_inplace(a, "stop", 2)
        a.sort_values(by=["chr","start","stop"], inplace=True)
        a = a.fillna(np.NaN)
        a.to_csv("./output/%s" % trainFileName, header=None, sep='\t', mode = 'w', index = False, na_rep = np.NaN)
    else:
        # Make sure that the file contains all the clusters
        trainFile = pd.merge(clusters, trainFile, how = "left", left_index=True, right_index=True)
        trainFile_rmNA = trainFile.dropna(axis = 0)
        trainFile_rmNA.to_csv("./output/%s" % trainFileName, header=None, sep='\t', mode = 'w', index = False, na_rep = np.NaN)

    #print("The number of columns in the %sfile is: %i" %  (trainFileName,trainFile.shape[2]))
    #print("The number of columns in the %sfile after removing all NA values is: %i" %  (trainFileName,trainFile_rmNA.shape[2]))
