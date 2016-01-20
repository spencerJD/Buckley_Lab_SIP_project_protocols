# Quality Control of Sequences

In this walkthrough, we will go through the steps necessary to clean up your sequences and discard any that do not meet our quality control standards.

### Seting variables and Initializing Jupyter Notebook
* Start your notebook by setting your directory and file variables, and the threshold value for maximum expected error.

  ```python
workDir = '/var/seq_data/PROJECT_NAME/QC/'
seqFile = 'pear_merged-YYYY-MM-DD.assembled.dmult.fastq'

# number of processors
nprocs = 24

# max expected error
maxee = 1
  ```
  
* Load all libraries needed for this notebook.

  ```python
from screed.fasta import fasta_iter
from pandas import DataFrame
import os
import re
import pandas as pd
from cogent import DNA
from qiime.assign_taxonomy import UclustConsensusTaxonAssigner
from IPython.display import Image
  ```

  ```python
%load_ext rpy2.ipython
  ```

  ```python
%%R
library(ggplot2)
library(dplyr)
library(tidyr)
library(gridExtra)
  ```

  ```python
%pylab inline
  ```

  ```python
if not os.path.isdir(workDir):
    os.mkdir(workDir)
  ```

  ```r
%%bash -s "$workDir" "$seqFile"
  ```

  ```r
ln -s -f $1../$2 $1
  ```

### Discard Sequences That Exceed the Maximum Expected Error Rate
* We will cycle through all sequences with `USEARCH` to filter out those sequences with an error rate above our maximum expected rate.
* This process is performed in parallel, which the majority of the following code is structuring.

  ```r
%%bash -s "$workDir" "$seqFile" "$nprocs" "$maxee"

tmpdir1=`mktemp -d`
trap "rm -r $tmpdir1" 1 2 3 15
split -d -l 2000000 $1$2 $tmpdir1/Block

tmpdir2=`mktemp -d`
trap "rm -r $tmpdir2" 1 2 3 15
ls $tmpdir1/Block?? | parallel --gnu -j $3 -k "usearch -fastq_filter {} -fastq_maxee $4 \
-fastaout $tmpdir2/{#}.fasta >/dev/null 2>&1 && cat $tmpdir2/{#}.fasta" > $1$2_maxee1.fasta
rm -r $tmpdir2 $tmpdir1
  ```







