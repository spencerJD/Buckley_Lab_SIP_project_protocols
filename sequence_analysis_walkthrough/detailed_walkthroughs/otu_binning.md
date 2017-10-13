# OTU Binning

### Initialize your Notebook and Set Variables
* First we must initialize your notebook, set important variables, and load all necessary libraries.

  ```r
workDir = '/home/USERNAME/NOTEBOOK/PROJECT_NAME/SUBDIRECTORY_WHERE_DATA_IS_LOCATED'
baseDir = '/home/USERNAME/NOTEBOOK/PROJECT_NAME/'
databaseDir = '/home/USERNAME/NOTEBOOK/PROJECTNAME/data/databases/'
nprocs = 24
  ```
  * NOTE: Everyone organizes their data  bit differently, but in this case, we will assume that you have a general directory for this particular project, and a subdirectory within that directory where you are performing this OTU binning. Also, after you have performed the SILVA database formatting in the quality control section, you can add that formatted database to another subdirectory under your main project directory.

  ```r
import os
import re
import glob
import pandas as pd
from qiime.assign_taxonomy import UclustConsensusTaxonAssigner
from IPython.display import Image
from cogent.app.usearch import clusters_from_blast_uc_file
from cogent.parse.fasta import MinimalFastaParser as parse
  ```
  
  ```r
%load_ext rpy2.ipython
  ```
  
  ```r
%%R
library(ggplot2)
library(dplyr)
library(tidyr)
library(gridExtra)
  ```
  ```r
if not os.path.isdir(workDir):
    os.makedirs(workDir)
  ```
  
### Concatenate All Relevant QC Sequence Files
* If your full data set contains samples from multiple sequencing runs, you will need to concatenate your post-QC fasta files together before doing any OTU picking. 
* In the following example, sequence files from multiple runs will be combined together for downstream analysis. The file structure for this dataset contains multiple libraries, some of which have been run multiple times, all of which have passed through quality filtering before this step.
  * If you do not have a similar file structure to this (in terms of directory and file naming schemes), you can manually concatenate all `finalQC.fasta` files from your individual libraries together.

  ```r
p = os.path.join(baseDir, '*V4_Lib*', 'run*', 'QC', 'finalQC.fasta')
fileList = glob.glob(p)
fileList
  ```
  * Example Output:
  
    ```
'/home/nick/notebook/fullCyc/data/MiSeq_16S/515f-806r/150620_V4_Lib1/run1/QC/finalQC.fasta',
'/home/nick/notebook/fullCyc/data/MiSeq_16S/515f-806r/150811_V4_Lib6/run1/QC/finalQC.fasta',
'/home/nick/notebook/fullCyc/data/MiSeq_16S/515f-806r/150715_V4_Lib3/run1/QC/finalQC.fasta',
'/home/nick/notebook/fullCyc/data/MiSeq_16S/515f-806r/150710_V4_Lib2/run1/QC/finalQC.fasta',
'/home/nick/notebook/fullCyc/data/MiSeq_16S/515f-806r/150710_V4_Lib2/run2/QC/finalQC.fasta',
'/home/nick/notebook/fullCyc/data/MiSeq_16S/515f-806r/150727_V4_Lib5/run1/QC/finalQC.fasta',
'/home/nick/notebook/fullCyc/data/MiSeq_16S/515f-806r/150721_V4_Lib4/run1/QC/finalQC.fasta'
    ```

* Next we will print the number of sequences in each file.

  ```r
# number of sequences per file
for f in fileList:
    ret = !cd $workDir; grep -c ">" $f
    print '{} : {}'.format(f, ret[0])
  ```
  * Example Output:

    ```
/home/nick/notebook/fullCyc/data/MiSeq_16S/515f-806r/150620_V4_Lib1/run1/QC/finalQC.fasta : 2660350
/home/nick/notebook/fullCyc/data/MiSeq_16S/515f-806r/150811_V4_Lib6/run1/QC/finalQC.fasta : 3879125
/home/nick/notebook/fullCyc/data/MiSeq_16S/515f-806r/150715_V4_Lib3/run1/QC/finalQC.fasta : 6770040
/home/nick/notebook/fullCyc/data/MiSeq_16S/515f-806r/150710_V4_Lib2/run1/QC/finalQC.fasta : 6463568
/home/nick/notebook/fullCyc/data/MiSeq_16S/515f-806r/150710_V4_Lib2/run2/QC/finalQC.fasta : 3178580
/home/nick/notebook/fullCyc/data/MiSeq_16S/515f-806r/150727_V4_Lib5/run1/QC/finalQC.fasta : 2134224
/home/nick/notebook/fullCyc/data/MiSeq_16S/515f-806r/150721_V4_Lib4/run1/QC/finalQC.fasta : 4237418
    ```

* These files are then concatenated together into the final QC fasta file that will be used in OTU binning.

  ```r
files = ' '.join(fileList)
!cd $workDir; \
    cat $files > finalQC.fasta
  ```
* The total number of sequences can then be printed.

  ```r
ret = !cd $workDir; grep -c ">" finalQC.fasta
print 'Number of sequences: {}'.format(ret[0])
  ```
  * Example Output:
  
    ```
Number of sequences: 29323305
    ```

### Generate a List of Only Unique Sequences
* Like we did in the previous quality control steps, we now want to generate a list of only the unique sequences in our dataset. This will reduce the time and computing power needed for OTU binning.

  ```r
!cd $workDir; \
    mothur "#unique.seqs(fasta=finalQC.fasta)" | head -n 50
  ```
  
* We can then print the total number of unique sequences in our dataset.

  ```r
!printf "Number of unique sequences: "
!cd $workDir; grep -c ">" finalQC.unique.fasta
  ```
  
  * Example Output:

    ```
Number of unique sequences: 16068344
    ```

* Next we must reformat the names of our sequences from the Mothur format in order to work with them in USEARCH.

  ```r
counts = {}
inFile = os.path.join(workDir, 'finalQC.names')
with open(inFile) as iFH:
    for line in iFH:
        seedID, seqIDs = line.split("\t")
        count = len(seqIDs.split(","))
        counts[seedID] = count 
  ```
  
  ```r
outFile = os.path.join(workDir, 'finalQC.unique.usearch_names.fasta')
with open(outFile, 'w') as oFH:
    inFile = os.path.join(workDir, 'finalQC.unique.fasta')
    for n, s in parse(open(inFile)):
        if counts[n] > 1:
            oFH.write(">%s;size=%s;\n%s\n"%(n,counts[n],s))
        else:
            continue
  ```
  
  ```r
!cd $workDir; \
    head -n 6 finalQC.unique.usearch_names.fasta
  ```
  
  * Example Output: 
  
    ```
>13C-Pal.D48.R2_F28_10;size=464;
TACAGAGGTCTCAAGCGTTGTTCGGATTCATTGGGCGTAAAGGGTGCGTAGGCGGCGCGGCAAGTCGGGTGTGAAATCTCGGAGCTTAACTCCGAAACTGCATTCGATACTGCCGTGCTTGAGGACTGGAGAGGAGACTGGAATTTACGGTGTAGCGGTGAAATGCGTAGATATCGTAAGGAAGACCAGTGGCGAAGGCGGGTCTCTGGACAGTTCCTGACGCTGAGGCACGAAGGTCAGGGGAGCAAACGGG
>13C-Van.D30.R1_F24_13;size=9669;
TACGGAGGGTGCAAGCGTTAATCGGAATTACTGGGCGTAAAGCGCACGCAGGCGGTCTGTCAAGTCGGATGTGAAATCCCCGGGCTCAACCTGGGAACTGCATTCGAAACTGGCAGGCTAGAGTCTTGTAGAGGGGGGTAGAATTCCAGGTGTAGCGGTGAAATGCGTAGAGATCTGGAGGAATACCGGTGGCGAAGGCGGCCCCCTGGACAAAGACTGACGCTCAGGTGCGAAAGCGTGGGGAGCAAACAGG
>13C-Van.D30.R1_F24_21;size=2639;
TACGTAGGTGGCAAGCGTTATCCGGAATTATTGGGCGTAAAGCGCGCGCAGGTGGTTTCTTAAGTCTGATGTGAAAGCCCACGGCTCAACCGTGGAGGGGCATTGGAAACTGGGAGACTTGAGTGCAGAAGAGGAAAGTGGAATTCCATGTGTAGCGGTGAAATGCGTAGAGATATGGAGGAACACCAGTGGCGAAGGCGACTTTCTGGTCTGTAACTGACACTGAGGCGCGAAAGCGTGGGGAGCAAACAGG
    ```

### `USEARCH` Pipeline
* In this section, we will perform OTU clustering with out unique sequences file.
* We will also remove any singleton sequences from our dataset.
* First, we will sort our sequences by size. In this case, 'size' is the total number of times that a particular unique sequence appears in our dataset. 

  ```r
!cd $workDir; \
    usearch \
    -sortbysize finalQC.unique.usearch_names.fasta \
    -output finalQC_uniques_sorted.fasta \
    -minsize 
  ```

* We can then print out the first and last few sequences in order to check that this sorting worked.

  ```r
!cd $workDir; \
    head finalQC_uniques_sorted.fasta; \
    tail finalQC_uniques_sorted.fasta
  ```
  
  * Example Output:
  
    ```
>13C-Xyl.D3.R2_F14_343;size=332811;
TACAGAGGGTGCAAGCGTTAATCGGAATTACTGGGCGTAAAGCGCGCGTAGGTGGTTCGTTAAGTTGGATGTGAAAGCCC
CGGGCTCAACCTGGGAACTGCATTCAAAACTGACGAGCTAGAGTATGGTAGAGGGTGGTGGAATTTCCTGTGTAGCGGTG
AAATGCGTAGATATAGGAAGGAACACCAGTGGCGAAGGCGACCACCTGGACTGATACTGACACTGAGGTGCGAAAGCGTG
GGGAGCAAACAGG
>13C-Xyl.D1.R2_F9_927;size=247102;
TACAGAGGGTGCAAGCGTTAATCGGAATTACTGGGCGTAAAGCGCGCGTAGGTGGTTCGTTAAGTTGGATGTGAAATCCC
CGGGCTCAACCTGGGAACTGCATCCAAAACTGGCGAGCTAGAGTAGGGCAGAGGGTGGTGGAATTTCCTGTGTAGCGGTG
AAATGCGTAGATATAGGAAGGAACACCAGTGGCGAAGGCGACCACCTGGGCTCATACTGACACTGAGGTGCGAAAGCGTG
GGGAGCAAACAGG
>13C-Gly.D6.R2_F17_1117646;size=2;
TACAGAGGGTGCAAGCGTTACTCGGAATTACTGGGCGTAAAGCGTGCGTAGGTGGTTGTTTAAGTCTGTTGTGAAAGCCC
TGGGCTCAACCTGGGAACTGCAGTGGAAACTGGACAACTAGAGTGTGGTAGAGGGTAGCGGAATTCCCGGTGTAGCAGTG
AAATGCGTAGAGATCGGGAGGAACATCCGTGGCGAAGGCGACTACCTGGACCAGCACTGACACTGAGGCACGAAAGCGTG
GGGAGCAAACAGG
>12C-Con.D14.R1_F10_1117656;size=2;
TACGTAGGGGGCAAGCGTTGTCCGGAATTATTGGGCGTAAAGAGCGTGTAGGCGGCCAGATAAGTCCGCTGTGAAAACTC
GAGGCTCAACTTCGAGCTGTCGGTGGAAACTATTTGGCTAGAGTCCGGAAGAGGAGAGTGGAATTCCTGGTGTAGCGGTG
AAATGCGCAGATATCAGGAAGAACACCCGTGGCGAAGGCGGCTCTCTGGTACGGAACTGACGCTGAGACGCGAAAGCGTG
GGGAGCGAACAGG
    ```

* Now we are able to perform OTU clustering using our sorted unique sequences file.

  ```r
!cd $workDir; \
    usearch \
    -cluster_otus finalQC_uniques_sorted.fasta \
    -otus otus.fasta
  ```

* At this point we want to rename these OTUs with simple OTU numbers for downstream applications.

  ```r
%%bash -s "$workDir"
cd $1

bioawk -c fastx '{print ">" "OTU" "." NR "\n" $seq}' otus.fasta > otusn.fasta
  ```
  
  ```r
!cd $workDir; head -n 4 otusn.fasta
  ```
  
  * Example Output:
   
    ```
>OTU.1
TACAGAGGGTGCAAGCGTTAATCGGAATTACTGGGCGTAAAGCGCGCGTAGGTGGTTCGTTAAGTTGGATGTGAAAGCCCCGGGCTCAACCTGGGAACTGCATTCAAAACTGACGAGCTAGAGTATGGTAGAGGGTGGTGGAATTTCCTGTGTAGCGGTGAAATGCGTAGATATAGGAAGGAACACCAGTGGCGAAGGCGACCACCTGGACTGATACTGACACTGAGGTGCGAAAGCGTGGGGAGCAAACAGG
>OTU.2
TACGTAGGTGGCAAGCGTTGTCCGGAATTATTGGGCGTAAAGCGCGCGCAGGTGGTTCCTTAAGTCTGATGTGAAAGCCCACGGCTCAACCGTGGAGGGTCATTGGAAACTGGGGAACTTGAGTGCAGAAGAGGAAAGTGGAATTCCAAGTGTAGCGGTGAAATGCGTAGAGATTTGGAGGAACACCAGTGGCGAAGGCGACTTTCTGGTCTGTAACTGACACTGAGGCGCGAAAGCGTGGGGAGCAAACAGG
    ```
  
* Next, we want to remove any Chloroplast, Eukaryal, and Archaeal sequences from our OTU list.
* We start by using the SILVA database to assign taxonomy to out OTUs.

  ```r
%%bash -s "$workDir" "$databaseDir" "$nprocs"

cd $1

parallel_assign_taxonomy_uclust.py \
    -r $2'97_Silva_111_rep_set_no_ambig.fasta' \
    -t $2'Silva_111_taxa_map_full.txt' \
    -O $3 \
    -i otusn.fasta \
    -o otusn_tax
  ```

* Then, we make a list of all those OTUs which fall into Chloroplast, Eukaryota, Archaea, unassigned, or Mitochondria classifications.

  ```r
%%bash -s "$workDir" "$databaseDir"

cd $1

egrep "Chloroplast|Eukaryota|Archaea|Unassigned|mitochondria" \
    otusn_tax/otusn_tax_assignments.txt |\
    awk '{print $1}' > to_remove_tax.accnos
  ```

* Next, we remove those sequences from our OTU file.

  ```r
!cd $workDir; \
    mothur "#remove.seqs(fasta=otusn.fasta, accnos=to_remove_tax.accnos)" | head -n 50
  ```

* We can then check to see how many sequences were filtered out this way.

  ```r
!printf "Pre-filter: number of sequences: "
!cd $workDir; grep -c ">" otusn.fasta
!printf "Post-filter: number of sequences: "
!cd $workDir; grep -c ">" otusn.pick.fasta
  ```
  
  * Example Output:
  
    ```
Pre-filter: number of sequences: 13790
Post-filter: number of sequences: 12163
    ```

### Mapping Reads
* Now that we have our OTUs identified, we can determine how many times these OTUs appear in each sample.
* First, we have to rename our sequences in our total sequence file into `USEARCH` formatting. 

  ```r
%%bash -s "$workDir"
cd $1

perl -pe 's/^>(.+)(_[^_]+)\n$/>$1$2\_$.;barcodelabel=$1\n/' finalQC.fasta > finalQC_usearchfmt.fasta
  ```

* We can check to make sure that these are formatted correctly as well.

  ```r
!cd $workDir; head -n 6 finalQC_usearchfmt.fasta
  ```

  * Example Output:
  
    ```
>13C-Pal.D14.R1_F8_0_1;barcodelabel=13C-Pal.D14.R1_F8
TACGGAGGGTGCAAGCGTTGTCCGGATTTATTGGGTTTACAGGGAGCGTAGGCGGTCTTTTAAGTCAGTGGTGAAATCCTCCAGCTTAACTGGAGAACTGCCATTGATACTGAAAGACTTGAGTACAGACGAGGTAGGCGGAATTGACAGTGTAGCGGTGAAATGCATAGATATTGGCAAGAACACCGATTGCGAAGGCAGCTTACTAGACTGTAACTGACGCTGAGGCTCGAAAGTGCGGGGATCAAACAGG
>13C-Pal.D14.R1_F7_1_3;barcodelabel=13C-Pal.D14.R1_F7
TACGGAGGGGGCTAGCGTTGTTCGGAATTACTGGGCGTAAAGCGCACGTAGGCCTTTGTAAGTTAGAGGTGAAAGCCCGGAGCTCAACTCCGGAACTGCCTTTAAGACTGCATCGCTTGAACGTCGGAGAGGTAAGTGGAATTCCGCGTGTAGAGGTGAAATTCGTAGATATTCGGAAGAACACCAGTGGCGAAGGCGACTTACTGGACGACTGTTGACGCTGAGGTGCGAAAGCGTGGGGAGCAAACAGG
>13C-Pal.D48.R2_F25_2_5;barcodelabel=13C-Pal.D48.R2_F25
TACCAGCACCCCGAGTGGTCGGGACGGTTATTGGGCCTAAAGCATCCGTAGCCGGTTCTACAAGTCTTCCGTTAAATCCACCTGCTTAACAGATGGGCTGCGGAAGATACTATAGAGCTAGGAGGCGGGAGAGGCAAGCGGGACTCGATGGGTAGGGGTAACATCCGTTGATCCATTGAAGACCGCCAGTGGCGAAGGCGGCTTGCCAGAACGCGCTCGACGGTGAGGGATGAAAGCTGGGGGAGCAAACCGG
    ```

* Our full sequence file is quite large, and must be split into multiple files in order to run efficiently.
  * We can check the file size using the following code:
  
    ```r
!cd $workDir; \
    du -h finalQC_usearchfmt.fasta
    ```
    * Example Output:
    
      ```
8.8G	finalQC_usearchfmt.fasta
      ```

  * Splitting the file for downstream applications.
  
    ```r
#spliting file
!cd $workDir; \
    pyfasta split -n 5 finalQC_usearchfmt.fasta
    ```

* Now we can create a file list of these split files.

  ```r
g = os.path.join(workDir, 'finalQC_usearchfmt.*.fasta')
fileList = glob.glob(g)
fileList  
  ```
  
  * Example Output:
  
    ```
'/home/nick/notebook/fullCyc/data/MiSeq_16S/515f-806r/V4_Lib1-6/OTU_binning/finalQC_usearchfmt.0.fasta',
'/home/nick/notebook/fullCyc/data/MiSeq_16S/515f-806r/V4_Lib1-6/OTU_binning/finalQC_usearchfmt.4.fasta',
'/home/nick/notebook/fullCyc/data/MiSeq_16S/515f-806r/V4_Lib1-6/OTU_binning/finalQC_usearchfmt.3.fasta',
'/home/nick/notebook/fullCyc/data/MiSeq_16S/515f-806r/V4_Lib1-6/OTU_binning/finalQC_usearchfmt.1.fasta',
'/home/nick/notebook/fullCyc/data/MiSeq_16S/515f-806r/V4_Lib1-6/OTU_binning/finalQC_usearchfmt.2.fasta'
    ```

* We can now run `USEARCH` on these split files in parallel. This step uses our OTU file generated earlier as a database for clustering.

  ```r
# running usearch on each split file
for f in fileList:
    sys.stderr.write('Processing {}\n'.format(f))
    
    ff,_ = os.path.splitext(f)
    _,i = os.path.splitext(ff)
    uc = 'readmap{}.uc'.format(i.lstrip('.')) 
    
    !cd $workDir; \
        usearch \
        -usearch_global $f \
        -db otusn.pick.fasta \
        -strand plus -id 0.97 \
        -uc $uc \
        -threads $nprocs
  ```

* After creating these readmap files for each split sequence file, we can concatenate them together.

  ```r
!cd $workDir; \
    cat readmap[0-9].uc > readmap_all.uc
  ```

* At this point, we must convert our readmap file from `USEARCH` format into biom format for downstream application in the phyloseq R package.
  * Converting from `USEARCH` format to OTU table format:
  
    ```r
!cd $workDir; \
    python /opt/bioinfo/edgar_python_scripts/uc2otutab.py readmap_all.uc > otu_table.txt
    ```
  
  * Converting from an OTU table into a biom table:
  
    ```r
%%bash -s "$workDir"

cd $1

if [ -f otu_table.biom ]; then
    rm otu_table.biom
fi 

biom convert -i otu_table.txt -o otu_table.biom --table-type "otu table"
    ```

  * Create and print a biom table summary:
  
    ```r
%%bash -s "$workDir"

cd $1

if [ -f otu_table_summary.txt ]; then
    rm otu_table_summary.txt
fi 

biom summarize-table -i otu_table.biom -o otu_table_summary.txt
    ```


     ```r
!cd $workDir; cat otu_table_summary.txt
     ```

* Finally, we will add taxonomy metadata to our biom table.

  ```r
!cd $workDir;\
    biom add-metadata -i otu_table.biom \
        -o otu_table_wtax.biom \
        --observation-metadata-fp otusn_tax/otusn_tax_assignments.txt \
        --sc-separated taxonomy \
        --float-fields consensus \
        --int-fields numhits \
        --observation-header OTUID,taxonomy,consensus,numhits
  ```

* Now we are ready for all downstream analyses!





