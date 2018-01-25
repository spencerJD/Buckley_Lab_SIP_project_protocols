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
* We can then print some summary information, such as the total number of sequences before and after filtering, and the first few filtered sequences.

  ```r
%%bash -s "$workDir" "$seqFile" "$nprocs" "$maxee"
​
cd $1
​
printf "Max expected error cutoff: "
echo $4
​
printf "Number of sequence pre-filter: "
grep -c "orig_name=" $2
​
printf "Number of sequences post-filter: "
grep -c ">" $2_maxee1.fasta
​
echo '' 
head -n 8 $2_maxee1.fasta
  ```

* This will result in something like the following:

  ![Example Filter Output](https://cloud.githubusercontent.com/assets/7449496/12463675/046bb3de-bf93-11e5-922f-800e167901df.png)

### Remove Sequences with 'N' Characters
* After filtering out sequences with high error rates, we want to remove all sequences with ambiguous 'N' characters.
* A summary of the total number of sequences before and after this quality control step will also be printed.

  ```r
%%bash -s "$workDir" "$seqFile"

cd $1

bioawk -c fastx '{if ($seq !~ /N/){print ">" $name " " $4 "\n" $seq}}' \
    $2_maxee1.fasta > $2_maxee1_noN.fasta

printf "Number of sequence pre-filter: "
grep -c ">" $2_maxee1.fasta

printf "Number of sequences post-filter: "
grep -c ">" $2_maxee1_noN.fasta
  ```
* Example Output:
  ```
Number of sequence pre-filter: 1667925

Number of sequences post-filter: 1667925
  ```

### Alignment-Based Quality Control Using Mothur
* In this step, we will align our sequences with the SILVA database using Mothur
* This allows us to be certain that our sequences are indeed from the region we expect, and allows for trimming down to our specified sequence length.
* The end result of this step will be our final quality filtered sequence file.

#####Create Unique Sequence File
* First, we will create a file of unique sequences from our dataset. This reduces the size of our sequence file for alignment by only including one copy of each unique sequence. 
* This process returns one file with only unique sequences, and another contining the names of all of the sequences in the total sequence files which are identical to each of the unique sequences.

  ```r
%%bash -s "$workDir" "$seqFile"

cd $1

mothur "#unique.seqs(fasta=$2_maxee1_noN.fasta)" 
  ```

* Then we will check the number of unique sequences in our dataset. This should be significantly less than the total number of sequences we checked previously.

  ```r
%%bash -s "$workDir" "$seqFile"

cd $1

printf "Number of unique sequences: "
grep -c ">" $2_maxee1_noN.unique.fasta
  ```
  * Example Output:
  
    ```
Number of unique sequences: 639570
    ```
  
#####Create Mothur Group File
* Next, we will create the Mothur Group File. This file is used to assign sequences to a specific group. This file consists of two columns, with the sequence name in the first column and the group name in the second column. In our case, the group file will be in the following format: SAMPLENAME_Seq#  SAMPLENAME

  ```r
%%bash -s "$workDir" "$seqFile"
# making mothur group file

cd $1

perl -ne 'if(/^>/){ s/>(.+)(_\d+) .+/$1$2\t$1/; print;}' $2_maxee1_noN.fasta > group_file.txt

head group_file.txt
  ```
  * Example Output:
  
    ```
13C-Oxa.D6.R3_F6_0	13C-Oxa.D6.R3_F6
    ```

#####Download and Format SILVA Database
* In order to make our sequence alignment, we must first download and format the SILVA sequence database to meet our needs.
* First, we will download the database.

  ```r
%%bash -s "$workDir"
# dowloading mothur silva databases

cd $1

if ! [ -d mothur_silva_db ]; then
    mkdir mothur_silva_db
fi

cd mothur_silva_db

if ! [ -e $1silva_ref_aln_mothur.fasta ]; then
    curl -o silva_Bac.zip http://www.mothur.org/w/images/9/98/Silva.bacteria.zip && unzip silva_Bac.zip
    curl -o silva_Euk.zip http://www.mothur.org/w/images/1/1a/Silva.eukarya.zip && unzip silva_Euk.zip
    curl -o silva_Arc.zip http://www.mothur.org/w/images/3/3c/Silva.archaea.zip && unzip silva_Arc.zip
fi
  ```
  ```r
%%bash -s "$workDir"

cd $1'mothur_silva_db'

cat silva.bacteria/silva.bacteria.fasta \
    silva.eukarya/silva.eukarya.fasta \
    Silva.archaea/silva.archaea.fasta \
    > silva_ref_aln_mothur.fasta

printf "Number of sequences in mothur silva fasta: "
grep -c ">" silva_ref_aln_mothur.fasta
  ```
  * Example Output:
    ```
Number of sequences in mothur silva fasta: 18491
    ```

* Then we will run a filter on this database before we make our alignment.

  ```r
%%bash -s "$workDir"
# filtering column positions in silva db

cd $1'mothur_silva_db'

mothur "#filter.seqs(vertical=t, \
    fasta=silva_ref_aln_mothur.fasta, \
    processors=24)" > /dev/null

printf "Number of sequences post-filter: "
grep -c ">" silva_ref_aln_mothur.filter.fasta
   ```
  * Example Output: 
    ```
Number of sequences post-filter: 18491
    ```

#####Align Sequences to the SILVA Database
* This next step will align our sequences with the representative sequences in the SILVA Database.

  ```r
%%bash -s "$workDir" "$seqFile"
# aligning sequences

cd $1

mothur "#align.seqs(candidate=pear_merged\-YYYY\-MM\-DD.assembled.dmult.fastq_maxee1_noN.unique.fasta, \
    template=mothur_silva_db/silva_ref_aln_mothur.filter.fasta, \
    processors=24, \
    flip=T)" | head -n 100
  ```

* Next, we will filter out any gaps in the alignment.

  ```r
%%bash -s "$workDir" "$seqFile"
# filtering out gap positions in the alignment

cd $1

mothur "#filter.seqs(vertical=t, \
    fasta=pear_merged\-YYYY\-MM\-DD.assembled.dmult.fastq_maxee1_noN.unique.align, \
    processors=24)" | head -n 50
  ```
  
* Then we will print a summary of the sequence alignment.

  ```r
%%bash -s "$workDir" "$seqFile"

cd $1

mothur "#summary.seqs(fasta=pear_merged-YYYY-MM-DD.assembled.dmult.fastq_maxee1_noN.unique.filter.fasta, \
    processors=24, \
    name=pear_merged-YYYY-MM-DD.assembled.dmult.fastq_maxee1_noN.names)" 
  ```
  * Example Output (_Note: there is a large amount of header text from Mothur, which I have omitted here_):
    ![Example Mothur Summary](https://cloud.githubusercontent.com/assets/7449496/12484414/36b5f54e-c028-11e5-91de-956d82221790.png)
  * The information from this summary is necessary for the sequence trimmed step after this.

#####Trimming Homopolymers and Sequences That Don't Align to our Amplicon Region
* In this step, we will remove any sequences containing homopolymers with a length greater than 8, and any sequences that do not lie in our amplicon region.
* In the previous step, the median 'start' and 'end' positions will be used to trim off any sequences that do not lie in our amplicon region. Every time you run this step, you must adjust the values in the following code to match these for your specific dataset.

  ```r
%%bash -s "$workDir" "$seqFile"

cd $1

mothur "#screen.seqs(fasta=pear_merged-YYYY-MM-DD.assembled.dmult.fastq_maxee1_noN.unique.filter.fasta, \
    processors=24, \
    name=pear_merged-YYYY-MM-DD.assembled.dmult.fastq_maxee1_noN.names, \
    group=group_file.txt, \
    start=104, \
    end=868, \
    maxhomop=8)" | tail -n 50
  ```

* We will then check the number of sequences remaining after this filtering step.

  ```r
%%bash -s "$workDir" "$seqFile"

cd $1

printf "Number of sequences post-filter: "
grep -c ">" pear_merged-YYYY-MM-DD.assembled.dmult.fastq_maxee1_noN.unique.filter.good.fasta
  ```
  * Example Output:
    ```
Number of sequences post-filter: 635911
    ```

* Next we will perform one last filter to remove gaps. When sequences are removed in the previous step, some new gaps can arise in the alignment. This step will remove those.

  ```r
%%bash -s "$workDir"

cd $1

mothur "#filter.seqs(fasta=pear_merged\-YYYY\-MM\-DD.assembled.dmult.fastq_maxee1_noN.unique.filter.good.fasta, \
    processors=24, \
    vertical=T)" | head -n 50
  ```

* Finally we will print a summary of our alignment one more time.

  ```r
%%bash -s "$workDir" 

cd $1

mothur "#summary.seqs(fasta=pear_merged-YYYY-MM-DD.assembled.dmult.fastq_maxee1_noN.unique.filter.good.filter.fasta, \
    name=pear_merged-YYYY-MM-DD.assembled.dmult.fastq_maxee1_noN.names, processors=24)"
  ```
  
  * Example Output:
  ![Mothur Post Filter Seq Summary](https://cloud.githubusercontent.com/assets/7449496/12486923/c52fc84a-c032-11e5-9ac5-24e801dd8cd3.png)

#####Deunique Sequence File
* Now we can 'deunique' our sequences. This process will return our total sequence file with all those sequences that did not pass our quality control removed. 

  ```r
%%bash -s "$workDir"

cd $1

mothur "#deunique.seqs(fasta=pear_merged-YYYY-MM-DD.assembled.dmult.fastq_maxee1_noN.unique.filter.good.filter.fasta, \
    name=pear_merged-YYYY-MM-DD.assembled.dmult.fastq_maxee1_noN.good.names)" 
  ```

### Create Final QC File
* Now we are ready to create our final quality filtered file for all downstream applications.

  ```r
!cd $workDir; \
    perl -pe 's/[-.]//g if ! /^>/' \
    pear_merged-YYYY-MM-DD.assembled.dmult.fastq_maxee1_noN.redundant.fasta \
    > finalQC.fasta
  ```

  ```r
!cd $workDir; \
    mothur "#summary.seqs(fasta=finalQC.fasta, processors=24)"
  ```
  * Example Output:
  ![Final QC File Summary](https://cloud.githubusercontent.com/assets/7449496/12487650/0f190c70-c036-11e5-808d-b1133bff26ac.png)

* Print the total number of sequences in your file.

  ```r
!printf "Total number of sequences: "
!cd $workDir; grep -c ">" finalQC.fasta
  ```
  * Example Output:

    ```
Total number of sequences: 1662124
    ```

### Summary of Number of Sequences Per Samples
* Before moving on to the OTU binning steps, you might be interested in the amount of sequences you have in each of your samples. 
* We can quickly and easily summarize this information in the manners below.
* First, we can create a table with sample names, and the total number of sequences in each.

  ```r
os.chdir(workDir)

inFile = 'finalQC.summary'

re1 = re.compile('_[0-9]+$')

samp_count = dict()
with open(inFile, 'r') as inf:
    for line in inf.readlines():
        if line.startswith('seqname'):
            continue
        line = line.rstrip().split('\t')
        line[0] = re.sub(re1, '', line[0])
        try:
            samp_count[line[0]] += 1
        except KeyError:
            samp_count[line[0]] = 1
  ```
  
  ```r
# converting to dataframe
df_seq_cnt = pd.DataFrame(samp_count.items(), columns=['Sample', 'seq_count'])
df_seq_cnt
  ```
  * Example Output:
  ![Seq Table](https://cloud.githubusercontent.com/assets/7449496/12488015/99147ac6-c037-11e5-9b56-4e2e42eac479.png)

* We can also graphically display the number of sequences per sample:

  ```r
%%R -i df_seq_cnt -w 900
# plotting all

df_seq_cnt$seq_count = as.numeric(df_seq_cnt$seq_count)


df_seq_cnt = df_seq_cnt %>% 
    mutate(H2O = grepl('_F0$', Sample)) %>%
    mutate(rank = min_rank(seq_count)) %>%
    arrange(desc(rank))

df_seq_cnt$Sample = factor(df_seq_cnt$Sample, levels=df_seq_cnt$Sample)

ggplot(df_seq_cnt, aes(Sample, seq_count, fill=H2O)) +
    geom_bar(stat='identity') +
    theme(
        axis.text.x = element_text(angle=90)
    )
  ```
  * Example Output:
  ![Seq Numbers Graph](https://cloud.githubusercontent.com/assets/7449496/12488076/d666e94a-c037-11e5-9f25-f2a21927d14b.png)










