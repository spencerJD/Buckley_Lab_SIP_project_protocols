# Demultiplexing Sequence Reads by Sample ID

This walkthrough will lead you through demultiplexing your sequences by index read and the creation of some basic summary information regarding the number of sequences in each sample.

### Initializing Your Jupyter Notebook
* Start by importing all libraries you will need for all steps involved in demultiplexing.

    ```python
import screed
import pandas as pd
from collections import Counter
import matplotlib.pyplot as plt
import numpy as np
from mpltools import style
import mpld3

style.use("ggplot")

%matplotlib inline
    ```

### Load Your Index Read Files
* Next, set variables to your index read files.

    ```python
ir1 = "/var/seq_data/PROJECT_NAME/index_read1.fastq.gz"
ir2 = "/var/seq_data/PROJECT_NAME/index_read2.fastq.gz"
    ```

* Create screed databases for these index read files.

    ```python
screed.read_fastq_sequences(ir1)
screed.read_fastq_sequences(ir2)

ir1db = screed.ScreedDB(ir1)
ir2db = screed.ScreedDB(ir2)
    ```

### Load Mapping File
* The mapping file is a critical component for demultiplexing your sequences.
* This file contains the sample names assigned to each index read, allowing our demultiplexing script to rename sequences according to the sample that they came from.
* Additional metadata describing each sample should be included in this file, for downstream application in the `phyloseq` R package.
* The format of this file is critical, with an example shown [here](INSERT_LINK_HERE).


  ```python
df_map = pd.read_csv("data/Mapping_file.txt", delimiter="\t", usecols=[0,1,2], names=["SampleID","i1","i2"])
df_map
  ```
  ![Example Table Output](https://cloud.githubusercontent.com/assets/7449496/12459259/d4d1a1cc-bf7b-11e5-88fa-4730b005a140.png)

  * NOTE: Using this code, the first three columns of your metadata file must be the SampleID, index read 1 sequence, and index read 2 sequence.

* Next, you will make a dictionary containing each Sample ID, and a concatenated string of the two index reads.

  ```python
map_d = dict(map(lambda x: (x[0] + x[1], x[2]), zip(df_map.i1, df_map.i2, df_map.SampleID)))
  ```

* It is a good idea to double check that this dictionary is defined exactly as it should be.

  ```python
print map_d.keys()[10], map_d[map_d.keys()[10]]
cgagcgactcgacgag PTR2
  ```

  * NOTE: You can also enter `print map_d` in order to check the entire dictionary for accuracy.

### Demultiplex Samples

* First, we will read in the screed database for our `PEAR` merged sequences.

  ```python
screed.read_fastq_sequences("/var/seq_data/PROJECT_NAME/pear_merged.assembled.fastq")
db = screed.ScreedDB("/var/seq_data/PROJECT_NAME/pear_merged.assembled.fastq_screed")
  ```

* Next, we will demultiplex by renaming sequences based on their index reads.

  ```python
counter = 0
unassigned = 0
cnt = Counter()
fn = "/var/seq_data/PROJECT_NAME/pear_merged_raw.assembled.dmult.fastq"
with open(fn, "w") as out:
    for rec in db.itervalues():
        index_concat = ir1db[rec["name"]]["sequence"] + ir2db[rec["name"]]["sequence"]
        try:
            new_name = map_d[index_concat.lower()] + "_" + str(counter) #case sensitive
        except KeyError:
            unassigned += 1
            continue
        counter += 1
        cnt[map_d[index_concat.lower()]] += 1 #case sensitive
        s, q = rec["sequence"], rec["accuracy"]
        out.write("@%s orig_name=%s\n%s\n+\n%s\n"%(new_name,rec["name"],s,q))
print "Unable to assign %s reads to samples"%unassigned
  ```

  * NOTE: Depending on the plate submitted for sequencing (especially if it was a shared plate) there will be a large number of sequences unassigned in this step. We are only concerned with the sequences assigned to your samples.

* Now we can print a quick summary graph showing the number of sequences in each sample

  ```python
X = np.arange(len(cnt.keys()))
Y = [v for k, v in cnt.most_common()]
fig = plt.figure()
fig.set_size_inches((18,8))
ax = fig.add_subplot(111)
bars = ax.bar(X, Y, width=0.50)
xlims = ax.set_xlim((X.min(), X.max()))
xt = ax.set_xticks(X)
xtl = ax.set_xticklabels([k for k, v in cnt.most_common()], rotation=90, ha="center", size=6)
  ```
  
  ![Sequence Counts](https://cloud.githubusercontent.com/assets/7449496/12461210/addc76e6-bf85-11e5-8378-1ff0b3147e86.png)









