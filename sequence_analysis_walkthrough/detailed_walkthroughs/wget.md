# Using `wget` to Download Sequence Files

### Notebook Setup

* First, create a directory for your project sequences in the `/var/seq_data/` directory on the server
  * Use the naming scheme `/var/seq_data/PROJECT_NAME/` with any subdirectories for individual sections or runs for your project
* Next, setup the `$seqDir` variable in your notebook and prepare for doanloading sequences.
   
  ```r
seqDir = '/var/seq_data/PROJECT_NAME/SUBDIRECTORY' 
  ```
   
  ```r
if not os.path.isdir(seqDir):  
     os.makedirs(seqDir)
  ```

### Download Sequence and Index Files
###### Read1
```r
!cd $seqDir; \
    wget -O read1.fq.gz \
    "LINK_TO_READ1_FILE_DOWNLOAD"
```

###### Read2
```r
!cd $seqDir; \
    wget -O read2.fq.gz \
    "LINK_TO_READ2_FILE_DOWNLOAD"
```

###### Index1
```r
!cd $seqDir; \
    wget -O index1.fq.gz \
    "LINK_TO_INDEX1_FILE_DOWNLOAD"
```

###### Index2
```r
!cd $seqDir; \
    wget -O index2.fq.gz \
    "LINK_TO_INDEX2_FILE_DOWNLOAD"
```

###Example Output
You should expect to see something similar to the following when you download each file.
```
--2015-07-02 11:21:22--  http://cbsuapps.tc.cornell.edu/Sequencing/showseqfile.aspx?mode=http&cntrl=1544200496&refid=58068
Resolving cbsuapps.tc.cornell.edu (cbsuapps.tc.cornell.edu)... 128.84.70.250
Connecting to cbsuapps.tc.cornell.edu (cbsuapps.tc.cornell.edu)|128.84.70.250|:80... connected.
HTTP request sent, awaiting response... 200 OK
Length: 1019103830 (972M) [application/octet-stream]
Saving to: `read1.fq.gz'

100%[====================================>] 1,019,103,830 51.6M/s   in 17s     

2015-07-02 11:21:40 (56.2 MB/s) - `read1.fq.gz' saved [1019103830/1019103830]
```
