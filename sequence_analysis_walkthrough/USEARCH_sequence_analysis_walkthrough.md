#Workflow for Sequencing Pipeline

##Downloading Links to Sequences

* First, you will receive an email from Peter Schweitzer with links to your sequences.
* Expect 4 files in Peter's email:
  * Forward Indices
  * Reverse Indicies
  * Forward Reads
  * Reverse Reads
* Using `wget`, download your files to the sequence directory on the server.
  * [Detailed Walkthrough for `wget`](./detailed_walkthroughs/wget.md)

NOTE: At this point, raw read QC analysis can be performed. Refer to the `raw_read_QC` Jupyter Notebook in the `fullCyc` dataset for these commands.

## Merge Sequences Together 

* `PEAR` is used to merge sequences from Read1 and Read2 together, while updating the quality scores of bases based on the overlapping region.
* First, we run `PEAR` in order to merge reads 1 and 2 together.
* Then we create a screed database to allow to fast access to information about the merged sequences (length, quality, etc.)
* Finally we run basic quality statistics regarding the quality scores of our sequences. 

[Detailed Walkthrough for `PEAR`](./detailed_walkthroughs/pear.md)

## Demultiplex Sequences

* At this stage, your sequences will be separated by their index reads.
* A metadata file is needed to assign sample names based on the index sequences. 
* First, screed databases are created for the two index read files.
* The merged sequence file is cycled through, replacing raw read names with sample IDs based on the metadata file loaded.
* Basic summary information such as number of sequences per sample can then be plotted.

[Detailed Walkthrough for Demultiplexing](./detailed_walkthroughs/demultiplexing.md)

## Quality Control of Sequences

* Now that we have demultiplexed our sequences, we need to run some quality control measures before downstream applications.
* First, we will remove any sequences that exceed our maximum allowed error thresholds.
* All sequences with unassigned 'N' bases are also removed.
* Remove sequences that do not align to the expected region (usually 16S or ITS regions).
* Trim any sequences that fall outside of our expected region length.
* Summarize the total number of sequences in each sample after these quality control measures.

[Detailed Walkthrough for Quality Control](./detailed_walkthroughs/qc.md)

## OTU Binning
* In this step, we will bin our quality filtered sequences into operational taxonomic units (OTUs). 
* It is recommended that all steps prior to OTU binning be performed on a per-run basis. If your full data set contains samples from more than one single run, you can concatenate your quality filtered sequence files here and proceed with OTU picking.
* First, we will reduce our sequence files to only contain unique sequences. 
* Then, we will create OTU centroids (>97% similarity) using this unique sequences file. 
* Singletons (OTUs which only appear once in the full dataset) will be removed here as well.
* OTUs clustered using `USEARCH`.
  * `USEARCH` clustering will also perform chimera checking at this stage.
* Rename sequences with OTU numbers.
* Assign taxonomy to OTUs using the SILVA database.
* Remove any non-target sequences (Chloroplast, Eukaryotic, etc.) from the dataset.
* Map reads back to the full sequence set.
* Create and OTU table and add taxonomy information to that table.

[Detailed Walthrough for OTU Binning](./detailed_walkthroughs/otu_binning.md)

