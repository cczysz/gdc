This description applies to each directory. Here I will use `thca` as an example.

# Step 1: Link count files with phenotypes

The `convert.sh` file prepares the text files necessary to create the count and phenotype Robjects.

First, create 3 files

File with absolute path to each count file

`ls /group/stranger-lab/grossman_cancer_data/thca/*.gz > thca_files.txt`

File provided by Allison to link GDC UUID with TCGA UUID

`thca_htseq.txt`

File of phenotype data, 3 columns with TCGA_ID, TCGA_UUID, and Gender

Created by:

```
grep bcr_patient_barcode All_CDEs.txt | tr '\t' '\n' > barcode.txt
grep bcr_patient_uuid All_CDEs.txt | tr '\t' '\n' > uuid.txt
grep gender All_CDEs.txt | tr '\t' '\n' > gender.txt

paste barcode.txt uuid.txt gender.txt > phens.txt
```

Run the `convert.sh` script with these input files

`sh convert.sh thca_files.txt thca_htseq.txt phens.txt`

This creates the `id_conversions.txt` file read by the `expr.R` script.

# Step 2 Generate count and phenotype R objects

First, create a file with duplicated IDs

`cut -f 1 id_conversions.txt | sort | uniq -d > duplicated.txt`

This was used previously to remove individuals with both tumor and normal samples.

Run the `expr.R` script, which filters out duplicated ids, and sorts both the phen and count Robjects to have the same individuals in the same order.

# Step 3 Run the differential expression analysis

The `count_DE.R` script takes the previously generated count and phenotype Robjects and runs the `DESeq2` sex-biased differential expression analysis.
