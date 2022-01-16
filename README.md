# SewerNew

### CoronaPipeline conda environment

* In order to install all dependencies using conda, run:
`conda env create -f <your/covid19/path>/env/CoronaPipeline.yml`
* make sure to install, in CoronaPipeline enviroment, openpyxl with: `pip install openpyxl` and pysam with `pip install pysam`
* Make sure you have samtools V 1.10 or above installed.

## Suggested Workflow Using new_sewer.py
1. Activate conda enviroment with `conda activate CoronaPipeline`.
2. Create pileup.csv and Monitured_mutations.csv, run:   
`sewer_new.py -i /path/to/input/bams/dir/ -o /path/to/output/dir/ -r /path/to/reference/corona.fasta -b /path/to/mutationsTable.xlsx`
3. Create just Monitured_mutations.csv (incase pileup.csv is already exist), run:  
`sewer_new.py pileup -i /path/to/input/pileup/dir/ -o /path/to/output/dir/ -b /path/to/mutationsTable.xlsx -t number_of_threads`  
   Or if you want to execute on specific NGS files run:  
`sewer_new.py pileup -i /path/to/input/pileup/dir/ -o /path/to/output/dir/ -b /path/to/mutationsTable.xlsx -n NGS_runs`

4. Optional queries:  
   * Create (variant)Monitured_mutations.csv for specific variant (or variants), run:  
`sewer_new.py query_var -i /path/to/input/pileup/dir/ -o /path/to/output/dir/ -b /path/to/mutationsTable.xlsx -v variant_name`

   * Create Variants_Mutations_In_Samples.csv, run:  
`sewer_new.py query_freqMut -i /path/to/input/pileup/dir/ -o /path/to/output/dir/ -b /path/to/mutationsTable.xlsx -f minimum_frequency -m minimum_depth`

## Additional info: new_sewer.py
## General pipeline
`sewer_new.py [-i INPUT_DIR_PATH] [-o OUTPUT_DIR_PATH] [-r CORONA_REFERANCE.fasta] [-b mutationsTable.xlsx] [optional: -t NUMBER_OF_THREDS] [optional: -n NGS_RUNS]`

Creates 2 types of files:
### pileup.csv
Creates pileup.csv file (or files) with information about all mapped positions in input bam files. for example:
| samplename | pos | ref | alt | count | freq |
| --- | --- | --- | --- | --- | --- |
| env613 | 123 | A | G | 10 | 0.25 |
| env613 | 123 | A | A | 30 | 0.75 |
| env614 | 125 | C | G | 1 | 0.001 |

This action goes recursively over all BAM directories in the NGS_run directories from the input argument, and for each NGS_run, creates 'result' directory with pileup.csv file.
The pileup file name will be unique for each NGS run. For example: NGS130_pileup.csv file will be created for NGS130_13122021. 
### Monitored_Mutations.csv
Creates Monitored_Mutations.csv file containing all mutations in mutationsTable.xlsx (all COVID19 variants) merged with pileup.csv file (or files). for example:
| index | cov_variant | Position | Reference | Mutation | protein | variant | Mutation type | annotation | varname | nuc sub | env613 | env614 | ... | env700 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 275T | A.2.5.1 | 275 | C | T | NSP1 | L4F | SNP | Leader protein | NSP1:L4F | C275T | 0 | 0 | ... | 0.8 |
| 10747T | A.23.1 | 10747 | C | T | NSP5 | N231N | SNP_silent | 3C-like proteinase | NSP5:N231N | C10747T | 0.5 | 0.1 | ... | NC |

This action creates one Monitored_Mutations.csv table for all exsit pileup.csv files from the input argument, unless NGS runs was specified in --ngs argunent. 
The Monitored_Mutations file name will be dated. For example: Monitored_Mutations_20220112.csv will be created in date 12/01/2022 . 

### Arguments
**-i | --input : input directory (required)**  
`/path/input/` - the path to the bam files.
To create one pileup table and one Monitored_Mutations table for specific NGS run, provide a path to NGS run directory. For example: -o /data3/sewer/NGS137_07012022
To create pileup tables for all NGS runs and one Monitored_Mutations for all those NGS runs, provide a path to parent directory of all NGS runs directories. For example: -o /data3/sewer.

**-o | --output : output directory (required)**  
`/path/output` - the path to the output directory for Monitored_Mutations.csv table. created if directory do not exist.

**-r | --ref : provide refseq path (required)**  
`/path/ref.fasta` - the path to reference sequence fasta file.

**-b | --bodek : provide mutationsTable file path (required)**  
`/path/mutationsTable.xlsx` - the path to mutationsTable.xlsx.

**-t | --threads : threads number**  
`int` - requested number of threads. can be used in case you have multiple bam files. default=1.

**-n | --ngs : NGS runs**
`string` - NGS run (or runs) to focus on in generating the Monitored_Mutations.csv table. provide the NGS runs seperated by comma, for example: -n 132,133,134 (for runs NGS132_14122021, NGS133_23122021, NGS134_30122021). Make sure that the asked NGS runs directories are children directories of the input path argument.

### Commands Examples:
1. `sewer_new.py -i /data/sewer/NGS134_30122021/ -o /data/sewer/NGS134_30122021/result/ -r /data/COVID19/REF_NC_045512.2.fasta -b /data/COVID19/mutationsTable.xlsx`
This command creates result directory with NGS134_pileup.csv and Monitored_Mutations_20220112.csv files.
2. `sewer_new.py -i /data/sewer/ -o /data/sewer/result_NGS130-134/ -r /data/COVID19/REF_NC_045512.2.fasta -b /data/COVID19/mutationsTable.xlsx -n 130,131,132,134`
This command creates 'result' directories for all NGS runs in /data/sewer/ that not contain pileup file. Afterwards the pipeline will create 'result_NGS130-134' directory with Monitored_Mutations_20220112.csv table for all samples in NGS130, NGS131, NGS132, NGS133 and NGS134.

&nbsp;
## pileup action
`sewer_new.py pileup [-i INPUT_DIR_PATH] [-o OUTPUT_DIR_PATH] [-b mutationsTable.xlsx]`

**pileup : action | action to be executed (required)**  
`pileup` - creates Monitored_Mutations.csv file containing all mutations in mutationsTable.xlsx (all COVID19 variants) merged with pileup.csv file. for example:
| index | cov_variant | Position | Reference | Mutation | protein | variant | Mutation type | annotation | varname | nuc sub | env613 | env614 | ... | env700 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 275T | A.2.5.1 | 275 | C | T | NSP1 | L4F | SNP | Leader protein | NSP1:L4F | C275T | 0 | 0 | ... | 0.8 |
| 10747T | A.23.1 | 10747 | C | T | NSP5 | N231N | SNP_silent | 3C-like proteinase | NSP5:N231N | C10747T | 0.5 | 0.1 | ... | NC |

**-i | --input : input directory (required)**  
`/path/input/` - the path to the pileup file (or files).
To create one Monitored_Mutations table for specific NGS run, provide a path to NGS run directory. make sure that pileup.csv table is already exist, if not run the general pipeline instead. For example: -o /data3/sewer/NGS137_07012022
To create one Monitored_Mutations for all NGS runs, provide a path to parent directory of all NGS runs directories. For example: -o /data3/sewer.

**-o | --output : output directory (required)**  
`/path/output` - the path to the output directory. created if directory do not exist 

**-b | --bodek : provide mutationsTable file path (required)**  
`/path/mutationsTable.xlsx` - the path to mutationsTable.xlsx.

**-t | --threads : threads number**  
`int` - requested number of threads. can be used in case you have multiple pileups files. default=1.

**-n | --ngs : NGS runs**
`string` - NGS run (or runs) to focuse on in generating the Monitored_Mutations.csv table. provide the NGS runs seperated by comma, for example: -n 132,133,134 (for runs NGS132_14122021, NGS133_23122021, NGS134_30122021). Make sure that the asked NGS runs directories are children directories of the input path argument. 

### Commands Examples:
1. `sewer_new.py pileup -i /data/sewer/NGS134_30122021/ -o /data/sewer/NGS134_30122021/result/ -r /data/COVID19/REF_NC_045512.2.fasta -b /data/COVID19/mutationsTable.xlsx`
This command creates 'result' directory with Monitored_Mutations_20220112.csv file. Important to make sure that NGS134_pileup.csv file exist in /data/sewer/NGS134_30122021/ .
2. `sewer_new.py pileup -i /data/sewer/ -o /data/sewer/result_NGS130-132/ -r /data/COVID19/REF_NC_045512.2.fasta -b /data/COVID19/mutationsTable.xlsx -n 130,131,132,134`
This command creates 'result_NGS130-132' directory with Monitored_Mutations_20220112.csv table for all samples in NGS130, NGS131, NGS132. Important to make sure that NGS130_pileup.csv, NGS131_pileup.csv, NGS132_pileup.csv files exist in /data/sewer/ .

&nbsp;
## query_var action
`sewer_new.py query_var [-i INPUT_DIR_PATH] [-o OUTPUT_DIR_PATH] [-b mutationsTable.xlsx] [-v VARIANT_NAME] [optional: -t NUMBER_OF_THREDS] [optional: -n NGS_RUNS]`

**query_var : action | action to be axecuted (required)**  
`query_var` - creates (variant)Monitored_Mutations.csv file containing all mutations in one asked variant from mutationsTable.xlsx merged with pileup.csv file for example:
| cov_variant | Position | Reference | Mutation | protein | variant | Mutation type | annotation | varname | nuc sub | env613 | env614 | ... | env700 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
A.2.5.1 | 275 | C | T | NSP1 | L4F | SNP | Leader protein | NSP1:L4F | C275T | 0 | 0 | ... | 0.8 |
A.2.5.1 | 10747 | C | T | NSP5 | N231N | SNP_silent | 3C-like proteinase | NSP5:N231N | C10747T | 0.5 | 0.1 | ... | NC |

**-i | --input : input diectory (required)**  
`/path/input/` - the path to the pileup file (or files).
To create one (variant)Monitored_Mutations.csv table for specific NGS run, provide a path to NGS run directory. make sure that pileup.csv table already exist. For example: -o /data3/sewer/NGS137_07012022
To create one (variant)Monitored_Mutations.csv for all NGS runs, provide a path to parent directory of all NGS runs directories. For example: -o /data3/sewer.

**-o | --output : output directory (required)**  
`/path/output` - the path to the output directory. created if directory do not exist 

**-b | --bodek : provide mutationsTable file path (required)**  
`/path/mutationsTable.xlsx` - the path to mutationsTable.xlsx.

**-v | --variant : provide variant name (required)**  
`variant` - name of variant (or variants) in the same format as mutationsTable. provide the variant names seperated by comma, for example: -v A.2.5.1,B.1.617.2 .

**-t | --threads : threads number**  
`int` - requested number of threads. can be used in case you have multiple pileups files. default=1.

**-n | --ngs : NGS runs**
`string` - NGS run (or runs) to focus on in generating the Monitored_Mutations.csv table. provide the NGS runs separated by comma, for example: -n 132,133,134 (for runs NGS132_14122021, NGS133_23122021, NGS134_30122021). Make sure that the asked NGS runs directories are children directories of the input path argument. 

### Commands Examples:
1. `sewer_new.py query_var -i /data/sewer/NGS134_30122021/ -o /data/sewer/NGS134_30122021/result/ -r /data/COVID19/REF_NC_045512.2.fasta -b /data/COVID19/mutationsTable.xlsx -v A.2.5.1,B.1.617.2`
This command creates 'result' directory (if not exist) with A.2.5.1,B.1.617.2Monitored_Mutations_20220111.csv table for all samples in NGS134_30122021. Important to make sure that NGS134_pileup.csv file exist in /data/sewer/NGS134_30122021/ .
2. `sewer_new.py query_var -i /data/sewer/ -o /data/sewer/result_NGS130-132/ -r /data/COVID19/REF_NC_045512.2.fasta -b /data/COVID19/mutationsTable.xlsx -n 130,131,132,134`
This command creates 'result_NGS130-132' directory (if not exist) with B.1.617.2Monitored_Mutations_20220111.csv table for all samples in NGS130, NGS131, NGS132. Important to make sure that NGS130_pileup.csv, NGS131_pileup.csv, NGS132_pileup.csv files exist in /data/sewer/ .

&nbsp;
### query_freqMut action
`sewer_new.py query_sam [-i INPUT_DIR_PATH] [-o OUTPUT_DIR_PATH] [-b mutationsTable.xlsx] [optional: -t NUMBER_OF_THREDS] [optional: -n NGS_RUNS]`

**query_freqMut : action | action to be axecuted (required)**  
`query_freqMut` - creates Variants_Mutations_In_Samples.csv file containing filtered mutations by frequency from pileups.csv with association to the COVID19 variants (1=mutation exist in vaiant, 0=otherwise). for example:
| samplename | pos | ref | alt | count | freq | variants | A | A.2.5.1 | ... | PDI 353 (B.1.637 Based) |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| env613 | 275 | C | T | 5 | 0.01 | PDI 353 (B.1.637 Based) | 0 | 0 | ... | 1 |
| env614 | 10747 | C | T | 111 | 1 | A;A.2.5.1 | 1 | 1 | ... | 0 |
| env614 | 10852 | G | T | 21 | 0.5 | | 0 | 0 | ... | 0 |

**-i | --input : input directory (required)**  
`/path/input/` - the path to the pileup file (or files).
To create one Monitored_Mutations table for specific NGS run, provide a path to NGS run directory. make sure that pileup.csv table is already exist, if not run the general pipeline instead. For example: -o /data3/sewer/NGS137_07012022
To create one Monitored_Mutations for all NGS runs, provide a path to parent directory of all NGS runs directories. For example: -o /data3/sewer.

**-o | --output : output directory (required)**  
`/path/output` - the path to the output directory. created if directory do not exist 

**-b | --bodek : provide mutationsTable file path (required)**  
`/path/mutationsTable.xlsx` - the path to mutationsTable.xlsx.

**-t | --threads : threads number**  
`int` - requested number of threads. can be used in case you have multiple pileups files. default=1.

**-v | --variant : provide variant name (required)**  
`variant` - name of variant (or variants) in the same format as mutationsTable. provide the variant names seperated by comma, for example: -v A.2.5.1,B.1.617.2 .

**-n | --ngs : NGS runs**
`string` - NGS run (or runs) to focuse on in generating the Monitored_Mutations.csv table. provide the NGS runs seperated by comma, for example: -n 132,133,134 (for runs NGS132_14122021, NGS133_23122021, NGS134_30122021). Make sure that the asked NGS runs directories are children directories of the input path argument. 

**-f | --frequency : minimum frequency**
`int` - filter the mutations in the output Variants_Mutations_In_Samples table by frequency. default=0.3.  
Notice: mutation that occure in the same position as variant's mutation, will be shown in the output table even if the mutation have low frequency.

**-m | --min_depth** 
`int` - filter the mutations in the output Variants_Mutations_In_Samples table by count depth. default=10.  
Notice: mutation that occure in the same position as variant's mutation, will be shown in the output table even if the mutation have low occurence.

### Commands Examples:
1. `sewer_new.py query_freqMut -i /data/sewer/NGS134_30122021/ -o /data/sewer/NGS134_30122021/result/ -r /data/COVID19/REF_NC_045512.2.fasta -b /data/COVID19/mutationsTable.xlsx -f 0.2 -m 15`
This command creates 'result' directory (if not exist) with Variants_Mutations_In_Samples_20220111.csv file for all samples in NGS134. All the mutations in the output file will occur in frequency of at least 0.2 and in depth of atleast 15. Important to make sure that NGS134_pileup.csv file exist in /data/sewer/NGS134_30122021/ .
2. `sewer_new.py query_freqMut -i /data/sewer/ -o /data/sewer/result_NGS130-132/ -r /data/COVID19/REF_NC_045512.2.fasta -b /data/COVID19/mutationsTable.xlsx -n 130,131,132,134 -v A.2.5.1,B.1.617.2`
This command creates 'result_NGS130-132' directory with Variants_Mutations_In_Samples_20220111.csv table for all samples in NGS130, NGS131, NGS132.  All the mutations in the output file will occur in frequency of at least 0.3(=default) and in depth of atleast 10(=default). Just A.2.5.1 and B.1.617.2 will be indicted. Important to make sure that NGS130_pileup.csv, NGS131_pileup.csv, NGS132_pileup.csv files exist in /data/sewer/ .

## Aditional outputs
1. log file - creats command__'date'.log file in output directory that contains the command arguments and the NGS runs.
2. parallel to the pipeline running, prints to screen the bam files paths / pileup files paths which currently execute. 
