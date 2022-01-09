# SewerNew

### CoronaPipeline conda environment

* In order to install all dependencies using conda, run:
`conda env create -f <your/covid19/path>/env/CoronaPipeline.yml`
* make sure to install, in CoronaPipeline enviroment, openpyxl with: `pip install openpyxl` and pysam with `pip install pysam`
* Make sure you have samtools V 1.10 or above installed.

## Suggested Workflow Using new_sewer.py
1. Activate conda enviroment with `conda activate CoronaPipeline`.
2. Create pileup.csv, run:   
`sewer_new.py bam -i /path/to/input/bams/dir/ -o /path/to/output/dir/ -r /path/to/reference/corona.fasta`
3. Create Monitured_mutations.csv, run:  
`sewer_new.py pileup -i /path/to/input/pileup/dir/ -o /path/to/output/dir/ -b /path/to/mutationsTable.xlsx`

4. Optional queries:  
   * Create (variant)Monitured_mutations.csv for specific variant, run:  
`sewer_new.py query_var -i /path/to/input/pileup/dir/ -o /path/to/output/dir/ -b /path/to/mutationsTable.xlsx -v variant_name`

   * Create Variants_Mutations_In_Samples.csv, run:  
`sewer_new.py query_sam -i /path/to/input/pileup/dir/ -o /path/to/output/dir/ -b /path/to/mutationsTable.xlsx`

## Additional info: new_sewer.py
### bam action
`sewer_new.py bam [-i INPUT_DIR_PATH] [-o OUTPUT_DIR_PATH] [-r CORONA_REFERANCE.fasta] [optional -t NUMBER_OF_THREDS]`

**bam : action | action to be axecuted (required)**  
`bam` - creates pileup.csv file with information on all mapped positions in input bam files. for example:  
| samplename | pos | ref | alt | count | freq |
| --- | --- | --- | --- | --- | --- | 
| env613 | 123 | A | G | 10 | 0.25 |
| env613 | 123 | A | A | 30 | 0.75 |
| env614 | 125 | C | G | 1 | 0.001 |

**-i | --input : input diectory (required)**  
`/path/input/` - the path to the bams files

**-o | --output : output directory (required)**  
`/path/output` - the path to the output directory. created if directory do not exist 

**-r | --ref : provide refseq path (required)**  
`/path/ref.fasta` - the path to reference sequence fasta file.

**-t | --threads : threads number**  
`int` - requested number of threads. can be used in case you have multiple bam files. default=1.

&nbsp;
### pileup action
`sewer_new.py pileup [-i INPUT_DIR_PATH] [-o OUTPUT_DIR_PATH] [-b mutationsTable.xlsx]` 
pileup : action | action to be axecuted (required)  
`pileup` - creates Monitored_Mutations.csv file containing all mutations in mutationsTable.xlsx (all COVID19 variants) merged with pileup.csv file. for example:
| index | cov_variant | Position | Reference | Mutation | protein | variant | Mutation type | annotation | varname | nuc sub | env613 | env614 | ... | env700 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 275T | A.2.5.1 | 275 | C | T | NSP1 | L4F | SNP | Leader protein | NSP1:L4F | C275T | 0 | 0 | ... | 0.8 |
| 10747T | A.23.1 | 10747 | C | T | NSP5 | N231N | SNP_silent | 3C-like proteinase | NSP5:N231N | C10747T | 0.5 | 0.1 | ... | NC |

**-i | --input : input diectory (required)**  
`/path/input/` - the path to the pileup file (or files).

**-o | --output : output directory (required)**  
`/path/output` - the path to the output directory. created if directory do not exist 

**-b | --bodek : provide refseq path (required)**  
`/path/mutationsTable.xlsx` - the path to mutationsTable.xlsx.

**-t | --threads : threads number**  
`int` - requested number of threads. can be used in case you have multiple pileups files. default=1.

&nbsp;
### query_var action
`sewer_new.py query_var [-i INPUT_DIR_PATH] [-o OUTPUT_DIR_PATH] [-b mutationsTable.xlsx]` 
query_var : action | action to be axecuted (required)  
`query_var` - creates (variant)Monitored_Mutations.csv file containing all mutations in one asked variant from mutationsTable.xlsx merged with pileup.csv file for example:
| index | cov_variant | Position | Reference | Mutation | protein | variant | Mutation type | annotation | varname | nuc sub | env613 | env614 | ... | env700 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 275T | A.2.5.1 | 275 | C | T | NSP1 | L4F | SNP | Leader protein | NSP1:L4F | C275T | 0 | 0 | ... | 0.8 |
| 10747T | A.23.1 | 10747 | C | T | NSP5 | N231N | SNP_silent | 3C-like proteinase | NSP5:N231N | C10747T | 0.5 | 0.1 | ... | NC |

**-i | --input : input diectory (required)**  
`/path/input/` - the path to the pileup file (or files).

**-o | --output : output directory (required)**  
`/path/output` - the path to the output directory. created if directory do not exist 

**-b | --bodek : provide refseq path (required)**  
`/path/mutationsTable.xlsx` - the path to mutationsTable.xlsx.

**-t | --threads : threads number**  
`int` - requested number of threads. can be used in case you have multiple pileups files. default=1.


&nbsp;
### query_sam action
`sewer_new.py query_sam [-i INPUT_DIR_PATH] [-o OUTPUT_DIR_PATH] [-b mutationsTable.xlsx]`
query_sam : action | action to be axecuted (required)  
`query_sam` - creates Variants_Mutations_In_Samples.csv file containing all mutations from pileups.csv with association to the COVID19 variants (1=mutation exist in vaiant, 0=otherwise). for example:
| index | samplename | pos | ref | alt | count | freq | A | A.2.5.1 | ... | PDI 353 (B.1.637 Based) |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 275T | env613 | 275 | C | T | 5 | 0.01 | 0 | 0 | ... | 1 |
| 10747T | env614 | 10747 | C | T | 111 | 1 | 1 | 1 | ... | 0 |

**-i | --input : input diectory (required)**  
`/path/input/` - the path to the pileup file (or files).

**-o | --output : output directory (required)**  
`/path/output` - the path to the output directory. created if directory do not exist 

**-b | --bodek : provide refseq path (required)**  
`/path/mutationsTable.xlsx` - the path to mutationsTable.xlsx.

**-t | --threads : threads number**  
`int` - requested number of threads. can be used in case you have multiple pileups files. default=1.
