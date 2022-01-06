# SewerNew

### CoronaPipeline conda environment

* In order to install all dependencies using conda, run:
`conda env create -f <your/covid19/path>/env/CoronaPipeline.yml`
* make sure to install, in CoronaPipeline enviroment, openpyxl with: `pip install openpyxl` and pysam with `pip install pysam`
* Make sure you have samtools V 1.10 or above installed.

## Suggested Workflow Using new_sewer.py
1. Create output directory.
2. Activate conda enviroment with `conda activate CoronaPipeline`.
3. Create pileup.csv, run:   
`sewer_new.py bam -i /path/to/input/bams/dir/ -o /path/to/output/dir/ -r /path/to/reference/corona.fasta`
4. Create Monitured_mutations.csv, run:  
`sewer_new.py bam -i /path/to/input/pileup/dir/ -o /path/to/output/dir/ -b /path/to/mutationsTable.xlsx`
`sewer_new.py bam [-i INPUT_DIR_PATH] [-o OUTPUT_DIR_PATH] [-b mutationsTable.xlsx]` 
5. Optional queries:  
   * Create (variant)Monitured_mutations.csv for specific variant, run:  
`sewer_new.py bam -i /path/to/input/pileup/dir/ -o /path/to/output/dir/ -b /path/to/mutationsTable.xlsx -v variant_name`
`sewer_new.py bam [-i INPUT_DIR_PATH] [-o OUTPUT_DIR_PATH] [-b mutationsTable.xlsx] [-v VARIANT]`
   * Create Variants_Mutations_In_Samples.csv, run:  
`sewer_new.py bam -i /path/to/input/pileup/dir/ -o /path/to/output/dir/ -b /path/to/mutationsTable.xlsx`

## Additional info: new_sewer.py
`sewer_new.py bam [-i INPUT_DIR_PATH] [-o OUTPUT_DIR_PATH] [-r CORONA_REFERANCE.fasta] [optional -t NUMBER_OF_THREDS]` 
