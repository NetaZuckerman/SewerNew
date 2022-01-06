
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec  7 14:13:05 2021

@author: omera
"""
import argparse
from collections import (defaultdict, Counter)
import multiprocessing as mp
from pathlib import Path
import re

import numpy as np
import pandas as pd
import pdb
import pysam


#%% Utils

def define_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        'action',
        type=str,
        choices=['bam', 'pileup', 'query_var', 'query_sam', 'move']
        )
    
    parser.add_argument(
        '-i', '--input',
        dest='input',
        type=Path,
        default='None',
        help="Input path"
        )
    
    parser.add_argument(
        '-o', '--output',
        dest='output',
        type=Path,
        default='None',
        help="Output path"
        )
    
    parser.add_argument(
        '-t', '--threads',
        dest='threads',
        type=int,
        default=1,
        help="Number of threads"
        )
    
    bam = parser.add_argument_group('BAM-specific arguments')
    bam.add_argument(
        '-r',
        '--ref',
        type=Path,
        dest="ref",
        default='None',
        help="Path to reference fasta"
        )
    
    bodek = parser.add_argument_group('Pileup arguments')
    bodek.add_argument(
        '-b',
        '--bodek',
        type=Path,
        default='None',
        help="Path to MutationsTable.xlsx"
        )
    variant = parser.add_argument_group('variant-specific arguments')
    variant.add_argument(
        '-v',
        '--variant',
        dest='query_var',
        type=str,
        default='',
        help="Name of variant, For example B.1.617.2"
        )
      
    return parser


def verify_path(args, arg_name):
    if arg_name not in args:
        raise ValueError(f"{arg_name} not provided.")
    _path = args.__dict__[arg_name]
    if not _path.exists():
        raise FileNotFoundError(f"{arg_name}: {_path} doesn't exist")
    return _path



#%% BAM

class BamConverter():
    
    def __init__(self, args):
        self.args = args

    '''
    generate pileup from bam
    for each bam file, creates pileup object and organize the data in data frame  
    input: a bam file path.
    output: dataframe with info about the positions in the input bam file.
    '''
    def _generate_pileups_from_bam(self, bam_file):
        
        with pysam.AlignmentFile(bam_file) as bamfile, pysam.FastaFile(self.args.ref) as fastafile:

            pileups = bamfile.pileup(stepper='samtools', fastafile=fastafile)
            counts = []
            pos = []
            ref = []
            # for each position in pileup (= each position in bam)
            for i, pos_pileup in enumerate(pileups):
                n = sum([pread.is_refskip for pread in pos_pileup.pileups]) 
                pos_count = Counter(map(str.upper, pos_pileup.get_query_sequences())) # count ach nucleutide depth
                pos_count['N'] = n
                pos.append(pos_pileup.reference_pos + 1)
                counts.append(pos_count)
                ref.append(fastafile.fetch(fastafile.references[0],pos_pileup.reference_pos,pos_pileup.reference_pos+1))
                
        sample_name = re.findall('(\w+).mapped.sorted.bam', bam_file.name)
        multi_index = pd.MultiIndex.from_product([sample_name, pos], names=['samplename', 'pos'])
        bamfile.close()
        fastafile.close()

        df = pd.DataFrame.from_dict(counts) \
            .set_index(multi_index) \
            .fillna(0) \
            .rename(columns={'': '-'}) \
            .astype(pd.SparseDtype(np.int16, 0))
        if '-' not in df.columns:
            df['-'] = 0 # for deletions
        #total_count = df.sum(axis=1)
        
        df['ref'] = ref
        df.set_index(['ref'],inplace = True,append = True,drop=True)
        
        #df = df.div(total_count, axis=0)
        #df['total'] = total_count #a change!
        
        return df
    
    '''
    main function of bam action.  
    with multyproccesing, send each bam file to '_generate_pileups_from_bam' func and creates united data frame. 
    input: path to dir with all bam files.
    output: one csv file named 'pileups' which contains information about the positions of all mapped bams.
    '''
    def convert(self):
        bam_path = verify_path(self.args, 'input')
        
        bam_files = bam_path.glob('*mapped.sorted.bam')
        
        
        with mp.Pool(self.args.threads) as pool:
            dfs = pool.map(self._generate_pileups_from_bam, bam_files)
        
        
        
        df = pd.concat(dfs, axis=0) \
            .melt(ignore_index=False, var_name='nucleotide', value_name='count') \
            .reset_index()
        #df = df.loc[df['count'].gt(0)] # not rm 0 - indication if mapped / not mutated.
        df['freq'] = df['count'] / df.groupby(['samplename','pos'])['count'].transform('sum') 
        df = df.dropna()
        
        out_path = self.args.output / 'pileups.csv'
        df.to_csv(out_path,index=False)

#%% Bodek-based table

class BodekMerge():
    
    def __init__(self, args):
        self.args = args
        
    '''
    merge pileup file with mutations_table
    for each pileup table, goes over the variants and merge the dataframes.
    input: pileup file path
           mutations table path
    output: dataframe with mutations frequencies in all samples for each variant.
    '''
    def _merge_pileups_bodek(self,pileup_files):

        pileup_df = pd.read_csv(pileup_files)
        bodek_dict = pd.read_excel(self.args.bodek, sheet_name=None,engine='openpyxl')

        pileup_df = pileup_df.pivot(index=['pos', 'nucleotide'], columns=['samplename'], values=['freq']) \
            .droplevel(0, axis=1) \
            .reset_index()
            
        pileup_df['pos_nuc'] = pileup_df['pos'].astype(str) + pileup_df['nucleotide']
        pileup_df = pileup_df.set_index('pos_nuc')

            
        env_cols = list(filter(lambda col: 'env' in col, pileup_df.columns))
        # for query_var extract only the asked variant
        if(self.args.action == 'query_var'):
            bodek_dict = {self.args.query_var: bodek_dict[self.args.query_var]}
            

        output_dfs = []
        # for each variant (sheet) in mutations table
        for variant, sheet in bodek_dict.items():
            sheet = sheet.set_index(sheet.Position.astype(str) + sheet['Mutation'])
            ind = pileup_df.index.intersection(sheet.index)
            sheet[env_cols] = 'NC' # if the mutations did not mapped to the reference
            sheet.loc[ind, env_cols] = pileup_df.loc[ind, env_cols]
            sheet.insert(0, 'cov_variant', variant)
            output_dfs.append(sheet)
        outputs_df = pd.concat(output_dfs, axis=0) 
        outputs_df = outputs_df.drop(["Unnamed: 9","nuc sub.1"], axis=1,errors='ignore')
        outputs_df = outputs_df.fillna(0) # if the mutation did not occured
        return outputs_df


    '''
    main function of pileup action.  
    with multyproccesing, send pilup file (or several in the same dir) to '_merge_pileups_bodek' func 
        and creates (united) data frame. 
    input: path to dir with pileup file (-or files).
           path to Mutations table
    output: one csv file named 'Monitored_Mutations.csv' of merged pileup-mutation_table
    
    
    main function of query_var action.  
    same as above, just focusing on one chosen variant from mutation table.
    '''      
    def merge(self):

        pileup_path = verify_path(self.args, 'input')
        
        pileup_files = pileup_path.glob('*pileups.csv')

        
        #for pileup_file in pileup_files:
            # print(pileup_file)
            # pileups_df = pd.read_csv(pileup_file)
            # merged_pileups = self._merge_pileups_bodek((bodek_dict, pileups_df[['samplename','pos','nucleotide','freq']]))
        
        with mp.Pool(self.args.threads) as pool:

            merged_pileups = pool.map(self._merge_pileups_bodek,pileup_files ) 
                
        merged_pileup = pd.concat(merged_pileups, axis=1) \
            .reset_index()
        merged_pileup = merged_pileup.loc[:,~merged_pileup.columns.duplicated()]

        out_path = args.output / ('%sMonitored_Mutations.csv' % (self.args.query_var))
        merged_pileup.to_csv(out_path,index=False)


#%% Pileup-based table

class PileupMerge():
    
    def __init__(self, args):
        self.args = args
        
    '''
    merge pileup file with mutations_table
    for each sample, goes over all the mutations and charactrize each mutation by the variant it belong to.
    input: pileup file path
           mutations table path
    output: dataframe of mutations indexed by variants.
    ''' 
    def _merge_pileups_bodek(self,pileup_files):
        pileup_df = pd.read_csv(pileup_files)
        bodek_dict = pd.read_excel(self.args.bodek, sheet_name=None,engine='openpyxl')

        pileup_df = pileup_df.drop(pileup_df[(pileup_df.ref == pileup_df.nucleotide)].index)
        pileup_df = pileup_df.drop(pileup_df[pileup_df.freq == 0].index)

        pileup_df = pileup_df.set_index(pileup_df.pos.astype(str) + pileup_df['nucleotide'])


        variants = []
        for variant, sheet in bodek_dict.items():
            ind = sheet.Position.astype(str) + sheet['Mutation']
            pileup_df.loc[:, variant] = 0
            pileup_df.loc[pileup_df.index.isin(ind),variant] = 1 # index 1 if mutation contained in variant
            variants.append(variant)
            #rows_to_rm = [pileup_df.loc[:,7:pileup_df.shape[1]].sum(axis=1)]

        pileup_df['variants_count'] = pileup_df[variants].sum(axis=1) 
        position = pileup_df.loc[pileup_df['variants_count'] != 0 ,'pos']
        pileup_df = pileup_df.loc[pileup_df.pos.isin(position)]
        pileup_df.drop('variants_count', inplace=True, axis=1)
        return pileup_df

    '''
    main function of query_sam action.  
    with multyproccesing, send pilup file (or several in the same dir) to '_merge_pileups_bodek' func 
        and creates (united) data frame. 
    input: path to dir with pileup file (-or files).
           path to Mutations table
    output: one csv file named 'Variants_Mutations_In_Samples.csv' of merged pileup-mutation_table, for all samples.
    
    '''  
    def merge(self):

        pileup_path = verify_path(self.args, 'input')
        
        pileup_files = pileup_path.glob('*pileups.csv')


        with mp.Pool(self.args.threads) as pool:

            merged_pileups = pool.map(self._merge_pileups_bodek,pileup_files ) 
         
        merged_pileup = pd.concat(merged_pileups, axis=0) \
            .reset_index()
        #merged_pileup = merged_pileup.loc[:,~merged_pileup.columns.duplicated()]

        out_path = args.output / ('Variants_Mutations_In_Samples.csv')
        merged_pileup.to_csv(out_path,index=False)

#%% Main
if __name__ == "__main__":

    parser = define_parser()
    args = parser.parse_args()

    if args.action == 'bam':
        bam_converter = BamConverter(args)
        bam_converter.convert()
    elif args.action == 'query_sam':
        pileup_merge = PileupMerge(args)
        pileup_merge.merge()
    elif args.action == 'pileup' or 'query_var':
        bodeq_merge = BodekMerge(args)
        bodeq_merge.merge()

    else:
        raise NotImplementedError("Those parts aren't implemented yet :(")



