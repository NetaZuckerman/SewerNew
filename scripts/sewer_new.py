
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec  7 14:13:05 2021

@author: orz
"""
import argparse
from collections import (defaultdict, Counter)
import multiprocessing as mp
from pathlib import Path
import re
import os
import numpy as np
import pandas as pd
import pysam
import glob
from datetime import datetime
import itertools
import logging
import sys
from pandas import ExcelWriter
logger = logging.getLogger(__name__)


#%% Utils

def define_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        'action',
        default = 'all',
        const='all',
        nargs='?',
        type=str,
        choices=['bam', 'pileup', 'query_var', 'query_freqMut', 'all']
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
        dest='variant',
        type=str,
        default='',
        help="Name of variant, For example B.1.617.2"
        )
    
    ngs = parser.add_argument_group('ngs-specific arguments')
    variant.add_argument(
        '-n',
        '--ngs',
        dest='ngs',
        type=str,
        default=None,
        help="Name of variant, For example B.1.617.2"
        )
    
    frq = parser.add_argument_group('Pileup arguments')
    bodek.add_argument(
        '-f',
        '--frequency',
        dest='frq',
        type=float,
        default=0.03,
        help="Path to MutationsTable.xlsx"
        )
    min_depth = parser.add_argument_group('min_depth arguments')
    variant.add_argument(
        '-m',
        '--min_depth',
        dest='min_depth',
        type=int,
        default=10,
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
        print(bam_file, flush=True)
        with pysam.AlignmentFile(bam_file) as bamfile, pysam.FastaFile(self.args.ref) as fastafile:
            
            pileups = bamfile.pileup(stepper='samtools', fastafile=fastafile)
            counts = []
            pos = []
            ref = []
            # for each position in pileup (= each position in bam)
            for i, pos_pileup in enumerate(pileups):
                n = sum([pread.is_refskip for pread in pos_pileup.pileups])
                pos_count = Counter(map(str.upper, pos_pileup.get_query_sequences()))
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
            df['-'] = 0
        #total_count = df.sum(axis=1)
        
        df['ref'] = ref
        df.set_index(['ref'],inplace = True,append = True,drop=True)
        
        #df = df.div(total_count, axis=0)
        #df['total'] = total_count #a change!
        
        return df
    
    '''
    send bam files to create pileup.  
    with multyproccesing, from /PATH/NGSxx_DATE dir send each bam file to '_generate_pileups_from_bam' func and creates united data frame. 
    input: path to dir with all bam files.
    output: one csv file named 'pileups' which contains information about the positions of all mapped bams.
    '''
    def create_output(self,path):
        ngs_dir = os.path.basename(path).split('_')[0]
        pileup_file = ngs_dir + '_pileup.csv'
        #pileup_file = os.path.join(path, pileup_file)
        listOfFiles = list()
        for (dirpath, dirnames, filenames) in os.walk(path):
            listOfFiles += [file for file in filenames]
        if pileup_file not in listOfFiles:   
            bam_files = Path(path).rglob('*mapped.sorted.bam')
            with mp.Pool(self.args.threads) as pool:
                dfs = pool.map(self._generate_pileups_from_bam, bam_files)

            df = pd.concat(dfs, axis=0) \
                .melt(ignore_index=False, var_name='alt', value_name='count') \
                .reset_index()
            #df = df.loc[df['count'].gt(0)] # not rm 0 - indication if mapped / not mutated.
            df['freq'] = df['count'] / df.groupby(['samplename','pos'])['count'].transform('sum') 
            df = df.dropna()
                
            out_dir = os.path.join(path, 'result')
            isExist = os.path.exists(out_dir)
            if not isExist: 
                os.makedirs(out_dir)
            out_path = os.path.join(out_dir, pileup_file)
                
            df.to_csv(out_path,index=False)

            
    '''
    main function of bam action.  
    sends /PATH/NGSxx_DATE directories extracted from arg.input to create_output func .
    '''                
    def convert(self):
        
        basepath = verify_path(self.args, 'input')
        print("starts creating pileup files", flush=True)
        # if args.input dir is /PATH/NGSxx_DATE directory
        if os.path.basename(basepath).startswith("NGS"):
            self.create_output(basepath) 
        # if args.input dir contains /PATH/NGSxx_DATE directories 
        else:
            for fname in os.listdir(basepath):
                path = os.path.join(basepath, fname)
                if os.path.isdir(path) and fname.startswith('NGS'):
                    self.create_output(path)
                    
        print("done creating pileup files", flush=True)


#%% Bodek-based table

class BodekMerge():
    
    def __init__(self, args):
        self.args = args
    
    def combine_indels(self,sheet,env_cols):
        df_dup = sheet[sheet.duplicated(['variant'], keep=False)]
        df_uniq = sheet[~sheet.isin(df_dup)].dropna(how = 'all')
        df_dup[env_cols] = df_dup[env_cols].apply(pd.to_numeric, errors='coerce', axis=1)
        df_new = pd.DataFrame(columns=df_dup.columns)
                
        for var in np.unique(df_dup['variant']):
            temp_df = df_dup[df_dup['variant'] == var]
            temp_df.sort_values(by=['Position'])
            p = str(temp_df.iloc[0]['Position']) + '-' + str(temp_df.iloc[-1]['Position'])
            r = ''.join(temp_df['Reference'].tolist())
            df_new = df_new.append(temp_df.iloc[0], ignore_index=True)
            df_new.at[len(df_new)-1, 'Position'] = p
            df_new.at[len(df_new)-1, 'Reference'] = r
                    
        env_cols.append('variant')   
        df_dup = df_dup[env_cols].groupby("variant").mean()
        env_cols.remove('variant')
        df_new[env_cols] = df_dup[env_cols].values
            
        sheet = pd.concat([df_new,df_uniq]).reset_index(drop=True)
        return(sheet)

    def write_in_sheets(self,merged_pileup,day):
        variants_name = self.args.variant.replace(' ','') 
        variants_name = variants_name.replace(',','_')
        variants_name = '_' + variants_name
        out_path = args.output / ('Monitored_Mutations_%s%s.xlsx' % (day,variants_name))
        writer = ExcelWriter(out_path)
        variants_names = np.unique(merged_pileup['cov_variant']) 
        for row in variants_names:
            merged_pileup[merged_pileup['cov_variant'] == row].to_excel(writer, row,index=False)
        writer.save() 
        
    '''
    merge pileup file with mutations_table
    for each pileup table, goes over the variants and merge the dataframes.
    input: pileup file path
           mutations table path
    output: dataframe with mutations frequencies in all samples for each variant.
    '''
    def _merge_pileups_bodek(self,pileup_files):
        print(pileup_files, flush=True)
        logger.info(os.path.dirname(pileup_files))
        pileup_df = pd.read_csv(pileup_files)
        bodek_dict = pd.read_excel(self.args.bodek, sheet_name=None,engine='openpyxl')
        
        # filter count up to 10 - wiil be NC
        position_count = pileup_df.groupby(['samplename','pos'],as_index = False).agg({'count':'sum'})
        position_count = position_count.loc[position_count['count'] >= self.args.min_depth ,['samplename','pos']]
        pileup_df = pd.merge(pileup_df,position_count,how='inner',left_on=['samplename','pos'],right_on=['samplename','pos']) 
        
        pileup_df = pileup_df.pivot_table(index=['pos', 'alt'], columns=['samplename'], values=['freq']) \
            .droplevel(0, axis=1) \
            .reset_index()
            
        pileup_df['pos_nuc'] = pileup_df['pos'].astype(str) + pileup_df['alt']
        pileup_df = pileup_df.set_index('pos_nuc')

            
        env_cols = list(filter(lambda col: 'env' in col, pileup_df.columns))
        # for query_var extract only the asked variant
        if(self.args.action == 'query_var'):
            list_var = self.args.variant.split(',')
            bodek_dict = {i: bodek_dict[i] for i in list_var}
            

        output_dfs = []
        # for each variant (sheet) in mutations table
        for variant, sheet in bodek_dict.items():
            sheet = sheet.set_index(sheet.Position.astype(str) + sheet['Mutation'])
            ind = pileup_df.index.intersection(sheet.index)
            sheet[env_cols] = '' # if the mutations did not mapped to the reference
            sheet.loc[ind, env_cols] = pileup_df.loc[ind, env_cols]
            sheet.insert(0, 'cov_variant', variant)
            
            if(self.args.action == 'query_var'):
                sheet = self.combine_indels(sheet,env_cols)
                
            output_dfs.append(sheet)
        outputs_df = pd.concat(output_dfs, axis=0)
        outputs_df = outputs_df.drop(["Unnamed: 9","nuc sub.1"], axis=1,errors='ignore')
        outputs_df = outputs_df.fillna('NC') # if the mutation did not occured
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
        print("starts creating Monitored_Mutations table", flush=True)
        
        #pileup_files = pileup_path.glob('*pileups.csv')

        if self.args.ngs:
            ngs_list = self.args.ngs.split(',')
            ngs_list = ['NGS' + s + '_pileup.csv' for s in ngs_list]
            pileup_files = [os.path.join(root, file) for root, dirs, files in os.walk(pileup_path) for file in files if file in ngs_list]

        else:
            pileup_files = Path(pileup_path).rglob('NGS*_pileup.csv')
    

        with mp.Pool(self.args.threads) as pool:

            merged_pileups = pool.map(self._merge_pileups_bodek,pileup_files ) 
                
        merged_pileup = pd.concat(merged_pileups, axis=1) \
            .reset_index()
        merged_pileup = merged_pileup.loc[:,~merged_pileup.columns.duplicated()]
        merged_pileup = merged_pileup.drop(columns='index')
        merged_pileup.replace('','NC',inplace=True)
        merged_pileup = merged_pileup.rename({'variant': 'aa_mut'}, axis=1)
        day = datetime.now().strftime("%Y%m%d")
        
        if(self.args.action == 'query_var'):
            self.write_in_sheets(merged_pileup,day)
        
        else:
            out_path = args.output / ('Monitored_Mutations_%s.csv' % (day))
            merged_pileup.to_csv(out_path,index=False)
        print("done creating Monitored_Mutations table", flush=True)


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
        print(pileup_files, flush=True)
        logger.info(os.path.dirname(pileup_files))
        pileup_df = pd.read_csv(pileup_files)
        bodek_dict = pd.read_excel(self.args.bodek, sheet_name=None,engine='openpyxl')

        pileup_df = pileup_df.drop(pileup_df[pileup_df['count'] == 0].index)
                
        # depyh in postion       
        position_count = pileup_df.groupby(['samplename','pos'],as_index = False).agg({'count':'sum'})
        position_count = position_count.loc[position_count['count'] >= self.args.min_depth ,['samplename','pos']]

        pileup_df = pd.merge(pileup_df,position_count,how='inner',left_on=['samplename','pos'],right_on=['samplename','pos'])   
        pileup_df = pileup_df.drop(pileup_df[(pileup_df.ref == pileup_df.alt)].index)

        # if exist mutation in the bodek and there is another mutation in the same position
        pileup_df_low = pileup_df.loc[pileup_df['freq'] < self.args.frq]
        pileup_df = pileup_df.loc[pileup_df['freq'] >= self.args.frq] 
        
        pileup_df = pileup_df.set_index(pileup_df.pos.astype(str) + pileup_df['alt'])

        variants = []
        for variant, sheet in bodek_dict.items():
            ind = sheet.Position.astype(str) + sheet['Mutation']
            pileup_df.loc[:, variant] = 0
            pileup_df.loc[pileup_df.index.isin(ind),variant] = 1
            variants.append(variant)

        pileup_df['variants_count'] = pileup_df[variants].sum(axis=1)
        position = pileup_df.loc[pileup_df['variants_count'] != 0 ,'pos'] # for enother mutation in the same position
        pileup_df_low = pileup_df_low.loc[pileup_df_low.pos.isin(position)]
        pileup_df.drop('variants_count', inplace=True, axis=1)

        pileup_df = pd.merge(pileup_df, pileup_df_low,  how='outer', left_on=['samplename','pos','ref','alt','count','freq'], right_on = ['samplename','pos','ref','alt','count','freq'])
        pileup_df = pileup_df.fillna(0)
        pileup_df.reset_index(drop=True,inplace=True)

        # join variants to on column
        pileup_df['variants'] = ''
        col = pileup_df.columns.get_loc('variants')
        for index, row in pileup_df.iterrows():
            d = itertools.compress(variants,row[6:len(row)])
            pileup_df.iloc[index,col] = ';'.join([str(i) for i in list(d)])
          
        pileup_df.insert(6, 'variants', pileup_df.pop('variants')) # replace columns
        return pileup_df

    '''
    main function of query_freqMut action.  
    with multyproccesing, send pilup file (or several in the same dir) to '_merge_pileups_bodek' func 
        and creates (united) data frame. 
    input: path to dir with pileup file (-or files).
           path to Mutations table
    output: one csv file named 'Variants_Mutations_In_Samples.csv' of merged pileup-mutation_table, for all samples.
    
    '''  
    def merge(self):

        pileup_path = verify_path(self.args, 'input')
        print("starts creating Variants_Mutations_In_Samples table", flush=True)
        
        if self.args.ngs:
            ngs_list = self.args.ngs.split(',')
            ngs_list = ['NGS' + s + '_pileup.csv' for s in ngs_list]
            pileup_files = [os.path.join(root, file) for root, dirs, files in os.walk(pileup_path) for file in files if file in ngs_list]

        else:
            pileup_files = Path(pileup_path).rglob('NGS*_pileup.csv')
            
        with mp.Pool(self.args.threads) as pool:

            merged_pileups = pool.map(self._merge_pileups_bodek,pileup_files ) 
         
        merged_pileup = pd.concat(merged_pileups, axis=0) \
            .reset_index()
        #merged_pileup = merged_pileup.loc[:,~merged_pileup.columns.duplicated()]
        merged_pileup = merged_pileup.drop(columns='index')
        merged_pileup.drop_duplicates(inplace=True)
        day = datetime.now().strftime("%Y%m%d")
        
        variants_name = ''
        if(self.args.variant != ''):
           variants_name = self.args.variant.replace(' ','') 
           variants_name = self.args.variant.replace(',','_')
           variants_name = '_' + variants_name
          
           variants = self.args.variant.split(',')
           cols_to_keep = ['samplename','pos','ref','alt','count','freq','variants'] + variants
           merged_pileup = merged_pileup[cols_to_keep]
           
        
        out_path = args.output / ('Variants_Mutations_In_Samples_%s%s.csv' % (day,variants_name))
        merged_pileup.to_csv(out_path,index=False)
        print("done creating Variants_Mutations_In_Samples table", flush=True)

#%% Main
if __name__ == "__main__":

    parser = define_parser()
    args = parser.parse_args()
    
    # Check whether the specified path exists or not
    isExist = os.path.exists(args.output)
    if not isExist: 
      # Create a new directory because it does not exist
      os.makedirs(args.output)
      
    day = datetime.now().strftime("%Y%m%d")
    out_path = args.output / ('command__%s.log' % (day)) 
    logging.basicConfig(filename = os.path.join(out_path), level = logging.INFO,
                        format = '%(asctime)s %(message)s',
                        datefmt = '%Y-%m-%d %H:%M:%S')

    logger.info(sys.argv)
    if args.action == 'all':
        bam_converter = BamConverter(args)
        bam_converter.convert()
        bodeq_merge = BodekMerge(args)
        bodeq_merge.merge()
    elif args.action == 'query_freqMut':
        pileup_merge = PileupMerge(args)
        pileup_merge.merge()
    elif args.action == 'pileup' or 'query_var':
        bodeq_merge = BodekMerge(args)
        bodeq_merge.merge()

    else:
        raise NotImplementedError("Those parts aren't implemented yet :(")



