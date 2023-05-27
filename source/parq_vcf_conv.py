#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Module for conversion of PARQUET FILE to VCF and vice versa and obtain embryo imputed correctly

__author__ = Marco Reverenna, Ivan Vogel
__copyright__ = Copyright 2022-2023
__version__ = 1.0
__maintainer__ = Marco Reverenna
__email__ = mreverenna@sund.ku.dk
__status__ = Dev
"""

import os
import pandas as pd
import subprocess
import snparray_to_vcf_original as snplib
from seqseek import Chromosome, BUILD37

sample_pq = 'PGD036'
chromosomes_int = list(range(1,23))
chromosomes_str = [str(chrom) for chrom in chromosomes_int]

def cleaning(parquet):
    """ Remove missing values, consider autosomes only, sort chromosomes in ascending order.
    Args
        parquet = raw parquet file containing all chromosomes
    """
    # exclude NC genotypes referring to mother, father, embryo
    parquet = parquet[~parquet[['gtype_reconstructed','mother_gtype','father_gtype']].apply(lambda x: x.str.contains('NC')).any(axis=1)]
    # consider only autosomes
    parquet = parquet[parquet['Chr'].isin(chromosomes_str)]
    # convert into integers
    parquet['Chr'] = parquet['Chr'].astype(int)
    # sort values in ascending referring to chromosomes
    parquet.sort_values(by='Chr', ascending=True, inplace=True)
    # reset index
    parquet.reset_index(inplace=True, drop=True)
    # save the file
    parquet.to_parquet(os.path.join(stp_dir, f'{sample_pq}_stage2_hg19_nomissing_sorted.parquet'), engine='pyarrow')
    return None

def check_genotypes():
    """ Filter genotypes considering mother, father and embryo data.
    """
    #parquet.to_parquet(os.path.join(stp_dir, f'{sample_pq}_stage2_hg19_nomissing_sorted.parquet'), engine='pyarrow')
    parquet_cleaned = pd.read_parquet(os.path.join(stp_dir, f'{sample_pq}_stage2_hg19_nomissing_sorted.parquet'), engine='pyarrow')

    parquet_cleaned['validation'] = ['correct' if \
        # checking all mother homo geneotypes AA with all possible combinations with the father                      
        (row['mother_gtype'] == 'AA' and row['father_gtype'] == 'AA' and row['gtype_reconstructed']== 'AA') or \
        (row['mother_gtype'] == 'AA' and row['father_gtype'] == 'AB' and row['gtype_reconstructed']== 'AA') or \
        (row['mother_gtype'] == 'AA' and row['father_gtype'] == 'AB' and row['gtype_reconstructed']== 'AB') or \
        (row['mother_gtype'] == 'AA' and row['father_gtype'] == 'BB' and row['gtype_reconstructed']== 'AB') or \
        
        # checking all mother het geneotypes AB with all possible combinations with the father
        (row['mother_gtype'] == 'AB' and row['father_gtype'] == 'AA' and row['gtype_reconstructed']== 'AA') or \
        (row['mother_gtype'] == 'AB' and row['father_gtype'] == 'AA' and row['gtype_reconstructed']== 'AB') or \
        (row['mother_gtype'] == 'AB' and row['father_gtype'] == 'AB' and row['gtype_reconstructed']== 'AA') or \
        (row['mother_gtype'] == 'AB' and row['father_gtype'] == 'AB' and row['gtype_reconstructed']== 'AB') or \
        (row['mother_gtype'] == 'AB' and row['father_gtype'] == 'AB' and row['gtype_reconstructed']== 'BB') or \
        (row['mother_gtype'] == 'AB' and row['father_gtype'] == 'BB' and row['gtype_reconstructed']== 'BB') or \
        (row['mother_gtype'] == 'AB' and row['father_gtype'] == 'BB' and row['gtype_reconstructed']== 'AB') or \
        
        # checking all mother homo geneotypes BB with all possible combinations with the father
        (row['mother_gtype'] == 'BB' and row['father_gtype'] == 'AA' and row['gtype_reconstructed']== 'AB') or \
        (row['mother_gtype'] == 'BB' and row['father_gtype'] == 'AB' and row['gtype_reconstructed']== 'AB') or \
        (row['mother_gtype'] == 'BB' and row['father_gtype'] == 'AB' and row['gtype_reconstructed']== 'BB') or \
        (row['mother_gtype'] == 'BB' and row['father_gtype'] == 'BB' and row['gtype_reconstructed']== 'BB') \
        
        else 'incorrect' for i, row in parquet_cleaned.iterrows()]
    # filter just genotypes validated using parental information
    parquet_valid = parquet_cleaned[parquet_cleaned['validation']=='correct']
    parquet_notvalid = parquet_cleaned[parquet_cleaned['validation']=='incorrect']
    
    # export correct embryo and uncorrect embryo reffering to parental information
    parquet_valid.to_parquet(os.path.join(stp_dir, f'{sample_pq}_stage2_hg19_validated_correct.parquet'), engine='pyarrow')
    parquet_notvalid.to_parquet(os.path.join(stp_dir, f'{sample_pq}_stage2_hg19_validated_incorrect.parquet'), engine='pyarrow')
    # you can specify or not
    return None


def filt_autosomes(str_file, kind):
    """ Filter the chromosome in parquet file
    """

    chr_dir = os.path.join(stp_dir, f'chromosomes_prq_{kind}')
    if not os.path.exists(chr_dir):
        os.makedirs(chr_dir)


    # choose which chromosome you want to filter
    for chrom in chromosomes_int:
        parquet_valid = pd.read_parquet(os.path.join(stp_dir, sample_pq + str_file), engine='pyarrow')
        # split each chromosome
        parquet_chr = parquet_valid[parquet_valid['Chr']==chrom]
        # save the file and create a directory whether not present
        
        parquet_chr_filt = parquet_chr[parquet_chr['validation']== {kind}]
        """
        chr_dir = os.path.join(stp_dir, 'chromosomes_prq')
        if not os.path.exists(chr_dir):
            os.makedirs(chr_dir)
        """
        parquet_chr_filt.to_parquet(os.path.join(stp_dir +f'chromosomes_prq_{kind}/', sample_pq + f'_chr{chrom}_{kind}.parquet'))

"""
def split_incorrect(str_file):
    # choose which chromosome you want to filter
    for chrom in chromosomes_int:
        parquet_valid = pd.read_parquet(os.path.join(stp_dir, sample_pq + str_file), engine='pyarrow')
        # split each chromosome
        parquet_chr = parquet_valid[parquet_valid['Chr']==chrom]
        
        parquet_chr.to_parquet(os.path.join(stp_dir +'chromosomes_prq_inc/', sample_pq + f'_chr{chrom}_incorrect.parquet'))
"""

def convert_vcf(kind):
    """ Convert parquet file modified into VCF format, one sample parquet file.
    Arg
        kind = string (cor if you want correct embryo genotypes VCF or inc if you want incorrect embryo genotypes VCF)
    """

    # load parquet chromosome splitted
    for chrom in chromosomes_int:
        parq_chr = pd.read_parquet(os.path.join(stp_dir + f'chromosomes_prq_{kind}/', sample_pq + f'_chr{chrom}_{kind}.parquet'), engine='pyarrow')
        parq_chr.rename(columns={'Chr': '#CHROM', 'Position': 'POS', 'Name': 'ID'}, inplace=True)
        # add VCF variables
        parq_chr['QUAL'] = '.'
        parq_chr['FILTER'] = 'PASS'
        parq_chr['INFO'] = '.'
        parq_chr['FORMAT'] = 'GT'
        parq_chr['REF'] = parq_chr['REFALT_DBSNP'].apply(lambda x: x[0])
        parq_chr['ALT'] = parq_chr['REFALT_DBSNP'].apply(lambda x: x[1])
        parq_chr['SAMPLE'] = parq_chr['gtype_vcf'].apply(lambda x: f"{x[0]}/{x[1]}")
        # reorder VCF file
        vcf = parq_chr.loc[:,['#CHROM','POS','ID','REF', 'ALT','QUAL','FILTER', 'INFO', 'FORMAT', 'SAMPLE']] # remind to select all columns

        # save the file and create a directory whether not present
        chrv_dir = os.path.join(stp_dir, f'chromosomes_vcf_{kind}')
        if not os.path.exists(chrv_dir):
            os.makedirs(chrv_dir)

        # perhaps beagle pretends the positions are in the correct order    
        vcf_sort = vcf.sort_values('POS').reset_index(drop=True)
        vcf_sort.to_csv(os.path.join(stp_dir + f'chromosomes_vcf_{kind}/', sample_pq + f'_chr{chrom}_{kind}.vcf' ), sep='\t', index=False)

"""
def beagle_imputation(chrom):

    imputed_dir = os.path.join(stp_dir, 'imputed_chromosomes')
    if not os.path.exists(imputed_dir):
        os.makedirs(imputed_dir)

    ref = f'/home/mreverenna/reference/chr{chrom}.1kg.phase3.v5a.vcf.gz '
    map = f'/home/mreverenna/map/plink.chr{chrom}.GRCh37.map '
    gt = f'/home/mreverenna/snakemake_workflow_snp_array_data/steps/chromosomes_vcf_cor/PGD036_chr{chrom}.vcf '
    out = f'/home/mreverenna/snakemake_workflow_snp_array_data/steps/imputed_chromosomes/PGD036_chr{chrom}_imputed'
    #positions = f'/home/mreverenna/analysis/vcf_experiments/masking_seeds/positions/file_diff_{chrom}.txt '

    print(f'Imputation og PGD036_chr{chrom}.vcf started...')
    command = "java -jar /home/mreverenna/programs/beagle.22Jul22.46e.jar " + f"ref={ref}" + f"gt={gt}" + f"map={map}" + f"out={out}"
    result = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()

"""


def beagle_imputation(chrom, kind='cor'):

    imputed_dir = os.path.join(stp_dir, 'imputed_chromosomes')
    if not os.path.exists(imputed_dir):
        os.makedirs(imputed_dir)

    ref = f'/home/mreverenna/reference/chr{chrom}.1kg.phase3.v5a.vcf.gz '
    map = f'/home/mreverenna/map/plink.chr{chrom}.GRCh37.map '
    gt = f'/home/mreverenna/snakemake_workflow_snp_array_data/steps/chromosomes_vcf_cor/PGD036_chr{chrom}_{kind}.vcf '
    out = f'/home/mreverenna/snakemake_workflow_snp_array_data/steps/imputed_chromosomes/PGD036_chr{chrom}_imputed'
    #positions = f'/home/mreverenna/analysis/vcf_experiments/masking_seeds/positions/file_diff_{chrom}.txt '

    print(f'Imputation og PGD036_chr{chrom}.vcf started...')
    command = "java -jar /home/mreverenna/programs/beagle.22Jul22.46e.jar " + f"ref={ref}" + f"gt={gt}" + f"map={map}" + f"out={out}"
    result = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()

def filtering():
    """ Take the same position which are resulted incorrected after the checking 
    """
    chr_vcf_imp_filt = os.path.join(stp_dir, 'chromosomes_vcf_imp_filt')
    if not os.path.exists(chr_vcf_imp_filt):
        os.makedirs(chr_vcf_imp_filt)

    chr_prq_filt = os.path.join(stp_dir, 'chromosomes_prq_filt')
    if not os.path.exists(chr_prq_filt):
        os.makedirs(chr_prq_filt)

    for chrom in chromosomes_int:
        # load imputed vcf 
        vcf_imputed = pd.read_csv(os.path.join(stp_dir + 'imputed_chromosomes', sample_pq + f'_chr{chrom}_imputed.vcf.gz'), header=8, sep='\t', dtype='object')
        prq_inc = pd.read_parquet(os.path.join(stp_dir + 'chromosomes_prq_inc', sample_pq + f'_stage2_hg19_chr{chrom}_inc.parquet'), engine = 'pyarrow')
        # remove duplicates (sometimes we have indel, rarely but could happen)
        vcf_imputed.POS.drop_duplicates(inplace=True)

        # avoid to create dataset which contain different number of snps intersecating the positions
        list_inc = prq_inc.Position.astype(str).tolist()
        vcf_imputed_filtered = vcf_imputed[vcf_imputed['POS'].isin(list_inc)]
        list_imp_filt = vcf_imputed_filtered.POS.astype(int).tolist()
        prq_inc_filt = prq_inc[prq_inc['Position'].isin(list_imp_filt)]

        # check each file contain the same number of positions
        print(vcf_imputed_filtered.shape[0], prq_inc_filt.shape[0])

        # save the files 
        prq_inc_filt.to_parquet(os.path.join(stp_dir + 'chromosomes_prq_filt/', sample_pq + f'_stage2_hg19_chr{chrom}_inc_filtered.parquet'), engine='pyarrow')
        vcf_imputed_filtered.to_csv(os.path.join(stp_dir + 'chromosomes_vcf_imp_filt/', sample_pq + f'_{chrom}_imputed_filtered.vcf.gz'), sep='\t', index=False)        


if __name__=='__main__':

    # remove warnings
    pd.options.mode.chained_assignment = None

    chromosomes_int = list(range(1,23))

    # set the paths
    res_dir = '../snakemake_workflow_snp_array_data/results/'
    out_dir = '../snakemake_workflow_snp_array_data/outputs/'
    plt_dir = '../snakemake_workflow_snp_array_data/plots/'
    stp_dir = '../snakemake_workflow_snp_array_data/steps/'
    dat_dir = '../snakemake_workflow_snp_array_data/data/'

    # import an example of vcf file to load the ref and alt allele as used in dbSNP database   
    vcf_ref = snplib.import_vcf_reference()
    vcf_ref['REFALT_DBSNP'] = vcf_ref['REF'] + vcf_ref['ALT_1']
    vcf_ref = vcf_ref.rename(columns = {'ID':'Name'})[['REFALT_DBSNP','Name']]

    # load manifest
    to_refalt = snplib.load_manifest('../snakemake_workflow_snp_array_data/data/HumanCytoSNP-12v2-1_hg19.bpm')    
    
    # load parquet file...use params['parquet'] if loading from command line argument
    df = pd.read_parquet(f'../snakemake_workflow_snp_array_data/data/{sample_pq}_stage2_nofilter_hg19.df.parquet',engine='pyarrow')

    # map SNP ID to ref and alt alleles...df['refalt']=df.Name.map(to_refalt)
    
    df = df.set_index('Name').join(to_refalt.set_index('Name'))
    df.rename(columns = {'variant/REFALT':'REFALT_MANIFEST'},inplace = True)
    
    df = df.join(vcf_ref.set_index('Name'))
        
    # convert to binary vcf format
    df['gtype_vcf'] = df.apply(lambda x: snplib.map_ab_to_bin(x['REFALT_MANIFEST'], x['REFALT_DBSNP'], x['gtype'], x['Chr'], x['Position'], BUILD37),axis=1)

    # # check if the conversion was correct by converting back to nucleotides ...
    #df['gtype_nucleotide']=df.apply(lambda x: map_to_nucleotide(x['gtype_vcf'],x['refalt'][0],x['refalt'][1]),axis=1)
    # ... and then to Illumina calls
    df['gtype_reconstructed'] = df.apply(lambda x: snplib.map_bin_to_ab(x['REFALT_MANIFEST'],x['REFALT_DBSNP'],x['gtype_vcf'],x['Chr'],x['Position'], BUILD37),axis=1)
    
    
    # check concordance
    print(df[['gtype','gtype_reconstructed']].value_counts())

    parq = df.reset_index().copy()
    #parq = df.reset_index().copy()
    parq.to_parquet(os.path.join(stp_dir, f'{sample_pq}_stage2_nofilter_hg19_processed.parquet'), engine='pyarrow')
    
    ####################
    # create directories
    print('STEP 1')
    chr_dir = os.path.join(stp_dir, 'chromosomes_prq')
    if not os.path.exists(chr_dir):
        os.makedirs(chr_dir)

    chr_dir_inc = os.path.join(stp_dir, 'chromosomes_prq_inc')
    if not os.path.exists(chr_dir_inc):
        os.makedirs(chr_dir_inc)

    # load the parquet file processed by snparray_to_vcf.py script
    parquet_file = pd.read_parquet(os.path.join(stp_dir, f'{sample_pq}_stage2_nofilter_hg19_processed.parquet'), engine='pyarrow')
    # cleaning the dataframe
    print('STEP 2')
    cleaning(parquet_file)
    # check the genotypes considering parental info
    print('STEP 3')
    check_genotypes()
    # filter each autosome
    print('STEP 4')
    filt_autosomes(str_file = '_stage2_hg19_validated_correct.parquet', kind='correct')
    filt_autosomes(str_file = '_stage2_hg19_validated_incorrect.parquet', kind ='incorrect')
    
    #!!!!!!!!
    #split_incorrect(str_file = '_stage2_hg19_validated_incorrect.parquet')
    #!!!!!!!!

    # convert from parquet to VCF
    print('STEP 5')
    convert_vcf(kind='correct')
    convert_vcf(kind='incorrect')

    print('STEP 6: imputation')
    # start the imputation
    for chrom in chromosomes_int:
        beagle_imputation(chrom, kind='cor')
    print('STEP 7')
    # filtering the final positions before obtain the imputed_gtype variable
    filtering()