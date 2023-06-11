import os
import pandas as pd
#import subprocess
#import snparray_to_vcf_original as snplib
#from seqseek import Chromosome, BUILD37
#import matplotlib.pyplot as plt 
#import parq_to_vcf_snake_functions as prqlib

CHROMOSOME_INT = list(range(1, 23))
CHROMOSOME_STR = [str(chrom) for chrom in CHROMOSOME_INT]
KINDS = ['correct', 'incorrect']

#SAMPLES = ['PGD036']

"""
Set your configfile
"""
#configfile: "config.yaml"


rule all:
    input:
        #'step/PGD036_stage2_hg19_nomissing_sorted.parquet'
        #'step/PGD036_stage2_hg19_validated_correct.parquet',
        #'step/PGD036_stage2_hg19_validated_incorrect.parquet',
        #expand('step/PGD036_chr{chrom}_correct.parquet', chrom = CHROMOSOME_STR)
        #expand('step/PGD036_chr{chrom}_incorrect.parquet', chrom = CHROMOSOME_STR)
        #expand('step/PGD036_chr{chrom}_correct.vcf', chrom = CHROMOSOME_STR)
        #expand('step/PGD036_chr{chrom}_correct_imputed.vcf.gz', chrom = CHROMOSOME_STR)
        

rule cleaning:
    input:
        parquet_proc = 'data/PGD036_stage2_nofilter_hg19_processed.parquet'
    output:
        parquet_sort = 'step/PGD036_stage2_hg19_nomissing_sorted.parquet'        
    run:
        dataframe = pd.read_parquet(input.parquet_proc, engine='pyarrow')  
        dataframe = dataframe[~dataframe[['gtype_reconstructed', 'mother_gtype', 'father_gtype']].apply(lambda x: x.str.contains('NC')).any(axis=1)]
        dataframe = dataframe[dataframe['Chr'].isin(CHROMOSOME_STR)]
        dataframe['Chr'] = dataframe['Chr'].astype(int)
        dataframe.sort_values(by='Chr', ascending=True, inplace=True)
        dataframe.reset_index(inplace=True, drop=True)
        dataframe.to_parquet(output.parquet_sort, engine='pyarrow')

rule check_genotypes:
    input:
        # pesca l'output dalla regola precedente
        parquet_sort = rules.cleaning.output
    output:
        'step/PGD036_stage2_hg19_validated_correct.parquet',
        'step/PGD036_stage2_hg19_validated_incorrect.parquet'
    run:
        dataframe_cleaned = pd.read_parquet(input.parquet_sort, engine='pyarrow')

        dataframe_cleaned['validation'] = ['correct' if \
            (row['mother_gtype'] == 'AA' and row['father_gtype'] == 'AA' and row['gtype_reconstructed'] == 'AA') or \
            (row['mother_gtype'] == 'AA' and row['father_gtype'] == 'AB' and row['gtype_reconstructed'] == 'AA') or \
            (row['mother_gtype'] == 'AA' and row['father_gtype'] == 'AB' and row['gtype_reconstructed'] == 'AB') or \
            (row['mother_gtype'] == 'AA' and row['father_gtype'] == 'BB' and row['gtype_reconstructed'] == 'AB') or \
            (row['mother_gtype'] == 'AB' and row['father_gtype'] == 'AA' and row['gtype_reconstructed'] == 'AA') or \
            (row['mother_gtype'] == 'AB' and row['father_gtype'] == 'AA' and row['gtype_reconstructed'] == 'AB') or \
            (row['mother_gtype'] == 'AB' and row['father_gtype'] == 'AB' and row['gtype_reconstructed'] == 'AA') or \
            (row['mother_gtype'] == 'AB' and row['father_gtype'] == 'AB' and row['gtype_reconstructed'] == 'AB') or \
            (row['mother_gtype'] == 'AB' and row['father_gtype'] == 'AB' and row['gtype_reconstructed'] == 'BB') or \
            (row['mother_gtype'] == 'AB' and row['father_gtype'] == 'BB' and row['gtype_reconstructed'] == 'BB') or \
            (row['mother_gtype'] == 'AB' and row['father_gtype'] == 'BB' and row['gtype_reconstructed'] == 'AB') or \
            (row['mother_gtype'] == 'BB' and row['father_gtype'] == 'AA' and row['gtype_reconstructed'] == 'AB') or \
            (row['mother_gtype'] == 'BB' and row['father_gtype'] == 'AB' and row['gtype_reconstructed'] == 'AB') or \
            (row['mother_gtype'] == 'BB' and row['father_gtype'] == 'AB' and row['gtype_reconstructed'] == 'BB') or \
            (row['mother_gtype'] == 'BB' and row['father_gtype'] == 'BB' and row['gtype_reconstructed'] == 'BB') \
            else 'incorrect' for _, row in dataframe_cleaned.iterrows()]

        dataframe_valid = dataframe_cleaned[dataframe_cleaned['validation'] == 'correct']
        dataframe_notvalid = dataframe_cleaned[dataframe_cleaned['validation'] == 'incorrect']
        dataframe_valid.to_parquet('step/PGD036_stage2_hg19_validated_correct.parquet', engine='pyarrow')
        dataframe_notvalid.to_parquet('step/PGD036_stage2_hg19_validated_incorrect.parquet', engine='pyarrow')


rule filt_autosomes_valid:
    input:
        'step/PGD036_stage2_hg19_validated_correct.parquet'

    output:
        expand('step/PGD036_chr{chrom}_correct.parquet', chrom=CHROMOSOME_STR)

        #expand('step/PGD036_chr{chrom}_{type}.parquet', chrom = CHROMOSOME_STR, type = KINDS)

    run:
        parquet_valid = pd.read_parquet(input[0], engine='pyarrow')

        for chrom in CHROMOSOME_STR:

            parquet_chr = parquet_valid[parquet_valid['Chr'] == int(chrom)]
            #parquet_chr_correct = parquet_chr[parquet_chr['validation'] == 'correct']
            output_file = 'step/PGD036_chr{chrom}_correct.parquet'.format(chrom=chrom)
            parquet_chr.to_parquet(output_file, engine='pyarrow')
            #parquet_chr_correct.to_parquet(output[0], engine='pyarrow')
            
            #parquet_chr_correct.to_parquet(f'step/PGD036_chr{chrom}_{type}.parquet', engine='pyarrow')
            
            #parquet_chr_correct.to_parquet(output[0], engine='pyarrow')


rule filt_autosomes_invalid:
    input:
        'step/PGD036_stage2_hg19_validated_incorrect.parquet'

    output:
        expand('step/PGD036_chr{chrom}_incorrect.parquet', chrom=CHROMOSOME_STR)
        #expand('step/PGD036_chr{chrom}_{type}.parquet', chrom = CHROMOSOME_STR, type = KINDS)

    run:
        parquet_invalid = pd.read_parquet(input[0], engine='pyarrow')
        for chrom in CHROMOSOME_STR:
            parquet_chr = parquet_invalid[parquet_invalid['Chr'] == int(chrom)]
            #parquet_chr_incorrect = parquet_chr[parquet_chr['validation'] == 'incorrect']
            output_file = 'step/PGD036_chr{chrom}_incorrect.parquet'.format(chrom=chrom)
            parquet_chr.to_parquet(output_file, engine='pyarrow')
            
                #parquet_chr_correct.to_parquet(f'step/PGD036_chr{chrom}_{type}.parquet', engine='pyarrow')
            
                #parquet_chr_correct.to_parquet(output[0], engine='pyarrow')


rule vcf_conversion:
    input:
        'step/PGD036_chr{chrom}_correct.parquet'
    output:
        'step/PGD036_chr{chrom}_correct.vcf'
    run:
        for chrom in CHROMOSOME_INT:
        
            parq_chr = pd.read_parquet(input[0], engine='pyarrow')    
            parq_chr.rename(columns={'Chr': '#CHROM', 'Position': 'POS', 'Name': 'ID'}, inplace=True)
            parq_chr['QUAL'] = '.'
            parq_chr['FILTER'] = 'PASS'
            parq_chr['INFO'] = '.'
            parq_chr['FORMAT'] = 'GT'
            parq_chr['REF'] = parq_chr['REFALT_DBSNP'].apply(lambda x: x[0])
            parq_chr['ALT'] = parq_chr['REFALT_DBSNP'].apply(lambda x: x[1])
            parq_chr['SAMPLE'] = parq_chr['gtype_vcf'].apply(lambda x: f"{x[0]}/{x[1]}")
            vcf_sort = parq_chr.loc[:, ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'SAMPLE']]
            vcf_sort.sort_values('POS', inplace=True)
            vcf_sort.reset_index(drop=True, inplace=True)
            vcf_sort.to_csv(output[0], sep='\t', index=False)


rule beagle_gt:
    input:
        gtfile = 'step/PGD036_chr{chrom}_correct.vcf',
        mapfile = 'map/plink.chr{chrom}.GRCh37.map',
        reffile = 'reference/chr{chrom}.1kg.phase3.v5a.vcf.gz'
    
    output:
        outfile = 'step/PGD036_chr{chrom}_correct_imputed.vcf.gz'
    
    params:    
        outname = lambda wildcards, output: output.outfile.replace('.vcf.gz', '')
    
    shell:
        """
        echo 'Imputation of PGD036_chr{wildcards.chrom}_correct.vcf started...'
        java -jar beagle.22Jul22.46e.jar gt={input.gtfile} ref={input.reffile} out={params.outname} map={input.mapfile}
        
        """