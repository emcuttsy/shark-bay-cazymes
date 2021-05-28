# -*- coding: utf-8 -*-
#import click
import logging
from pathlib import Path
import pandas as pd
import src.config as config
import make_dataset_helpers as helpers
import make_cazyme_tables
import make_bond_and_substrate_tables
from Bio import SeqIO
import os
#from dotenv import find_dotenv, load_dotenv


#@click.command()
#@click.argument('input_filepath', type=click.Path(exists=True))
#@click.argument('output_filepath', type=click.Path())
#def main(input_filepath, output_filepath):
def main():
    """ Runs data processing scripts to turn raw data from (../raw) into
        cleaned data ready to be analyzed (saved in ../processed).
    """

    logger = logging.getLogger(__name__)
    logger.info('making final data set from raw data')

    ####################################################################
    # Make dbcAN_signalp_summary.tsv table
    ####################################################################
    
    # get the overview dataframes of dbCAN results for each MAG
    overviews = {}
    for i in Path(config.data_dir / 'raw' / 'dbCAN').iterdir():
        idd = i.stem.split('.txt')[0]
        overview = helpers.get_mag_dbcan_overview(i)
        overviews[idd] = helpers.remove_HMMer_bounds(overview)
        overviews[idd] = overview[overview['#ofTools'] > 1]
    
    # assign CAZyme families based on dbCAN data and store signalp annos
    dbcan_families = {}
    assemblies = {}
    signalp_annos = {}
    for key in overviews.keys():
        df = overviews[key]
        for index, row in df.iterrows():
            families = helpers.decide_families(row)
            if families:
                dbcan_families[row['Gene ID']] = ','.join(families)
                assemblies[row['Gene ID']] = key
                signalp_annos[row['Gene ID']] = signalp[signalp['# ID'] == row['Gene ID']]['Prediction'].tolist()[0]            
    
    # create dataframe with gene ID, taxon id, dbCAN anno, and singlP anno
    gene_df = pd.DataFrame({'dbCAN': dbcan_families, 'signalp': signalp_annos,
                            'assembly':assemblies}, index=dbcan_families.keys())

    gene_df.to_csv(Path(config.data_dir / 'processed' / 'dbCAN_signalp_summary.tsv'))


    ####################################################################
    # Make fasta file of all annotated CAZymes
    ####################################################################

    # make fasta files and get gene counts
    counts = dict.fromkeys(gene_df['assembly'], 0)
    os.system("mkdir data/processed/dbCAN_signalp_genes")

    for genome in overviews.keys():
        to_write = []
        fasta = config.data_dir / 'raw' / 'fastas' / str(genome + '.genes.faa')
        with open(fasta, 'r') as ifile, open(Path(config.data_dir / 'processed' / 'dbCAN_signalp_genes' / str('dbCAN_signalp_' + genome + '_genes.faa')), 'w') as ofile:
            for record in SeqIO.parse(ifile, 'fasta'):
                if int(record.id) in gene_df.index:
                    counts[genome] += 1
                    record.description = record.description + ' |' + genome + '|'+\
                    gene_df.loc[int(record.id)]['dbCAN'] +'|'+ gene_df.loc[int(record.id)]['signalp']
                    to_write.append(record)
            SeqIO.write(to_write, ofile, 'fasta')
    os.system("awk 1 data/processed/dbCAN_signalp_genes/*.faa > data/processed/dbCAN_signalp_genes.faa")
    

    ####################################################################
    # Make table of counts in MAG vs family
    ####################################################################

    # make table of family tallies 
    families_d = dict.fromkeys(set(gene_df['assembly']))
    excreted_d = dict.fromkeys(set(gene_df['assembly']))
    excreted_ids = gene_df[gene_df['signalp'] != 'OTHER'].index
    gene_families = set(sum(gene_df['dbCAN'].str.split(',').tolist(),[]))

    for index, row in gene_df.iterrows():
        
        if not families_d[row['assembly']]:
            families_d[row['assembly']] = dict.fromkeys(gene_families, 0)
        if not excreted_d[row['assembly']]:
            excreted_d[row['assembly']] = dict.fromkeys(gene_families, 0)
        d = families_d[row['assembly']]
        e_d = excreted_d[row['assembly']]
        dbCAN_anno = row['dbCAN'].split(',')
        
        for anno in dbCAN_anno:
            d[anno] += 1
            if index in excreted_ids:
                e_d[anno] += 1
        
    families_df = pd.DataFrame.from_dict(families_d,orient='index')
    families_df = families_df.sort_index(axis=1)
    families_df.to_csv(config.data_dir / 'processed' / 'CAZyme_ct_vs_MAG.tsv', sep='\t')

    excreted_df = pd.DataFrame.from_dict(excreted_d, orient='index')
    excreted_df = excreted_df.sort_index(axis=1)
    excreted_df.to_csv(config.data_dir / 'processed' / 'CAZyme_ct_vs_MAG_ex.tsv', sep='\t')

    
    ####################################################################
    # Make total counts table
    ####################################################################

    # make summary table with gene and family counts (excreted and not)
    summary_df = helpers.make_summary_table(gene_df, families_df, excreted_df, counts)
    summary_df.to_csv(config.data_dir / 'processed' / 'dbCAN_signalp_counts.tsv', sep = '\t')

    ####################################################################
    # Make auto-annotated CAZyme summary table
    ####################################################################

    make_cazyme_tables.make()

    ####################################################################
    # Make substrate and bond dataframes for additional annotation
    ####################################################################

    make_bond_and_substrate_tables.make() # makes the empty tables in interim used for manual annotation


if __name__ == '__main__':
    log_fmt = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    logging.basicConfig(level=logging.INFO, format=log_fmt)

    # not used in this stub but often useful for finding various files
    project_dir = Path(__file__).resolve().parents[2]
    signalp = pd.read_table(Path(config.data_dir / 'raw' / 'signalp' / 'output_protein_type.txt'), skiprows=1)

    # find .env automagically by walking up directories until it's found, then
    # load up the .env entries as environment variables
    #load_dotenv(find_dotenv())

    main()