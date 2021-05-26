# -*- coding: utf-8 -*-
#import click
import logging
from pathlib import Path
import pandas as pd
import src.config as config
import make_dataset_helpers as helpers
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
        print(idd)
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