# -*- coding: utf-8 -*-
#import click
import logging
from pathlib import Path
#import pandas as pd
import src.config as config
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

    # get the overview dataframes of dbCAN results for each MAG
    overviews = {}
    for i in Path(config.data_dir / 'raw' / 'dbCAN').iterdir():
        print(i)

if __name__ == '__main__':
    log_fmt = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    logging.basicConfig(level=logging.INFO, format=log_fmt)

    # not used in this stub but often useful for finding various files
    project_dir = Path(__file__).resolve().parents[2]

    # find .env automagically by walking up directories until it's found, then
    # load up the .env entries as environment variables
    #load_dotenv(find_dotenv())

    main()

######################################################################################
# Functions used in data processing
######################################################################################

def get_mag_dbcan_overview(mag_dir):
    # locate the overview file in the dbCAN results for a MAG
    overview_file = mag_dir / 'overview.txt'
    overview = pd.read_table(overview_file)
    return overview


def remove_HMMer_bounds(df):
    # remove the start and stop position info from HMMER results
    df['HMMER'] = df['HMMER'].str.split('(', expand=True)
    df['HMMER'] = df['HMMER'].str.split('_', expand=True)
    return df


def decide_families(row):
    # assign CAZyme families to a gene based on its dbCAN results
    fams = [row['HMMER'], row['Hotpep'], row['DIAMOND']]

    GH_count = 0
    for i in fams:
        if 'GH' in i:
            GH_count += 1

    tally = list(fams)
    for c in fams:
        if c == 'N':
            tally.remove(c)
        if '+' in c:
            tally.remove(c)
            tally = tally + c.split('+')
                        
    families = []
    for f in tally:
        if tally.count(f) > 1:
            families.append(f)
    families = set(families)

    if not families and GH_count > 1:
        families = set(['GH0']) # add GH0 for ambiguous GH family. 
    if not len(families) == 0:
        return families
    

def get_family_counts(overview_df, genome):
    fam_counts = {}
    for i, row in overview_df[genome].iterrows():
        families = decide_families(row)
        if families:
            for fam in families:
                if fam not in fam_counts.keys():
                    fam_counts[fam] = 1
                else:
                    fam_counts[fam] += 1
    return fam_counts