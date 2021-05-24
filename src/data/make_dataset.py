# -*- coding: utf-8 -*-
#import click
import logging
from pathlib import Path
#import pandas as pd
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

    # get the overview dataframes of dbCAN results for each MAG
    overviews = {}
    for i in Path(config.data_dir / 'raw' / 'dbCAN').iterdir():
        try:
            overview = helpers.get_mag_dbcan_overview(i)
            print('success!')
        except:
            print(i)

        overviews[i] = helpers.remove_HMMer_bounds(overview)
        overviews[i] = overview[overview['#ofTools'] > 1] # remove entries with less than 2 tools reporting.

if __name__ == '__main__':
    log_fmt = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    logging.basicConfig(level=logging.INFO, format=log_fmt)

    # not used in this stub but often useful for finding various files
    project_dir = Path(__file__).resolve().parents[2]

    # find .env automagically by walking up directories until it's found, then
    # load up the .env entries as environment variables
    #load_dotenv(find_dotenv())

    main()