import pandas as pd
from pathlib import Path
from src.config import data_dir

def run():
    # polish the metadata table with small changes to the annotations
    # and reordering according to a combination of alphabetic order and phylogeny
    
    metadata = pd.read_table(data_dir / 'raw' / 'imgjgi-sb-mags-metadata.xls').iloc[: , :-1]
    
    meta_temp = metadata.copy()

    meta_temp.replace(to_replace='Candidatus Hydrogenedentes', value='Ca. Hydrogenedentes', inplace=True)
    meta_temp.replace(to_replace='Candidatus Sumerlaeota', value='BRC1', inplace=True)
    meta_temp.loc[meta_temp['Class'] == 'Gammaproteobacteria', 'Phylum'] = 'Gammaproteobacteria'
    meta_temp.loc[meta_temp['Class'] == 'Alphaproteobacteria', 'Phylum'] = 'Alphaproteobacteria'
    meta_temp.loc[meta_temp['Class'] == 'Deltaproteobacteria', 'Phylum'] = 'Myxococcota'
    meta_temp.loc[meta_temp['Class'] == 'Deltaproteobacteria', 'Class'] = 'unclassified'

    phylum_list = ['Alphaproteobacteria', 'Gammaproteobacteria', 'Myxococcota', 'Planctomycetes',
                'Chlamydiae', 'Verrucomicrobia', 'Bacteroidetes', 'BRC1',
                'Ca. Hydrogenedentes', 'Chloroflexi', 'Cyanobacteria']

    for i in ['Genus', 'Family', 'Order', 'Class']:
        meta_temp = meta_temp.sort_values(by=i, ascending=False)

    new_metadata = pd.DataFrame()
    for phylum in phylum_list:        
        phylum_df = meta_temp[meta_temp['Phylum'] == phylum]
        for i in ['Genus', 'Family', 'Order', 'Class']:
            i_unc = phylum_df[phylum_df[i] == 'unclassified']
            i_cls = phylum_df[phylum_df[i] != 'unclassified']
            phylum_df = i_cls.append(i_unc)
        new_metadata = new_metadata.append(phylum_df)
    
    new_metadata['MAG_no'] = new_metadata['Genome Name / Sample Name'].str.split('SB-').str[1]
    new_metadata['shortname'] = new_metadata['Phylum'] + ' ' + new_metadata['MAG_no']

    new_metadata.to_csv(Path(data_dir / 'processed' / 'metadata.tsv'), sep='\t')