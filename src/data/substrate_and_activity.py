import pandas as pd
import re
from src.config import data_dir
import substrate_and_activity_helpers as helpers

def make():
    
    CAZy_table = pd.read_table(data_dir / 'raw' / 'CAZyDB.07302020.fam-activities.txt', skiprows=1, names=['family', 'activity'])

    finds = []
    ECs = []
    for activity in CAZy_table['activity']:
        found = re.findall('\((EC.*?)\)',activity)
        ECs.append(found)
        finds.append("; ".join(found))

    CAZy_table['ECs'] = finds
    CAZy_table.to_csv(data_dir / 'interim' / 'CAZyme_family_activities.tsv', sep='\t')
    ECs_list = pd.DataFrame(set(sum(ECs, [])))
    ECs_list.to_csv(data_dir / 'interim' / 'CAZyme_ECs.tsv', sep='\t')

    ##############

    cazyfam_df = pd.read_table(data_dir / 'interim' / 'CAZyme_family_activities.tsv', sep='\t')
    EC_df = pd.read_table(data_dir / 'raw' / 'final_manual_annos' / 'activities.tsv')
    sub_df = pd.read_table(data_dir / 'raw' / 'final_manual_annos' / 'substrates.tsv')
    familes_df = pd.read_table(data_dir / 'processed' / 'CAZyme_ct_vs_MAG.tsv')
    families_ex_df = pd.read_table(data_dir / 'processed' / 'CAZyme_ct_vs_MAG_ex.tsv')
    

    # map substrates and activity to CAZymes
    cazyfam_df = helpers.map_substrates_and_activity(cazyfam_df, sub_df, EC_df)

    # add information from a column in carbo_df to a dataframe based on its substrates column.
    cazyfam_df = helpers.add_substrate_metadata_cols(cazyfam_df, sub_df)

    # add count info
    def get_count_from_fam(fam, familes_df):
        if fam in familes_df.columns:
            return(familes_df[fam].sum())
        else:
            return(0)

    cazyfam_df['count'] = cazyfam_df['family'].apply(get_count_from_fam, args=(familes_df,))
    cazyfam_df['ex_count'] = cazyfam_df['family'].apply(get_count_from_fam, args=(families_ex_df,))

    cazyfam_df.to_csv(data_dir / 'interim' / 'CAZy_fams_autoannos_all.tsv', sep='\t')
    cazyfam_df[cazyfam_df['count'] > 0].to_csv(data_dir / 'interim' / 'CAZy_fams_autoannos_MAGs.tsv', sep='\t')