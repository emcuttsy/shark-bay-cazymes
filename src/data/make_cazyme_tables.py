import pandas as pd
import re
from src.config import data_dir
import make_cazyme_tables_helpers as helpers

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

    cazyfam_df = pd.read_table(data_dir / 'interim' / 'CAZyme_family_activities.tsv', sep='\t') #reimport introduces NaNs, which are needed.
    EC_df = pd.read_table(data_dir / 'raw' / 'final_manual_annos' / 'activities.tsv')
    sub_df = pd.read_table(data_dir / 'raw' / 'final_manual_annos' / 'substrates.tsv')
    families_df = pd.read_table(data_dir / 'processed' / 'CAZyme_ct_vs_MAG.tsv', index_col = 0)
    families_ex_df = pd.read_table(data_dir / 'processed' / 'CAZyme_ct_vs_MAG_ex.tsv', index_col = 0)
    

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

    cazyfam_df['count'] = cazyfam_df['family'].apply(get_count_from_fam, args=(families_df,))
    cazyfam_df['ex_count'] = cazyfam_df['family'].apply(get_count_from_fam, args=(families_ex_df,))
    cazyfam_df = cazyfam_df.iloc[: , 1:]

    cazyfam_df.to_csv(data_dir / 'interim' / 'CAZy_fams_autoannos_all.tsv', sep='\t')
    cazyfam_df[cazyfam_df['count'] > 0].to_csv(data_dir / 'interim' / 'CAZyme_fams_autoannos_MAGs.tsv', sep='\t')

    # make final tables based on manually refined file
    cazyfam_df = pd.read_table(data_dir / 'raw' / 'final_manual_annos' / 'CAZyme_fams_refined_MAGs.tsv')
    cazyfam_df = helpers.add_substrate_metadata_cols(cazyfam_df, sub_df)
    combined_cazyfam_df = helpers.make_combined_cazyfam_df(cazyfam_df, families_df, families_ex_df)
    combined_cazyfam_df.to_csv(data_dir / 'interim' / 'CAZyme_table_interim.tsv', sep='\t')