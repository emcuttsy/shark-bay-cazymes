import pandas as pd
import make_cazyme_tables_helpers as helpers
from src.config import data_dir
import collections

def make():

    cazyfam_df = pd.read_table(data_dir / 'interim' / 'CAZyme_table_interim.tsv', index_col = 0)
    sub_df = pd.read_table(data_dir / 'raw' / 'final_manual_annos' / 'substrates.tsv')
    families_df = pd.read_table(data_dir / 'processed' / 'CAZyme_ct_vs_MAG.tsv', index_col = 0)
    families_ex_df = pd.read_table(data_dir / 'processed' / 'CAZyme_ct_vs_MAG_ex.tsv', index_col = 0)


    ghpl_df = cazyfam_df[(cazyfam_df['family'].str.contains('GH') | cazyfam_df['family'].str.contains('PL'))]
    substrates = collections.Counter(helpers.get_list_from_col(ghpl_df, 'substrates', count_col='count'))
    ex_substrates = collections.Counter(helpers.get_list_from_col(ghpl_df, 'substrates', count_col='ex_count'))

    # add in missing (0 count) substrates from ex_counts
    for key in substrates.keys():
        if key not in ex_substrates.keys():
            ex_substrates[key] = 0
            
    # add make substrates dataframe and export

    substrates = dict(sorted(substrates.items(), key=lambda item: item[1], reverse=True))
    ex_substrates = dict(sorted(ex_substrates.items(), key=lambda item: item[1], reverse=True))
    sb_substrates = pd.DataFrame({'count': substrates.values(), 'ex_count': ex_substrates.values()}, index=substrates.keys())
    sb_substrates['delta'] = sb_substrates['count'] - sb_substrates['ex_count']

    mags = []
    ex_mags = []
    for substrate in substrates.keys(): 
        families = helpers.get_cazyfams_from_substrate(substrate, ghpl_df, 'substrates')
        magids = []
        ex_magids = []
        for i in families:
            if i in families_df.columns:
                magids += list(families_df[families_df[i] != 0].index)
            if i in families_ex_df.columns:
                ex_magids += list(families_ex_df[families_ex_df[i] != 0].index)
        mags.append(len(set(magids)))
        ex_mags.append(len(set(ex_magids)))

    sb_substrates['mags'] = mags
    sb_substrates['ex_mags'] = ex_mags
    sb_substrates['superfamily'] = pd.Series(sb_substrates.index).apply(helpers.get_substrate_info, args=('superfamily',sub_df,)).tolist()
    sb_substrates['family'] = pd.Series(sb_substrates.index).apply(helpers.get_substrate_info, args=('family',sub_df,)).tolist()
    sb_substrates['subfamily'] = pd.Series(sb_substrates.index).apply(helpers.get_substrate_info, args=('subfamily',sub_df,)).tolist()
    sb_substrates.to_csv( data_dir / 'interim' / 'substrates_inMAGs.tsv', sep='\t')



    # make bonds dataframe and export

    sb_activities = ghpl_df['bond activity'].str.split(' ')
    sb_activities_df = pd.DataFrame({1:sb_activities.str[0], 2:sb_activities.str[1],
                                    3:sb_activities.str[2], 4:sb_activities.str[3]})

    sb_activities_df = sb_activities_df[sb_activities_df[3].notna()]
    sb_activities_df[4] = sb_activities_df[4].str.split(';').str[0]
    sb_activities_df = sb_activities_df[sb_activities_df[3] != 'to']
    sb_bonds = set(sb_activities_df[2] + " " + sb_activities_df[3] + " " + sb_activities_df[4])

    ex_cazyfam_df = cazyfam_df[cazyfam_df['ex_count'] != 0]

    combined_mag_ghpl_df = cazyfam_df[cazyfam_df['family'].str.contains('GH') | cazyfam_df['family'].str.contains('PL')]
    ex_combined_mag_ghpl_df = combined_mag_ghpl_df[combined_mag_ghpl_df['ex_count'] != 0]
    mags = []
    ex_mags = []
    count = [] # COUNT HAS A PROBLEM 
    ex_count = []
    sb_bonds = [x for x in sb_bonds if str(x) != 'nan']
    for i in sb_bonds:
        rows_with_bond = combined_mag_ghpl_df[combined_mag_ghpl_df['bond activity'].str.contains(i, regex=False) == True]
        ex_rows_with_bond = ex_combined_mag_ghpl_df[ex_combined_mag_ghpl_df['bond activity'].str.contains(i, regex=False) == True]
        mags.append(len(set(rows_with_bond['accession'])))
        ex_mags.append(len(set(ex_rows_with_bond['accession'])))
        count.append(rows_with_bond['count'].sum())
        ex_count.append(ex_rows_with_bond['ex_count'].sum())

    sb_bonds_df = pd.DataFrame({'count': count, 'ex_count': ex_count, 'mags': mags, 'ex_mags': ex_mags}, index=sb_bonds)
    sb_bonds_df.to_csv(data_dir / 'interim' / 'bond_targets_inMAGs.tsv', sep='\t')


def reimport_and_update_cazyme_table():
    
    cazyfam_df = pd.read_table(data_dir / 'interim' / 'CAZyme_table_interim.tsv', index_col = 0)

    sb_bonds_df = pd.read_table(data_dir / 'raw' / 'final_manual_annos' / 'bond_targets_inMAGs_categories.tsv')
    sb_subs_df = pd.read_table(data_dir / 'raw' / 'final_manual_annos' / 'substrate_inMAGs_categories.tsv')

    cazyfam_df['substrate cat1'] = cazyfam_df['substrates'].apply(helpers.get_categories_from_activity_string, args=('cat1',sb_subs_df, ), return_strlist=True)
    cazyfam_df['substrate cat2'] = cazyfam_df['substrates'].apply(helpers.get_categories_from_activity_string, args=('cat2',sb_subs_df, ), return_strlist=True)
    cazyfam_df['activity cat1'] = cazyfam_df['bond activity'].apply(helpers.get_categories_from_activity_string, args=('cat1',sb_bonds_df, ), return_strlist=True)
    cazyfam_df['activity cat2'] = cazyfam_df['bond activity'].apply(helpers.get_categories_from_activity_string, args=('cat2',sb_bonds_df, ), return_strlist=True)

    cazyfam_df.to_csv(data_dir / 'processed' / 'CAZyme_table.tsv', sep='\t')