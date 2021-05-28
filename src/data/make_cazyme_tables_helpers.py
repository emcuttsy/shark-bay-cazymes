import pandas as pd
import numpy as np
import matplotlib.colors
import matplotlib.pyplot as plt
import seaborn as sns

def get_substrate_info(substrate_string, colname, carbo_df):
    """Get values in a column of the carbohydrates spreadsheet based on a string-list of substrates.
    
    Parameters:
    substrate_string (str): list of substrates represented as a string with each value separated by "; "
    colname (str): name of the column in the carbohydrates spreadsheet to access
    carbo_df: dataframe of carbohydrates

    Returns:
    str: "; "-separated set (no repeats) of items in the column specified by colname for the rows specified by substrate_string
    float: np.nan is returned if th
    """

    if not pd.isna(substrate_string):
        substrates = substrate_string.split('; ')
        substrate_info = []
        for substrate in substrates:
            info = carbo_df[carbo_df['name'] == substrate][colname].values[0]
            if not pd.isna(info):
                info = info.split('; ')
                substrate_info += info
        if substrate_info:
            substrate_info = '; '.join(set(substrate_info))
    else:
        substrate_info = np.nan
    return(substrate_info)


def get_substrates_and_activities_from_ECs(ECs, EC_df):
    """Get substrates and activities from a "; "-separated string list of EC numbers

    Parameters:
    ECs (str): list of ECs represented as a string with each value separated by "; "
    EC_df (DataFrame): table of EC values with substrates and activities

    Returns:
    tuple: substrates and activities, each either a "; "-separated string list or np.nan if no substrate/activity was available
    """
    substrates = []
    activities = []

    if not pd.isna(ECs):
        ECs = ECs.split('; ')
        for EC in ECs:
            EC_substrates = EC_df.loc[EC_df['EC'] == EC]['substrates'].values[0]
            EC_activities = EC_df.loc[EC_df['EC'] == EC]['activity'].values[0]
            if not pd.isna(EC_substrates):
                EC_substrates = EC_substrates.split(', ')
                substrates += EC_substrates
            if not pd.isna(EC_activities):
                EC_activities = EC_activities.split('; ')
                activities += EC_activities
                
    if not substrates:
        substrates = np.nan
    else:
        substrates = '; '.join(set(substrates))
    if not activities:
        activities = np.nan
    else:
        activities = '; '.join(set(activities))

    return(substrates, activities)

def get_cazyme_summary(family, df, meta_df):
    # returns count, counts df rows not equal to 0, and metadata rows of taxa with at least 1 gene
    return df[df[family]!=0], meta_df[meta_df['taxon_oid'].isin(df[df['CE11'] != 0].index)]




## DEALING WITH GH CATEGORIES
## Used for working with the manually-made bond and substrate category

def get_category(row, target_col, split_glucosidase = False):
    col_order = ['cat3', 'cat2', 'cat1']
    target_i = col_order.index(target_col)
    to_return = None
    if (row['cat1'] == 'other glycosidase'):
        to_return = row['cat1']
    if not pd.isna(row[target_col]):
        to_return = row[target_col]
    else:
        for i in col_order[target_i+1::]:
            if not pd.isna(row[i]) and 'name origin' in row.index:
                if row['name origin'] != 'polymer':
                    to_return = 'other ' + row[i]
            elif not pd.isna(row[i]):
                to_return = row[i]
    if split_glucosidase == True and to_return == 'glucosidase':
        to_return = row['cat2']
    if to_return:
        return(to_return)

def get_cat_list(catcol, bonds_df, split_glucosidase = False):
    catlist = []
    for index, row in bonds_df.iterrows():
        catlist.append(get_category(row, catcol, split_glucosidase=split_glucosidase))
    return(catlist)

def get_categories_from_activity_string(string, cat, bonds_df, return_strlist=False):
    if pd.isna(string):
        return
    else:
        spl = string.split('; ')
        key = bonds_df.columns[0]
        if key == 'bond':
            activities = []
            for i in spl:
                activities.append(' '.join(i.split(' ')[1:4]))
        else:
            activities = spl
        bonds = bonds_df[bonds_df[key].apply(lambda x: True if any(i in x for i in activities) else False)]
        categories = []
        for index, row in bonds.iterrows():
            categories.append(get_category(row, cat))
        categories = list(set(categories))
        categories = [x for x in categories if x != None]
        if not return_strlist:
            return(list(set(categories)))
        else:
            return('; '.join(list(set(categories))))


def get_cat_dict(exclude_length, bonds_df, combined_df, split_glucosidase = True):
    # get category list
    cat = 'cat1'
    cat_list = set(get_cat_list(cat, bonds_df, split_glucosidase=split_glucosidase))
    cat_dict = dict.fromkeys(cat_list)

    # get target column, either bond activity or substrates
    if bonds_df.columns[0] == 'bond':
        target_col = 'bond activity'
    else:
        target_col = 'substrates'

    # drop rows of combined_ghpl_df for which the activities column is nan
    no_nans = combined_df[combined_df[target_col].notna()]

    # loop through categories in cat list
    for c in cat_dict.keys():
        key = bonds_df.columns[0]
        # get the bond types associated with the category:
        if split_glucosidase == True and (c == 'alpha glucosidase' or c == 'beta glucosidase'):
            associated_bonds = bonds_df[bonds_df['cat2'] == c][key]
        else:
            associated_bonds = bonds_df[bonds_df[cat] == c][key]
        # get genes with an activity associated in that category
        genes = no_nans.loc[no_nans[target_col].apply(lambda x: True if any(i in x for i in associated_bonds) else False)]
        # exclude genes with activities in more categories than the cutoff
        genes = genes[genes[target_col].apply(get_categories_from_activity_string, args=(cat,bonds_df)).str.len() <= exclude_length]
        # add the sum of the counts column
        count = genes['count'].sum()
        num_mags = len(set(genes['accession']))
        cat_dict[c] = [count, num_mags]
    # remove none from cat dict
    cat_dict = {k:v for k,v in cat_dict.items() if v is not None}
    return(cat_dict)


def get_cazyfams_from_substrate(substrate, df, substratecol):
    substratestr = set([substrate])
    issubstratestr = substratestr.issubset
    return(df[df[substratecol].apply(get_strlist).map(issubstratestr)]['family'])



## AUTO-MAPPING
## Map substrates and activities to cazyme families based on carbohydrates.txt and CAZyme_EC_activities.txt


def map_substrates_and_activity(cazyfam_df, carbo_df, EC_df):
    # maps substrates and activity to cazyfam_df based on ECs. 
    substrates_list = []
    activities_list = []

    for index, row in cazyfam_df.iterrows():
        substrates, activities = get_substrates_and_activities_from_ECs(row['ECs'], EC_df)
        substrates_list.append(substrates)
        activities_list.append(activities)

    cazyfam_df['substrates'] = substrates_list
    cazyfam_df['bond activity'] = activities_list
    return(cazyfam_df)


def add_substrate_metadata_cols(df, carbo_df, colnames=['superfamily', 'family', 'subfamily', 'size']):
    # adds data from list of substrate metadata columns to given df
    if isinstance(colnames, str):
        colnames = [colnames]
    for name in colnames:
        df['substrate_' + name] = df['substrates'].apply(get_substrate_info, args=(name, carbo_df))
    return(df)


def make_cazyfam_rows(fam, magid, cazyfam_df, cazyme_ex_df):
    cazyfam_row = cazyfam_df[cazyfam_df['family'] == fam.name].copy()
    cazyfam_row['count'] = fam.tolist()[0]
    cazyfam_row['ex_count'] = cazyme_ex_df.loc[magid][fam.name]
    cazyfam_row['accession'] = magid
    return(cazyfam_row)
    

def make_sub_df_from_row(row, cazyfam_df, cazyme_ex_df):
    accession = row.name
    nonzero = row[row !=0]
    returns = nonzero.to_frame(0).apply(make_cazyfam_rows, axis=1, args=(row.name,cazyfam_df, cazyme_ex_df,))
    return(pd.concat(returns.tolist(), axis=0))


def make_combined_cazyfam_df(cazyfam_df, cazyme_df, cazyme_ex_df):
    return(pd.concat(cazyme_df.apply(make_sub_df_from_row, axis=1, args=(cazyfam_df, cazyme_ex_df)).tolist(), axis=0))


## HELPERS

def get_strlist(strlist):
    # takes a '; '-separated list as a string and returns the list of strings
    if isinstance(strlist, str):
        strs = strlist.split('; ')
    else:
        strs = []
    return strs


def filter_by_len(ls, exclude_len):
    # returns empty list if length of list is over exclude length
    if isinstance(ls, list):
        if len(ls) > exclude_len and len(ls) != 0:
            return([])
        else:
            return(ls)
    else:
        return(ls)

def get_list_from_col(df, col, count_col='count', exclude_len=None):
    # returns list with each item appearing the number of times in the count column
    strlists = df[col]
    col = strlists.apply(get_strlist)*df[count_col]
    if exclude_len != None:
        col = col.apply(filter_by_len, args=(exclude_len,))
    col = col.tolist()
    col = sum(col, [])
    return col


def get_cazyfams_from_substrate(substrate, df, substratecol):
    substratestr = set([substrate])
    issubstratestr = substratestr.issubset
    return(df[df[substratecol].apply(get_strlist).map(issubstratestr)]['family'])


def get_rows_from_substrate(substrate, df):
    return(df[df['substrates'].apply(get_strlist).apply(lambda x: substrate in x)])


def get_rows_from_substrate_list(substrate_list, df):
    return(df[df['substrates'].apply(get_strlist).apply(lambda x: any(item in x for item in substrate_list))])