import pandas as pd
from src.config import data_dir
metadata = pd.read_table(data_dir / 'raw' / 'imgjgi-sb-mags-metadata.xls')

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


def make_summary_table(gene_df, families_df, excreted_df, counts):

    # make summary table with gene and family counts (excreted and not), and percent excreted
    summary_d = dict.fromkeys(set(gene_df['assembly']))
    family_counts = families_df.astype(bool).sum(axis=1)
    secreted_counts = excreted_df.sum(axis=1)
    secreted_family_counts = excreted_df.astype(bool).sum(axis=1)

    for key in summary_d.keys():
        count = counts[key]
        families = family_counts[key]
        secreted = secreted_counts[key]
        secreted_families = secreted_family_counts[key]
        summary_d[key] = {'ct': count, 'fams': families, 
                        'ex_ct':secreted, 'ex_fams': secreted_families}
        
    summary_df = pd.DataFrame.from_dict(summary_d, orient='index')
    
    phyla = [metadata[metadata['taxon_oid'] == int(idd)]['Phylum'].tolist()[0] for idd in summary_df.index]
    classes = [metadata[metadata['taxon_oid'] == int(idd)]['Class'].tolist()[0] for idd in summary_df.index]
    orders = [metadata[metadata['taxon_oid'] == int(idd)]['Order'].tolist()[0] for idd in summary_df.index]
    families = [metadata[metadata['taxon_oid'] == int(idd)]['Family'].tolist()[0] for idd in summary_df.index]
    genera = [metadata[metadata['taxon_oid'] == int(idd)]['Genus'].tolist()[0] for idd in summary_df.index]
    MAG = [metadata[metadata['taxon_oid'] == int(idd)]['Genome Name / Sample Name'].tolist()[0] for idd in summary_df.index]

    temp = []
    for i in range(len(phyla)):
        p = phyla[i]
        if p.startswith('Candidatus'):
            p = p.split('Candidatus ')[1]
        if p == 'Proteobacteria':
            p = classes[i]
        if p == 'Deltaproteobacteria':
            p = 'Myxococcota'
        temp.append(p)
    phyla=temp

    temp = []
    for i in range(len(classes)):
        p = phyla[i]
        if p == 'Deltaproteobacteria':
            p = 'Myxococcota'
        temp.append(p)
    classes=temp
        
    summary_df['MAG'] = MAG
    summary_df['phylum'] = phyla
    summary_df['class'] = classes
    summary_df['order'] = orders
    summary_df['family'] = families
    summary_df['genus'] = genera

    return summary_df
