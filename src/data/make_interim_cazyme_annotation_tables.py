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
    CAZy_table.to_csv(data_dir / 'interim' / 'CAZyme_fam_vs_activity', sep='\t')
    ECs_list = pd.DataFrame(set(sum(ECs, [])))
    ECs_list.to_csv(data_dir / 'interim' / 'CAZyme_ECs.tsv', sep='\t')
    