# python helper functions

def sort_dict(d, reverse=False):
    return(dict(sorted(d.items(), key=lambda item: item[1], reverse=reverse)))

def list_series_contains(s, search_string):
    string = set([search_string])
    isstring = string.issubset
    return (s.map(isstring))