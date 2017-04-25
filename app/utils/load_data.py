from __future__ import division
import csv
import pandas as pd
import numpy as np

def get_cummulatives():
    gene_counts_df = pd.read_csv('../data/res_finder_gene_count.csv', delimiter='\t')
    gene_counts = gene_counts_df.ix[:, gene_counts_df.columns[1:]]

    cummulative = gene_counts.sum(axis=0)
    metadata = pd.read_csv('../data/metadata.csv', delimiter='\t')
    metadata = pd.read_csv('../data/metadata.csv', delimiter='\t')

    country_sample = {}
    country_sample['absolute'] = sorted([(metadata.ix[i]['country'], cummulative.loc[str(metadata.ix[i]['sample_code'])] ) for i in metadata.index])
    country_sample['normalized'] = sorted([(metadata.ix[i]['country'], cummulative.loc[str(metadata.ix[i]['sample_code'])] / metadata.ix[i]['norm_Bacteria_pairs'] ) for i in metadata.index])

    options = ['normalized','absolute']
    cummulative_per_country = dict((option,{}) for option in options)

    for option in options:
        for x, y in country_sample[option]:
            cummulative_per_country[option].setdefault(x, []).append(y)

    return cummulative_per_country
