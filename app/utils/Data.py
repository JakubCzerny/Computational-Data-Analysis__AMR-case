from __future__ import division
import csv
import pandas as pd
import numpy as np


class Data(object):
    def __init__(self):
        self.load()

    def load(self):
        self.gene_counts_df = pd.read_csv('../data/res_finder_gene_count.csv', delimiter='\t')
        self.metadata = pd.read_csv('../data/metadata.csv', delimiter='\t')
        self.gene_raw_df = pd.read_csv('../data/res_finder_raw_count.csv', delimiter='\t')
        self.bacteria_count_df = pd.read_csv('../data/bacteria_genera_count.csv', delimiter='\t')
        self.country_drug_use_df = pd.read_csv('../data/country_drug_use.csv', delimiter='\t')
        self.flagged_samples_df = pd.read_csv('../data/flagged_samples.csv', delimiter='\t')

    def get_cummulatives(self):
        gene_counts = self.gene_counts_df.ix[:, self.gene_counts_df.columns[1:]]
        cummulative = gene_counts.sum(axis=0)

        country_sample = {}
        country_sample['absolute'] = sorted([(self.metadata.ix[i]['country'], cummulative.loc[str(self.metadata.ix[i]['sample_code'])] ) for i in self.metadata.index])
        country_sample['normalized'] = sorted([(self.metadata.ix[i]['country'], cummulative.loc[str(self.metadata.ix[i]['sample_code'])] / self.metadata.ix[i]['norm_Bacteria_pairs'] ) for i in self.metadata.index])

        options = ['normalized','absolute']
        cummulative_per_country = dict((option,{}) for option in options)

        for option in options:
            for x, y in country_sample[option]:
                cummulative_per_country[option].setdefault(x, []).append(y)

        return cummulative_per_country
