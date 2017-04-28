from __future__ import division
import csv
import pandas as pd
import numpy as np
from sklearn.preprocessing import MaxAbsScaler

class Data(object):
    def __init__(self):
        self.load()


    def load(self):
        self.gene_counts_df = pd.read_csv('../data/res_finder_gene_count.csv', delimiter='\t')
        self.metadata = pd.read_csv('../data/metadata.csv', delimiter='\t')
        # self.metadata.set_index('sample_code')
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


    def get_gene_count_with_drugs(self, cut_off=None, normalized=True):
        columns = [col for col in self.country_drug_use_df if 'ldu' in col]
        drugs_per_country = dict((self.country_drug_use_df['country_2letter'].ix[i], self.country_drug_use_df[columns].ix[i]) for i in range(9))
        sample_country = dict((self.metadata.ix[i]['sample_code'], self.metadata.ix[i]['country']) for i in range(self.metadata.shape[0]))

        X = np.array([drugs_per_country[sample_country[int(code)]].values for code in self.gene_counts_df.columns[1:]])
        Y = self.gene_counts_df.ix[:, self.gene_counts_df.columns[1:]].T

        if normalized:
            meta = self.metadata.ix[:,['sample_code','norm_Bacteria_pairs']]
            meta = meta.set_index('sample_code')

            # Y = np.array([Y.ix[str(code)].values / meta.ix[int(code)].ix['norm_Bacteria_pairs'] for code in Y.index])
            for code in Y.index:
                Y.ix[code] = Y.ix[code].apply(lambda x: np.divide(float(x),meta.ix[int(code),'norm_Bacteria_pairs']))

            # Y = Y * 1000000
            scaler = MaxAbsScaler()
            Y = scaler.fit_transform(Y)
        else:
            Y = Y.values

        normalizer = None

        # Y = Y.values
        # if normalized:
        #     normalizer = MaxAbsScaler()
        #     Y = normalizer.fit_transform(Y)
        #     Y = Y * 100
        # print Y.shape

        if cut_off:
            indices = np.where(Y.sum(axis=0) < cut_off)[0]
            X = np.delete(X, indices, axis=0)
            Y = np.delete(Y, indices, axis=0)

        return X,Y,normalizer
