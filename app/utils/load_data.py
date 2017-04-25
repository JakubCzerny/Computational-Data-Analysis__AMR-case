from __future__ import division
import csv
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np

gene_counts_df = pd.read_csv('../../data/res_finder_gene_count.csv', delimiter='\t')
gene_counts = gene_counts_df.ix[:, gene_counts_df.columns[1:]]

cummulative = gene_counts.sum(axis=0)
metadata = pd.read_csv('../../data/metadata.csv', delimiter='\t')
metadata = pd.read_csv('../../data/metadata.csv', delimiter='\t')

country_sample = {}
country_sample['absolute'] = [(metadata.ix[i]['country'], cummulative.loc[str(metadata.ix[i]['sample_code'])] ) for i in metadata.index]
country_sample['normalized'] = [(metadata.ix[i]['country'], cummulative.loc[str(metadata.ix[i]['sample_code'])] / metadata.ix[i]['norm_Bacteria_pairs'] ) for i in metadata.index]

options = ['normalized','absolute']
cummulative_per_country = dict((option,{}) for option in options)

for option in options:
    for x, y in country_sample[option]:
        cummulative_per_country[option].setdefault(x, []).append(y)


# ======================================= PLOTING =======================================
# f, ax = plt.subplots(2, sharex=True)
#
# for subplot, option in enumerate(options):
#     width = 0.75
#     ind = 0
#     gap = 5
#
#     for i, (country, samples) in enumerate(cummulative_per_country[option].iteritems()):
#         for sample in samples:
#             ax[subplot].bar(ind+i*gap, sample, width, color=cm.Vega20(i+1))
#             ind += 1
#
# ax[0].set_title('Normalized amount of resistant bacterias per sample per country')
# ax[1].set_title('Absolute number of resistant bacterias per sample per country')
# plt.show()
# ========================================================================================


# number of bins per country
no_bins = 2
binned_cum_per_country = dict((option,{}) for option in options)
for option in options:
    for country, samples in cummulative_per_country[option].iteritems():
        bins,edges = np.histogram(samples,no_bins)
        binned_cum_per_country[option][country] = sorted(bins)


# ==================================== PLOTING BINS =======================================
f, ax = plt.subplots(2, sharex=True)

for subplot, option in enumerate(options):
    width = 0.75
    ind = 0
    gap = 1

    for i, (country, samples) in enumerate(binned_cum_per_country[option].iteritems()):
        for sample in samples:
            ax[subplot].bar(ind+i*gap, sample, width, color=cm.Vega20(i+1))
            ind += 1

ax[0].set_title('Normalized amount of resistant bacterias per sample per country')
ax[1].set_title('Absolute number of resistant bacterias per sample per country')
plt.show()
# ========================================================================================
