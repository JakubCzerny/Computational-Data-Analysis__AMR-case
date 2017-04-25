import csv
import pandas as pd
import matplotlib.pyplot as plt

gene_counts_df = pd.read_csv('../../data/res_finder_gene_count.csv', delimiter='\t')
gene_counts = gene_counts_df.ix[:, gene_counts_df.columns[1:]]

cummulative = gene_counts.sum(axis=0)
metadata = pd.read_csv('../../data/metadata.csv', delimiter='\t')
country_sample = [(metadata.ix[i]['country'], cummulative[i]) for i in metadata.index]

cummulative_per_country = {}
for x, y in country_sample:
    cummulative_per_country.setdefault(x, []).append(y)


fig = plt.figure()
plot = fig.add_subplot(111)
colours = ['black', 'red', 'blue', 'green', 'purple', 'aqua','black', 'red', 'blue', 'green', 'purple', 'aqua']
width = 0.75
ind = 0
gap = 5

for i, (country, samples) in enumerate(cummulative_per_country.iteritems()):
    print samples, country
    colour = colours[i]

    for sample in samples:
        rect = plot.bar(ind+i*gap, sample, width, color=colour)
        ind += 1

plt.show()
