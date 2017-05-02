from utils.Data import Data
from sklearn.cluster import KMeans
from sklearn.cross_validation import train_test_split
from collections import Counter
from sklearn.metrics import silhouette_score, silhouette_samples
from sklearn.preprocessing import MaxAbsScaler
import numpy as np
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA

data = Data()

# change index and remap to string of metadata
data.metadata.set_index('sample_code', inplace=True)
data.metadata.index = data.metadata.index.map(unicode)
data.gene_counts_df.iloc[1:,1:] = data.gene_counts_df.iloc[1:,1:].astype(float)


# Normalize with respect to number of pairs
for individual in data.metadata.index:
    data.gene_counts_df.ix[1:, individual] = data.gene_counts_df.ix[1:, individual].apply(lambda x: np.divide(x, data.metadata.ix[individual, 'norm_Bacteria_pairs']))

scaler = MaxAbsScaler()
data.gene_counts_df[data.gene_counts_df.columns[1:]] = scaler.fit_transform(data.gene_counts_df[data.gene_counts_df.columns[1:]])
data.gene_counts_df.drop(['900118001','206217001','207618001'], 1, inplace=True)

pigs_df = data.gene_counts_df.T.loc[lambda df: data.metadata.ix[df.index,'type'] == 'Pooled pig feces', :]
poultry_df = data.gene_counts_df.T.loc[lambda df: data.metadata.ix[df.index,'type'] == 'Pooled poultry feces', :]
individuals_df = data.gene_counts_df.T

# Normalized data
pca = PCA(n_components=2)
# pca.fit(individuals_df[1:])
pca.fit(data.gene_counts_df.T[1:])
projection_pigs = pca.transform(pigs_df)
projection_poultry = pca.transform(poultry_df)
plot_pigs = plt.scatter(projection_pigs[:,0], projection_pigs[:,1], marker='x',color='red')
plot_poultry = plt.scatter(projection_poultry[:,0], projection_poultry[:,1],color='green')
plt.legend((plot_pigs, plot_poultry),
           ('Swine', 'Poultry'),
           scatterpoints=1,
           loc='high right')
print 'Variance explained: ', sum(pca.explained_variance_ratio_)


n = 10
avg_pigs_top = pigs_df.mean().sort_values()[-n:]
pigs_df[avg_pigs_top.index].plot(kind='box')
plt.title(str(n) + ' most popular bacteria genes in swine samples')
# plt.plot(range(0,(n+1)),[pigs_df.mean().mean()]*(n+1), '-r')
plt.xticks(range(1,n+1), data.gene_counts_df['clust_name'][avg_pigs_top.index])
print 'Normalizd average value of bacterias for pigs ', pigs_df.mean().mean()


avg_poultry_top = poultry_df.mean().sort_values()[-n:]
poultry_df[avg_poultry_top.index].plot(kind='box')
plt.title(str(n) + ' most popular bacteria genes in poultry samples')
# plt.plot(range(0,(n+1)),[poultry_df.mean().mean()]*(n+1), '-r')
plt.xticks(range(1,n+1), data.gene_counts_df['clust_name'][avg_poultry_top.index])
print 'Normalizd average value of bacterias for poultry ', poultry_df.mean().mean()

plt.show()
