from utils.Data import Data
from sklearn.cluster import KMeans
from sklearn.cross_validation import train_test_split
from collections import Counter
from sklearn.metrics import silhouette_score, silhouette_samples
from sklearn.preprocessing import MaxAbsScaler
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm



data = Data()

# change index and remap to string of metadata
data.metadata.set_index('sample_code', inplace=True)
data.metadata.index = data.metadata.index.map(unicode)
data.gene_counts_df.iloc[1:,1:] = data.gene_counts_df.iloc[1:,1:].astype(float)

# data.gene_counts_df.loc[:,data.metadata.index].astype(float, copy=False)

# Normalize with respect to number of pairs
for individual in data.metadata.index:
    data.gene_counts_df.ix[1:, individual] = data.gene_counts_df.ix[1:, individual].apply(lambda x: np.divide(x, data.metadata.ix[individual, 'norm_Bacteria_pairs']))

scaler = MaxAbsScaler(copy=False)
scaler.fit_transform(data.gene_counts_df.iloc[1:, 1:])

pigs_df = data.gene_counts_df.T.loc[lambda df: data.metadata.ix[df.index,'type'] == 'Pooled pig feces', :]
poultry_df = data.gene_counts_df.T.loc[lambda df: data.metadata.ix[df.index,'type'] == 'Pooled poultry feces', :]
individuals_df = data.gene_counts_df.T

c_individuals_total, c_individuals_pigs, c_individuals_poultry = {},{},{}
c_info_total, c_info_pigs, c_info_poultry = {},{},{}
n_clusters = range(3,13)
final_results = []
cluster_info, clustered_individuals = {},{}

for individuals, category in zip([individuals_df, pigs_df, poultry_df], ['total individuals', 'pigs', 'poultry']):

    for n_cluster in n_clusters:

        fig, ax1 = plt.subplots(1)
        ax1.set_xlim([-0.2, 1])
        ax1.set_ylim([0, len(individuals) + (len(n_clusters)+1)*10])

        kmeans = KMeans(n_clusters=n_cluster, n_jobs=3, n_init=50)
        # train_individuals, test_individuals = train_test_split(individuals.iloc[1:,:], train_size = 0.8)
        # kmeans.fit(train_individuals)
        #
        # result = kmeans.predict(test_individuals)
        # No train/test split for unsupervised?
        kmeans.fit(individuals.iloc[1:,:].values)
        cluster_labels = kmeans.predict(individuals.iloc[1:,:].values)

        silhouette_average = silhouette_score(individuals.iloc[1:,:].values, cluster_labels)
        silhouette_values = silhouette_samples(individuals.iloc[1:,:].values, cluster_labels)
        print 'Average score for ' + category +  ' with n_clusters ' + str(n_cluster) + ': ', silhouette_average

        # Prepare plots
        y_lower = 10
        for i in range(len(n_clusters)):
            ith_cluster_silhouette_values = silhouette_values[cluster_labels == i]
            ith_cluster_silhouette_values.sort()
            size_cluster_i = ith_cluster_silhouette_values.shape[0]
            y_upper = y_lower + size_cluster_i

            color = cm.spectral(float(i) / len(n_clusters))
            ax1.fill_betweenx(np.arange(y_lower, y_upper), 0, ith_cluster_silhouette_values,
                                facecolor=color, edgecolor=color, alpha=0.7)

            # Plot labels
            ax1.text(-0.05, y_lower + 0.5 * size_cluster_i, str(i))
            y_lower = y_upper + 10


        #let's determine which individuals belong to each cluster
        #for further exploration
        # clustered_individuals = {}
        for idx,cluster in enumerate(cluster_labels):
            if not clustered_individuals.has_key(cluster):
                clustered_individuals[cluster] = []

            # Test individuals
            # print individuals.index[idx]
            if individuals.index[idx] == 'clust_name':
                clustered_individuals[cluster].append(individuals.index[1:][idx])
            else:
                clustered_individuals[cluster].append(individuals.index[idx])


        #WE CAN USE SILHOUETTE METRIC FOR INTRA/INTER CLUSTER DISTANCE. CHECKIT IT OUT IN SKLEARN DOCS
        # http://scikit-learn.org/stable/auto_examples/cluster/plot_kmeans_silhouette_analysis.html#sphx-glr-auto-examples-cluster-plot-kmeans-silhouette-analysis-py

        # For each cluster, let's answer some questions.

        # cluster_info = {}
        # Predominant country
        # print clustered_individuals
        for cluster, cluster_members in clustered_individuals.iteritems():
            if not cluster_info.has_key(cluster):
                cluster_info[cluster] = {}
            country_list = []
            for member in cluster_members:
                # print member
                country_list.append(data.metadata.ix[member, 'country'])
            c = Counter(country_list)
            cluster_info[cluster]['predominant_country'] = c.most_common(2)
            cluster_info[cluster]['num_members'] = len(cluster_members)
        final_results.append((individuals, clustered_individuals, cluster_info, silhouette_average, silhouette_values, kmeans))
        clustered_individuals, cluster_info = {}, {}

        # Plotting
        ax1.set_title("Silhouette plot for clusters")
        ax1.set_xlabel("Silhouette coefficient values")
        ax1.set_ylabel("Cluster label")

        ax1.axvline(x=silhouette_average, color="red", linestyle='--')

        ax1.set_yticks([])
        ax1.set_xticks(np.arange(-0.5, 1, 0.1))
        plt.suptitle(("Silhouette analysis for KMeans on " + category + " with num_clusters = ", n_cluster),
                        fontsize=12, fontweight='bold')


        # 2nd plot showing actual clusters formed
        # Can't make it unless is 2D PCAed
        # colors = cm.spectral(cluster_labels.astype(float) / len(n_clusters))
        # ax2.scatter(individuals[:,0], X[:,1])

plt.show()
