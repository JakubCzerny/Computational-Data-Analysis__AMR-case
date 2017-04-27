from utils.Data import Data
from sklearn.cluster import KMeans
from sklearn.cross_validation import train_test_split


data = Data()

#individuals. Transposed for having observations as rows
individuals = data.gene_counts_df.iloc[:,1:].T
kmeans = KMeans(n_clusters=9, n_jobs=3)
train_individuals, test_individuals = train_test_split(individuals, train_size = 0.5)
kmeans.fit(train_individuals)

result = kmeans.predict(test_individuals)

#let's determine which individuals belong to each cluster
#for further exploration
clustered_individuals = {}
for idx,cluster in enumerate(result):
    if not clustered_individuals.has_key(cluster):
        clustered_individuals[cluster] = []
    clustered_individuals[cluster].append(test_individuals[idx].name)

#WE CAN USE SILHOUETTE METRIC FOR INTRA/INTER CLUSTER DISTANCE. CHECKIT IT OUT IN SKLEARN DOCS

# For each cluster, let's answer some questions.

# Predominant country
