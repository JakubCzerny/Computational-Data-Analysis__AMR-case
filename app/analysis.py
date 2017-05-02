import numpy as np
import matplotlib.pyplot as plt
import utils.plots as plots
from utils.Data import Data
from sklearn.linear_model import ElasticNet, ElasticNetCV, MultiTaskElasticNetCV
from sklearn.model_selection import train_test_split
import matplotlib.colors as mcolors
import pandas as pd
import seaborn as sns

#TODO We should look at metadata and discard wrong samples

#IDEA:
'''

- extract the baterias that are present in most of the samples and try to do analysis on them
- second option would be to select bacterias which amount varies a lot between the samples
- high in absolute value (negative or possitive) correlation between presence of drug and presence of the resistent bacteria
- cluster the samples based on the presence of bacterias (and their amount), we could use number of cluster equal to the #num_drugs*2 or *3 (should be still readable) => use colors for countries
- linear regression (input drugs and their doses => amount of each bacteria)

'''

# Load the data files
data = Data()

cummulative_per_country = data.get_cummulatives()

'''
Plots of cummulative number of resistent bacterias per sample
Grouped by country
2 variants - absolute & normalized
'''
plots.cummulative(cummulative_per_country)


'''te
Binned cummulative number of resistent bacterias per sample
Grouped by country
2 variants - absolute & normalized
'''
plots.binned_cummulative(cummulative_per_country, bins=1)
plots.binned_cummulative(cummulative_per_country, bins=2)
plots.binned_cummulative(cummulative_per_country, bins=4)

# Display all the plots together
plt.show()

# f = plt.figure()

X,Y,normalizer = data.get_gene_count_with_drugs(cut_off=0)

corrmat = pd.DataFrame(Y).corr()
max_corr = corrmat[abs(corrmat)>0.5]
sns.heatmap(max_corr, vmax=1., square=False).xaxis.tick_top()
plt.show()

'''
Build a model for predicting the amount of each bacteria based on the used drags
Apply sparse linear regression (ElasticNet) for easier analysis
Force the coefficients to be non-negative as none drug should increase the presence of the bacterias
'''

X_train, X_test, Y_train, Y_test = train_test_split(X, Y, test_size=0.25, random_state=1)

folds = 5
alphas = np.logspace(1,5,3)
l1_ratios = np.linspace(0,1,2,endpoint=True)

models = MultiTaskElasticNetCV(l1_ratio=l1_ratios, alphas=alphas, verbose=1, cv=folds, n_jobs=-1)
models.fit(X_train, Y_train)
models.score(X_test, Y_test)

print "Alpha: ", models.alpha_
print "L1 ratio: ", models.l1_ratio_
print "Score of Elastic-net on test data: ", models.score(X_test, Y_test)

model_EN = ElasticNet(l1_ratio=models.l1_ratio_, alpha=models.alpha_)
model_EN.fit(np.concatenate((X_train,X_test)), np.concatenate((Y_train,Y_test)))

test = np.rint(models.predict(X_test)).astype('int16')
coeff = model_EN.coef_.T
# coeff = models.coef_.T

# high=1.0
# low=0.0
# mins = np.min(coeff, axis=0)
# maxs = np.max(coeff, axis=0)
# rng = maxs - mins
# table = (high - (((high - low) * (maxs - coeff)) / rng))

coeff_differentiated = coeff.copy()
neg_ind = np.where(coeff<0)
pos_ind = np.where(coeff>=0)
for i in range(neg_ind[0].shape[0]):
    coeff_differentiated[neg_ind[0][i],neg_ind[1][i]] = coeff[neg_ind[0][i],neg_ind[1][i]] * 2
for i in range(pos_ind[0].shape[0]):
    coeff_differentiated[pos_ind[0][i],pos_ind[1][i]] = coeff[pos_ind[0][i],pos_ind[1][i]] * 2

fig, ax = plt.subplots(1)
# x = ax[0].imshow(table, cmap=plt.cm.seismic)
# plt.colorbar(mappable=x, ax=ax[0])
# ax[0].grid(False)

x = ax.imshow(coeff_differentiated, cmap=plt.cm.seismic)
plt.colorbar(mappable=x, ax=ax)
ax.grid(False)
plt.show()
