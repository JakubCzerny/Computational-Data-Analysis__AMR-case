import numpy as np
import matplotlib.pyplot as plt
import utils.plots as plots
from utils.Data import Data
from sklearn.linear_model import ElasticNet, ElasticNetCV, MultiTaskElasticNetCV
from sklearn.model_selection import train_test_split

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


'''
Binned cummulative number of resistent bacterias per sample
Grouped by country
2 variants - absolute & normalized
'''
plots.binned_cummulative(cummulative_per_country, bins=1)
plots.binned_cummulative(cummulative_per_country, bins=2)
plots.binned_cummulative(cummulative_per_country, bins=4)

# Display all the plots together
# plt.show()



'''
Build a model for predicting the amount of each bacteria based on the used drags
Apply sparse linear regression (ElasticNet) for easier analysis
Force the coefficients to be non-negative as none drug should increase the presence of the bacterias
'''

# Use cut_off to get rid off samples that have less than cut_off amount of bacterias
X,Y = data.get_gene_count_with_drugs(cut_off=30)
print X.shape
print Y.shape
X_train, X_test, Y_train, Y_test = train_test_split(X, Y, test_size=0.25, random_state=1)

folds = 5
alphas = np.logspace(-8,3,10)
l1_ratios = np.linspace(0,1,10,endpoint=True)

# Elastic-net

# NO IDEA WHY IT TELLS ME TO USE MultiTaskElasticNetCV

# models = ElasticNetCV(l1_ratio=l1_ratios, alphas=alphas, verbose=1, cv=folds, positive=True, n_jobs=-1)
models = MultiTaskElasticNetCV(l1_ratio=l1_ratios, alphas=alphas, verbose=1, cv=folds, n_jobs=-1)

models.fit(X_train, Y_train)
model.score(X_test, Y_test)

model_EN = ElasticNet(l1_ratio=models.l1_ratio_, alpha=models.alpha_)
model_EN.fit([X_train,X_test], [Y_train,Y_test])

print "Alpha: ", models.alpha_
print "L1 ratio: ", models.l1_ratio_
print "Score of Elastic-net on test data: ", model.score(X_test, Y_test)
