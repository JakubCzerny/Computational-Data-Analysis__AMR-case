import matplotlib.pyplot as plt
from utils.Data import Data
import utils.plots as plots


# We should look at metadata and discard wrong samples


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
plt.show()
