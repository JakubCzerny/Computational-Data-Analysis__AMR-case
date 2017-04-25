import matplotlib.pyplot as plt
import utils.load_data as load_data
import utils.plots as plots


cummulative_per_country = load_data.get_cummulatives()

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
