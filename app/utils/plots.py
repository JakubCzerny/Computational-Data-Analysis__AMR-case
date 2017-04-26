from __future__ import division
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np

options = ['normalized','absolute']

def cummulative(cummulative_per_country, show=False):
    f, ax = plt.subplots(2)
    ticks_labels = cummulative_per_country['absolute'].keys()
    ticks_positions = [[] for i in range(len(ticks_labels))]

    for subplot, option in enumerate(options):
        width = 0.75
        ind = 0
        gap = 5

        for i, (country, samples) in enumerate(cummulative_per_country[option].iteritems()):
            for j, sample in enumerate(samples):
                ax[subplot].bar(ind+i*gap, sample, width, color=cm.Vega20(i+1))
                ticks_positions[i].append(ind+i*gap)
                ind += 1

        centers = [(ticks_positions[i][0] + ticks_positions[i][-1])/2 for i in range(len(ticks_labels))]
        ax[subplot].set_xticks(centers)
        ax[subplot].set_xticklabels(ticks_labels)

    ax[0].set_title('Normalized amount of resistant bacterias per sample per country')
    ax[1].set_title('Absolute number of resistant bacterias per sample per country')
    plt.tight_layout()

    if show:
        plt.show()


def binned_cummulative(cummulative_per_country, bins=4, show=False):
    binned_cum_per_country = dict((option,{}) for option in options)
    for option in options:
        for country, samples in cummulative_per_country[option].iteritems():
            hist = np.histogram(samples,bins,weights=samples)[0] / np.histogram(samples,bins)[0]
            binned_cum_per_country[option][country] = sorted(hist)

    f, ax = plt.subplots(2)
    ticks_labels = cummulative_per_country['absolute'].keys()
    ticks_positions = [[] for i in range(len(ticks_labels))]

    for subplot, option in enumerate(options):
        width = 0.75
        ind = 0
        gap = 1

        for i, (country, samples) in enumerate(binned_cum_per_country[option].iteritems()):
            for sample in samples:
                ax[subplot].bar(ind+i*gap, sample, width, color=cm.Vega20(i+1))
                ticks_positions[i].append(ind+i*gap)
                ind += 1

        centers = [(ticks_positions[i][0] + ticks_positions[i][-1])/2 for i in range(len(ticks_labels))]
        ax[subplot].set_xticks(centers)
        ax[subplot].set_xticklabels(ticks_labels)

    ax[0].set_title('Normalized binned amount of resistant bacterias per sample per country')
    ax[1].set_title('Absolute binned number of resistant bacterias per sample per country')
    plt.tight_layout()

    if show:
        plt.show()
