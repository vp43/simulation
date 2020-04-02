import sys
import pandas as pd

import matplotlib.pyplot as plt


def make_histograms(occupied, detected):
    """
    makes histograms of the occupation and detection probabilities

    :param occupied: array of generated occupied values
    :param detected: array of generated detection values
    :return:
    None
    """
    print("exit to continue")
    fig, axs = plt.subplots(1, 2, sharey=True, tight_layout=True)

    axs[0].hist(occupied, bins=10)
    axs[1].hist(detected, bins=10)
    plt.show()


def show_histogram_file(f_in):

    df = pd.read_csv(f_in)

    occupied_prob = df.occupied_prob
    detected_prob = df.species_observed_syn_prob
    make_histograms(occupied_prob, detected_prob)

    occupied = df.occupied
    detected = df.species_observed_syn
    make_histograms(occupied, detected)


if __name__ == "__main__":
    syn_spec = sys.argv[1]
    show_histogram_file(syn_spec)
