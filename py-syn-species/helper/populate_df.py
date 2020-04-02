import os
import sys
from copy import copy

from pandas import read_csv
from numpy import random


def populate_df(df, occupied, detected, f_out):
    """
    appends occupation and detection probabilities of the synthetic species to each existing checklist
    :param df:          the dataframe of interest
    :param occupied:    a list of calculated occupied probabilities. length(occupied) = length(df)
    :param detected:    a list of calculated detected probabilities. length(detected) = length(df)
    :param f_out:       the file name to save the newly populated dataframe
    :return:
    the populated dataframe
    """

    df = df.assign(occupied_prob=occupied)
    df = df.assign(species_observed_syn_prob=detected)
    print(df.head(5))
    df.to_csv(f_out)
    return df


def zero_or_one_df(df, f_out):
    """
    this function appends two columns to the inputted dataframe: occupied and species_observed_syn

    these columns are generated by calculating a 0 or 1 value for occupation and detection
    given the detection and occupation probabilities for every checklist in the dataframe

    :param df:      dataframe containing occupied and detection probabilities stored in
                    the following columns: occupied_prob, species_observed_syn_prob
    :param f_out:   the file name to save the newly populated dataframe
    :return:
    the newly updated
    """

    occupied = copy(df.occupied_prob)
    detected = copy(df.species_observed_syn_prob)

    for idx in range(len(occupied)):
        occupied[idx] = random.binomial(1, occupied[idx], size=None)
        detected[idx] = random.binomial(1, detected[idx], size=None)

    df = df.assign(occupied=occupied)
    df = df.assign(species_observed_syn=detected)
    print(df.head(5))
    df.to_csv(f_out)
    return df


if __name__ == "__main__":
    syn_spec_file_name = sys.argv[1]
    df = read_csv(syn_spec_file_name)
    i = 0
    while os.path.exists(f"data/zero_one_syn_spec/syn_species_{i}.csv"):
        i += 1
    file_name = f"data/zero_one_syn_spec/syn_species_{i}.csv"
    zero_or_one_df(df, f_out=file_name)