import random
from typing import List
import os

import pandas as pd
import numpy as np
from pandas import DataFrame
from scipy.special import expit

from helper.formula import SynFormula
from helper import populate_df, show_hist

LINEAR = "linear"
QUADRATIC = "quadratic"
SQUARE_ROOT = "sqrt"
EXPONENTIAL = "exponential"


def gen_new_data(num_occ_covariates=8, og_file="weta_data_baseline.csv"):
    """
    generates a synthetic species given a csv of checklists and writes
    the new species to a file.

    after creating occ/det probabilities with the given formulas,
    the probability distributions are shown in a histogram.  The
    user is asked if they want to save these formulas

    :param num_occ_covariates:  the number of occupancy covariates to include.
                                    observation; < 10 results in highly separated
                                    distributions
    :param og_file:             original file containing checklists of a particular
                                species and environment variables appended to each
                                checklist
    :return:
    the new dataframe
    """

    df = pd.read_csv(og_file)
    while True:
        occ_obj = gen_occupation_data(df, num_occ_covariates)
        det_obj = gen_detection_data(df)

        both_form = occ_obj.formula + det_obj.formula

        show_hist.make_histograms(occupied=occ_obj.gen_values, detected=det_obj.gen_values)
        to_save = input("save these formulas? (y/n) ")

        if to_save == 'y':
            i = 0
            while os.path.exists(f"data/syn_spec/syn_species_{i}.csv"):
                i += 1
            file_name = f"data/syn_spec/syn_species_{i}.csv"
            formula_file = f"data/syn_spec/syn_species_{i}_formula.csv"
            save_formula(both_form, formula_file)
            print(f"saving dataframe to {file_name}")
            prob_df = populate_df.populate_df(
                df=df,
                occupied=occ_obj.gen_values,
                detected=det_obj.gen_values,
                f_out=file_name
            )

            final_df_txt = f"data/zero_one_syn_spec/syn_species_{i}.csv"
            populate_df.zero_or_one_df(prob_df, final_df_txt)
            break
        else:
            print("data looks bad :( ... generating new formulas!!")

    return df


def gen_occupation_data(df, num_occ_covariates: int):
    """
    generates a formula for a synthetic species by randomly assigning beta
    coefficients to num_occ_covariates covariates.
    these covariates are from

    :param df: dataframe of checklists
    :param num_occ_covariates: the number of occupany covariates
    :return:
    a SynFormula object with occupation information
    """
    numerics = ['int16', 'int32', 'int64', 'float16', 'float32', 'float64']
    num_df = df.select_dtypes(include=numerics)

    # TODO: good practice: check for highly correlated values
    occ_covariates = random.choices(num_df.columns[30:], k=num_occ_covariates)
    print(f"covariates: {occ_covariates}")
    occ_formula = ""
    occ_coefficients = None
    occ_prob = None
    occupied = None

    rerun = True
    while rerun:
        occ_formula = "OCCUPANCY FORMULA IS: \n\n logit(o_i) = 0"
        occ_formula, occ_coefficients = make_coefficients(occ_formula, occ_covariates)
        print(occ_formula)

        occupied = gen_single_values(df, occ_covariates, occ_coefficients)

        # loose checks to ensure the data is reasonable. if its not, it will generate a new formula
        num1s = 0
        num90s = 0
        num10s = 0
        size = len(occupied)
        for i in occupied:
            if float(i) == 1:
                num1s += 1
            if .9 <= float(i) < 1:
                num90s += 1
            if 0 <= float(i) <= .1:
                num10s += 1
        occ_prob = (num1s / size)

        print(f"num 10s: {num10s}")
        print(f"num 90s: {num90s}")
        print(f"num 1s: {num1s}")

        if num1s < (size * .7):
            # if the number of prob that are between .9 and 1 is greater than 70%, rerun
            if num90s > (size * .7):
                print(f"num90s: {num90s}")
                print("rerunning (prob is bad)")
                rerun = True
            # if the number of prob that are between 0 and .1 is greater than 70%, rerun
            elif num10s > (size * .7):
                print(f"num10s: {num10s}")
                print("rerunning (prob is bad)")
                rerun = True
            else:
                os.system('clear')
                print(occ_formula)
                rerun = False
        else:
            rerun = True
            print("generating new (occ) formula!")

    occupied_object = SynFormula(
        formula=occ_formula,
        covariate_list=occ_covariates,
        coefficient_dict=occ_coefficients,
        new_values=occupied,
        probability=occ_prob
    )

    return occupied_object


def gen_detection_data(df: DataFrame):
    """
    generates a detection formula for a synthetic species by randomly
    assigning beta coefficients to 5 covariates.
    these covariates are:
        1. day_of_year
        2. time_observations_started
        3. duration_minutes
        4. effort_distance_km
        5. number_observers

    :param df:  the dataframe of interest
    :return:
    a SynFormula object with detection information
    """
    det_covariates = [
        'day_of_year', 'time_observations_started', 'duration_minutes', 'effort_distance_km', 'number_observers'
    ]
    det_formula = ""
    detected = None
    det_coefficients = None
    det_prob = None

    rerun = True
    while rerun:
        det_formula = "\nDETECTION FORMULA IS: \n\n logit(o_i) = 0"
        det_formula, det_coefficients = make_coefficients(det_formula, det_covariates)
        detected = gen_single_values(df, det_covariates, det_coefficients)

        num1s = 0
        num10s = 0
        num90s = 0
        size = len(detected)
        for i in detected:
            if int(i) == 1:
                num1s += 1
            if .9 <= float(i) < 1:
                num90s += 1
            if 0 <= float(i) < .1:
                num10s += 1

        det_prob = num1s / size
        print(f"num 10s: {num10s}")
        print(f"num 90s: {num90s}")
        print(f"num 1s: {num1s}")

        if num1s < (size * .7):
            # if the number of prob that are between .9 and 1 is greater than 20%, rerun
            if num90s > (size * .7):
                print(f"num90s: {num90s}")
                print("rerunning (prob is bad)")
                rerun = True
                # if the number of prob that are between 0 and .1 is greater than 70%, rerun
            elif num10s > (size * .7):
                print(f"num10s: {num10s}")
                print("rerunning (prob is bad)")
                rerun = True
            else:
                print(det_formula)
                rerun = False
        else:
            rerun = True
            print("generating new (det) formula!")

    detection_obj = SynFormula(
        formula=det_formula,
        covariate_list=det_covariates,
        coefficient_dict=det_coefficients,
        new_values=detected,
        probability=det_prob
    )

    return detection_obj


def make_coefficients(formula: str, covariate_list: List[str]):
    """
    randomly generates a formula

    :param formula:         string notating what type of formula this is
    :param covariate_list:  list of covariates
    :return:
    the updated formula (string)
    a dictionary of coefficients where
        key: covariate string
        value: (coefficient float, variable transformation string)
    """

    coefficients = {}
    for cov in covariate_list:
        beta_i = np.random.normal(1, .5, 1)[0]
        if random.randint(0, 1):
            beta_i *= -1
        form = random.randint(0, 100)
        # linear
        if form < 25:
            formula += f" + {beta_i} * {cov}"
            coefficients[cov] = (beta_i, LINEAR)
        # sqrt
        elif 25 <= form < 50:
            formula += f" + {beta_i} * {cov} ^(1/2)"
            coefficients[cov] = (beta_i, SQUARE_ROOT)
        # quadratic
        elif 50 <= form < 75:
            formula += f" + {beta_i} * {cov} ^2"
            coefficients[cov] = (beta_i, QUADRATIC)
        # exponential
        else:
            formula += f" + {beta_i} * e^(-{cov})"
            coefficients[cov] = (beta_i, EXPONENTIAL)

    return formula, coefficients


def gen_single_values(df, covariate_list: List[str], coefficient_dict: dict):
    """
    generates probabilities of a single column given a list of covariates
    and the corresponding beta coefficients

    :param df:                  dataframe of interest
    :param covariate_list:      list of covariates as strings
    :param coefficient_dict:    dictionary of coefficients where
                                    key: covariate string
                                    value: (coefficient float, variable transformation string)
    :return:
    list of values between 0 and 1
    """
    col_values = []
    for index, row in df.iterrows():
        col_sum = 0
        for covariate in covariate_list:
            val = float(row[covariate])
            # normalizing step ... logical?
            while val > 100:
                val = val / 100
            coefficient_type = coefficient_dict[covariate][1]
            coefficient = coefficient_dict[covariate][0]
            if coefficient_type == LINEAR:
                col_sum += coefficient * val
            elif coefficient_type == SQUARE_ROOT:
                if val < 0:
                    val = -val
                col_sum += coefficient * np.math.pow(val, .5)
            elif coefficient_type == QUADRATIC:
                col_sum += coefficient * np.math.pow(val, 2)
            elif coefficient_type == EXPONENTIAL:
                col_sum += np.math.exp(-val)
            else:
                print(f"covariate: {covariate} coefficient {coefficient_type} is not linear, quad, or exp ...")
                assert 1 == 0

        col_values.append(col_sum)

    return expit(col_values)


def save_formula(form: str, f_in: str):
    """
    writes a the formula as a string to a file
    :param form:    formula to write
    :param f_in:    file to write formula to
    :return:
    None
    """

    with open(f_in, 'a') as f_save:
        to_write = form + "\n"
        f_save.write(to_write)
        f_save.close()


if __name__ == "__main__":
    gen_new_data()

