import sys
from sklearn import metrics

from pandas import read_csv
from imblearn.over_sampling import SMOTENC, ADASYN


def resample(df):
    """
    calculates a suite of statistics on a dataframe before and after a given sampling technique
    :param df:              dataframe with occupied and species_observed_syn columns
    :param sampling_tech:   the sampling technique
    :return:

    None
    """

    # df = df.head(50)
    print("~~~~BEFORE SAMPLING~~~~")
    print(calc_stats(df))

    # creating an array of indices specifying the categorical features;
    numerics = ['int16', 'int32', 'int64', 'float16', 'float32', 'float64']
    df = df.fillna(value=0)
    num_df = df.select_dtypes(include=numerics)

    numeric_column_list = num_df.columns

    categorical_column_idx = []
    for i in range(len(df.columns)):
        if df.columns[i] not in numeric_column_list:
            categorical_column_idx.append(i)

    # print(numeric_column_list)
    # print(categorical_column_idx)

    smote = SMOTENC(categorical_features=categorical_column_idx)
    sampled_data, axis = smote.fit_resample(df, df.species_observed_syn)

    print("~~~~AFTER SAMPLING~~~~")
    print(calc_stats(sampled_data))


def calc_stats(df):

    occupied = df.occupied
    detected = df.species_observed_syn

    # auc_score = metrics.roc_auc_score(occupied, detected)

    precision, recall, occ_0s, occ_1s, det_0s, det_1s = calc_prec_rec(occupied, detected)

    occ_mean = occ_1s / (occ_0s + occ_1s)
    det_mean = det_1s / (det_0s + det_1s)
    return \
        f"""
        
        number of occupations:      {occ_1s}
        number of non-occupations:  {occ_0s}
        occupation mean:            {occ_mean}
        
        
        number of detections:       {det_1s}
        number of non-detections:   {det_0s}
        detection mean:             {det_mean}          
        
        precision:                  {precision}
        recall:                     {recall}

        """
#         roc_auc_score:              {auc_score}

def calc_prec_rec(occupied, detected):
    false_p = 0
    false_n = 0
    true_p = 0
    true_n = 0

    occ_0s = 0
    occ_1s = 0

    det_0s = 0
    det_1s = 0

    # calc stats
    for idx in range(len(occupied)):
        if occupied[idx] == 1.0:
            occ_1s += 1

            if detected[idx] == 1.0:
                det_1s += 1
                true_p += 1
            elif detected[idx] == 0.0:
                det_0s += 1
                false_n += 1

        elif occupied[idx] == 0.0:
            occ_0s += 1

            if detected[idx] == 1.0:
                det_1s += 1
                false_p += 1
            elif detected[idx] == 0.0:
                det_0s += 1
                true_n += 1

    precision = true_p / (true_p + false_p)
    recall = true_p / (true_p + false_n)
    return precision, recall, occ_0s, occ_1s, det_0s, det_1s


if __name__ == "__main__":
    f_in = sys.argv[1]
    syn_df = read_csv(f_in)

    # smotenc = SMOTENC
    resample(syn_df)
