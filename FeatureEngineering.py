__author__ = 'jinghe'

'''
This script is used to process the predictors for prediction on lactate clearance/normalization
cd /media/sf_Box_Sync/Hyperlactemia sepsis project_jinghe/Data/

'''

# import packages
import pandas as pd
import numpy as np
from datetime import datetime
import cPickle as pickle
import PopulationExtraction_hospitalAdm as populations
import matplotlib.pyplot as plt
from scipy.stats import mode


def lastsValues(d, n=3):
    f = lambda d: d[-1]
    return f


def renameCols(cols, s):
    f = lambda x: str(x) + '-' + str(s)
    cols2 = map(f, cols)
    return cols2


def draw_histograms(df, variables, n_rows, n_cols, names):
    fig = plt.figure(figsize=(18, 15), dpi=100)
    for i, var_name in enumerate(variables):
        ax = fig.add_subplot(n_rows, n_cols, i+1)
        df[var_name].hist(bins=20, ax=ax, figsize=(15, 15))
        ax.set_title(names[var_name], fontsize=14)
        ax.set_xticklabels(map(int, ax.get_xticks()), fontsize=12)
        ax.set_yticklabels(map(int, ax.get_yticks()), fontsize=12)
    fig.tight_layout()  # Improves appearance a bit.
    plt.show()


def removeCols(df, colnames, threshold):
    cols_2 = [x for x in colnames if df[x].count() > threshold]
    return cols_2


def naImputation(df, colnames):
    for c in colnames:
        imputevalue = df[c].median()
        df[c].fillna(imputevalue, inplace=True)
    return df


def removeRows(data, thres):
    dmiss = data.count(axis=1)
    f = lambda x, y: 1 if x > y else 0
    data2 = map(f, dmiss, [thres]*len(dmiss))
    return data2


def variableNames():
    # codes = [u'52.0-1', u'198.0-1', u'211.0-1', u'455.0-1', u'618.0-1', u'678.0-1', u'4552.0-1', u'50002.0-1', u'50009.0-1',
    #          u'50016.0-1', u'50018.0-1', u'50019.0-1', u'50025.0-1', u'50061.0-1', u'50062.0-1', u'50068.0-1', u'50073.0-1',
    #          u'50079.0-1', u'50090.0-1', u'50112.0-1', u'50140.0-1', u'50148.0-1', u'50149.0-1', u'50159.0-1', u'50170.0-1',
    #          u'50172.0-1', u'50177.0-1', u'50383.0-1', u'50386.0-1', u'50399.0-1', u'50428.0-1', u'50439.0-1', u'50440.0-1',
    #          u'50468.0-1']
    # names = ['MAP', 'GCS', 'HR', 'SystolicBP', 'RR', 'Temp', 'DiastolicBP', 'BaseExcess', 'Potassium', 'PCO2', 'pH', 'PO2',
    #          'Bicarb', 'ALP', 'ALT (GPT)', 'AnionGap', 'AST (GOT)', 'Calcium', 'Creatinine', 'Glucose', 'Magnesium', 'Phosphorus',
    #          'Sodium', '']
    codes = [u'52.0-1', u'198.0-1', u'211.0-1', u'455.0-1', u'618.0-1', u'678.0-1', u'4552.0-1', u'52.0-2', u'198.0-2',
             u'211.0-2', u'455.0-2', u'618.0-2', u'678.0-2', u'4552.0-2', u'52.0-3', u'198.0-3', u'211.0-3', u'455.0-3',
             u'618.0-3', u'678.0-3', u'4552.0-3', u'new_id', u'Response']
    names = ['MAP1', 'GCS1', 'HR1', 'SystolicBP1', 'RR1', 'Temp1', 'DiastolicBP1',
             'MAP2', 'GCS2', 'HR2', 'SystolicBP2', 'RR2', 'Temp2', 'DiastolicBP2',
             'MAP3', 'GCS3', 'HR3', 'SystolicBP3', 'RR3', 'Temp3', 'DiastolicBP3', 'new_id', 'Response']

    code2name_dict = dict.fromkeys(codes)
    for i in range(len(codes)):
        code2name_dict[codes[i]] = names[i]
    return code2name_dict



def UseNames(names_dict, cols):
    cols2 = []
    for c in cols:
        cols2.append(names_dict[c])
    return cols2



if __name__ == '__main__':


    # convert the predictor table to a data matrix: rows are patients and columns are variables
    with open('predictors_sepsis3.pickle', 'rb') as f:
        charts_labs = pickle.load(f)

    with open('charts_variables_sepsis3.pickle', 'rb') as f:
        charts_variables = pickle.load(f)

    # in general, vitals are recorded every one hour; labs are every 1 hours
    # groupby to create a subset of data for each patient


    # charts_labs = charts_variables
    charts_labs_table_1 = charts_labs.pivot_table(values='valuenum', index='new_id', columns='itemid', aggfunc=lambda d: d[-1:])
    charts_labs_table_2 = charts_labs.pivot_table(values='valuenum', index='new_id', columns='itemid', aggfunc=lambda d: d[-2:-1])
    charts_labs_table_3 = charts_labs.pivot_table(values='valuenum', index='new_id', columns='itemid', aggfunc=lambda d: d[-3:-2])

    charts_labs_table_1.columns = renameCols(charts_labs_table_1.columns, 1)
    charts_labs_table_2.columns = renameCols(charts_labs_table_2.columns, 2)
    charts_labs_table_3.columns = renameCols(charts_labs_table_3.columns, 3)

    charts_labs_table = pd.concat([charts_labs_table_1, charts_labs_table_2], axis=1)
    charts_labs_table = pd.concat([charts_labs_table, charts_labs_table_3], axis=1)

    charts_labs_table = charts_labs_table_1
    # remove the columns with more than half missing values
    cols_new = removeCols(charts_labs_table, charts_labs_table.columns, np.round(len(charts_labs_table)/2))
    charts_labs_table_b = charts_labs_table[cols_new] # 21 columns left
    # remove the variables with at least half columns are missing
    charts_labs_table_b['<halfMissing?'] = removeRows(charts_labs_table_b, 15)
    charts_labs_table_c = charts_labs_table_b[charts_labs_table_b['<halfMissing?'] == 1]
    del charts_labs_table_c['<halfMissing?']






    names = variableNames()

    draw_histograms(charts_labs_table_c, charts_labs_table_c.columns, 3, 3, names)

    #============missing value imputation=========================================
    charts_labs_table_d = naImputation(charts_labs_table_c, charts_labs_table_c.columns)
    draw_histograms(charts_labs_table_d, charts_labs_table_d.columns, 3, 3, names)
    # draw_histograms(charts_labs_table_c[charts_labs_table_c.columns[:12]], charts_labs_table_c.columns[:12], 3, 4, names)
    # draw_histograms(charts_labs_table_c[charts_labs_table_c.columns[12:24]], charts_labs_table_c.columns[12:24], 3, 4, names)
    # draw_histograms(charts_labs_table_c[charts_labs_table_c.columns[24:]], charts_labs_table_c.columns[24:], 3, 4, names)

    # transformation on the predictors



    # #prepare the sepsis patients' response variables
    #
    with open('sepsis_lactate_infos.pickle', 'rb') as f:
        sepsis_lactate_infos_pd = pickle.load(f)

    # logistic regression within 24 hours for whether patients can
    response_1 = sepsis_lactate_infos_pd[(sepsis_lactate_infos_pd['cleared?'] == 1) & (sepsis_lactate_infos_pd['clear_duration'] <= 24)]['new_id'].values
    response_0 = sepsis_lactate_infos_pd[((sepsis_lactate_infos_pd['cleared?'] == 1) & (sepsis_lactate_infos_pd['clear_duration'] > 24)) | ((sepsis_lactate_infos_pd['cleared?'] == 0) & (sepsis_lactate_infos_pd['last_duration'] > 24))]['new_id'].values

    charts_labs_table_d['new_id'] = charts_labs_table_d.index
    Xs_1 = charts_labs_table_d[charts_labs_table_d['new_id'].isin(response_1)]
    Xs_1['Response'] = 1
    Xs_0 = charts_labs_table_d[charts_labs_table_d['new_id'].isin(response_0)]
    Xs_0['Response'] = 0

    matrix_data = pd.concat([Xs_1, Xs_0])
    matrix_data.columns = UseNames(names, matrix_data.columns)
    with open('matrix_data_sepsis3_clear24_vitals.pickle', 'wb') as f:
        pickle.dump(matrix_data, f)

    matrix_data.to_csv('matrix_data_sepsis3_clear24_vitals.csv', header=True, index=False)