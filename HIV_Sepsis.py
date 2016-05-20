__author__ = 'jinghe'
'''
Analysis on patients with HIV and sepsis, where SOFA score is the max SOFA >= 2
'''


# import packages
import pandas as pd
import numpy as np
from datetime import datetime
import PopulationExtraction_hospitalAdm as extraction
import cPickle as pickle
import FeatureEngineering as Features
from sklearn import preprocessing
from sklearn import cross_validation, svm, linear_model, metrics
import copy

def Convert2Time(data, colname):
    f = lambda x: datetime.strptime(x, '%Y-%m-%d %H:%M:%S')
    times = map(f, data[colname])
    return times


def timeDiff(times1, times2):
    f = lambda x, y: (x - y).total_seconds() / 3600.
    diffs = map(f, times1, times2)
    return diffs



def mortalityComputation(data):
    result = {}
    # ICU mortality
    icu_mortality = len(data[data['icustay_expire_flg'] == 'Y'])
    icu_mortality_rate = icu_mortality / float(len(data))
    result['icu_mortality'] = (icu_mortality, icu_mortality_rate)

    # hospital mortality
    f1 = lambda x, y: str(x) + '#%#' + str(y)
    hospital_ids = set(map(f1, data['subject_id'], data['hospital_seq']))
    hospital_mortality = len(data[data['hospital_expire_flg'] == 'Y'])
    hospital_mortality_rate = icu_mortality / float(len(hospital_ids))
    result['hospital_mortality'] = (hospital_mortality, hospital_mortality_rate)

    # 30-day mortality from time of event
        # remove the patients still alive
    data2 = data[data['dod'].notnull()]
        # convert the dod into datetime
    f2 = lambda x: datetime.strptime(x, '%Y-%m-%d %H:%M:%S')
    dod_times = map(f2, data2['dod'])
    event_times = data2['charttime'].values
        # compute the dod days from time of event
    f3 = lambda x, y: (x - y).total_seconds()
    mortality_30d = map(f3, dod_times, event_times)
        # patients dod in 30 days
    day30_mortality = len([x for x in mortality_30d if x >= 30 * 24 * 3600])
    subject_ids = set(data['subject_id'].values)
    day30_mortality_rate = day30_mortality / float(len(subject_ids))
    result['30day_mortality'] = (day30_mortality, day30_mortality_rate)

    return result


def variableExtract(data, events, ptids):
    data_new = []
    for i in range(len(ptids)):
        pid = ptids[i]
        eventtime = events[events['new_id'] == pid]['charttime'].values[0]
        f = lambda x: datetime.strptime(x, '%Y-%m-%d %H:%M:%S')
        times = data[data['new_id'] == pid]['charttime'].values
        times2 = map(f, times)
        variables = []
        data2 = data[data['new_id'] == pid]
        for m in range(len(times2)):
            diff_hours = (eventtime - times2[m]).total_seconds() / 3600.
            # if (diff_hours <= 12) and (diff_hours >= -12):
            if diff_hours >= -12:
                variables.append(data2.iloc[m])
        data_new += variables
    return data_new

def computeMEWS_SBP(x):
    if x <= 70:
        return 3
    elif x<= 100 or x >= 200:
        return 2
    elif x <= 199:
        return 0


def computeMEWS_HR(x):
    if x <= 40:
        return 2
    elif x<= 50:
        return 1
    elif x <= 100:
        return 0
    elif x <= 110:
        return 1
    elif x <= 129:
        return 2
    elif x >= 130:
        return 3

def computeMEWS_RR(x):
    if x < 9:
        return 2
    elif x<= 14:
        return 0
    elif x <= 20:
        return 1
    elif x <= 29:
        return 2
    elif x >= 30:
        return 3

def computeMEWS_Temp(x):
    if x < 95:
        return 2
    elif x<= 101.1:
        return 0
    elif x > 101.1:
        return 2


def computeMEWS_AVPU(x):
    if x > 14:
        return 0
    elif x > 9:
        return 1
    elif x > 6:
        return 2
    else:
        return 3


def computeMEWS(data):
    AVPU = map(computeMEWS_AVPU, data['GCS'])
    HR = map(computeMEWS_HR, data['HR'])
    SBP = map(computeMEWS_SBP, data['SBP'])
    RR = map(computeMEWS_RR, data['RR'])
    Temp = map(computeMEWS_Temp, data['Temp'])
    f_sum = lambda a, b, c, d, e: np.sum([a, b, c, d, e])
    score = map(f_sum, AVPU, HR, SBP, RR, Temp)
    return score


def computeNEWS_SBP(x):
    if x < 90:
        return 3
    elif x<= 100 and x >= 91:
        return 2
    elif x <= 110 and x >= 101:
        return 1
    elif x >= 219 and x <= 111:
        return 0
    else:
        return 3


def computeNEWS_HR(x):
    if x <= 40:
        return 3
    elif x <= 50:
        return 1
    elif x <= 90:
        return 0
    elif x <= 110:
        return 1
    elif x <= 130:
        return 2
    else:
        return 3

def computeNEWS_RR(x):
    if x <= 8:
        return 3
    elif x <= 11:
        return 1
    elif x <= 20:
        return 0
    elif x <= 24:
        return 2
    else:
        return 3

def computeNEWS_Temp(x):
    if x <= 95:
        return 3
    elif x <= 96.8:
        return 1
    elif x <= 100.4:
        return 0
    elif x <= 102.2:
        return 1
    else:
        return 2


def computeNEWS_AVPU(x):
    if x >= 15:
        return 0
    else:
        return 3

def computeNEWS_Oxygen(x):
    if x <= 91:
        return 3
    elif x >= 92 and x <= 93:
        return 2
    elif x >= 94 and x <= 95:
        return 1
    else:
        return 0

def computeNEWS(data):
    AVPU = map(computeNEWS_AVPU, data['GCS'])
    HR = map(computeNEWS_HR, data['HR'])
    SBP = map(computeNEWS_SBP, data['SBP'])
    RR = map(computeNEWS_RR, data['RR'])
    Temp = map(computeNEWS_Temp, data['Temp'])
    Oxyg = map(computeNEWS_Oxygen, data['Oxyg'])
    f_sum = lambda a, b, c, d, e, g: np.sum([a, b, c, d, e, g])
    score = map(f_sum, AVPU, HR, SBP, RR, Temp, Oxyg)
    return score



def string2bin(data):
    f = lambda x: 1 if x == 'Y' else 0
    data2 = map(f, data)
    return data2


def otherMetrics(total, pos, tpr, fpr):
    def measures(tpr, fpr, total=total, pos=pos):
        neg = total - pos
        tp = round(tpr * pos)
        fp = fpr * neg
        tn = neg - fp
        tnr = tn / float(neg)
        precision = 0
        if (tp + fp > 0):
            precision = tp / float(tp + fp)
        acc = float(tp + tn) / total
        f1 = 2 * precision * tpr / float(precision + tpr)
        f2 = 5 * precision * tpr / float(4 * precision + tpr)
        return [tnr, precision, acc, f1, f2]
    res = map(measures, tpr, fpr)
    return res


def temp_c2f(df):
    df.loc[df.itemid == 676, 'valuenum'] = df.loc[df.itemid == 676, 'valuenum'] * 9. / 5 + 32
    return df

def consistentItemid_2(labs_x):
    labs_x2 = labs_x
    labs_x2['itemid'] = labs_x2['itemid'].replace([1332], 211)
    labs_x2['itemid'] = labs_x2['itemid'].replace([1341], 211)
    labs_x2['itemid'] = labs_x2['itemid'].replace([3603], 618)
    labs_x2['itemid'] = labs_x2['itemid'].replace([8113], 618)
    labs_x2['itemid'] = labs_x2['itemid'].replace([818], 50010)
    labs_x2['itemid'] = labs_x2['itemid'].replace([1531], 50010)
    labs_x2['itemid'] = labs_x2['itemid'].replace([677], 676)
    labs_x2 = temp_c2f(labs_x2)
    labs_x2['itemid'] = labs_x2['itemid'].replace([679], 678)
    labs_x2['itemid'] = labs_x2['itemid'].replace([676], 678)
    return labs_x2


def imputeMedian(df, names):
    for n in names:
        df[n].fillna(df[n].median(), inplace=True)
    return df


if __name__ == '__main__':

    # ============================= get all HIV patients that meet the sepsis criteria ==============================
    # read hiv patient identified using ICD-9 codes and create patient ids

    hiv_pts = pd.read_csv('pts_hiv_icd9.csv')
    hiv_pts['new_id'] = extraction.createID(hiv_pts['subject_id'], hiv_pts['hospital_seq'])
    hiv_pts_id = hiv_pts['new_id'].values # 464 patients

    # patients with sepsis
    with open('ptids_sepsis3def.pickle', 'rb') as f:
        pts_abx_bld_sofa_id_all = pickle.load(f) # 2,947 patients

    pts_hiv_sepsis_id = set(hiv_pts_id).intersection(pts_abx_bld_sofa_id_all) # 52 patients

    # get the max sofa score of patients

    infos = pd.read_csv('pticu_infosb.csv')
    infos['new_id'] = extraction.createID(infos['subject_id'], infos['hospital_seq'])

    infos['hospital_disch_dt'] = Convert2Time(infos, 'hospital_disch_dt')
    infos['hospital_admit_dt'] = Convert2Time(infos, 'hospital_admit_dt')

    # infos['icustay_intime'] = Convert2Time(infos, 'icustay_intime')
    # infos['icustay_outtime'] = Convert2Time(infos, 'icustay_outtime')

    # infos['icustay_los'] = timeDiff(infos['icustay_outtime'], infos['icustay_intime'])
    # infos['hospital_los'] = timeDiff(infos['hospital_disch_dt'], infos['hospital_admit_dt'])
    # (pd.Series(infos['icustay_los'])).plot(kind='hist')

    # # select only patients with max SOFA >= 2
    infos_sofa2 = infos[infos['sofa_max'] >= 2]
    infos_sofa2_id = infos_sofa2['new_id'].values

    # get the pt ids of both max sofa >=2 and HIV

    pts_sofa_hiv_id = list(set(infos_sofa2_id).intersection(set(hiv_pts_id))) # 386 patients


    # select vital signs for these patients
    charts = pd.read_csv("pts_bld_sepsis3_vitals_hospitaladm.csv")
    charts_x = extraction.dataClean(charts, [], 'charttime',  [211, 618, 678, 455, 198], pts_sofa_hiv_id)
    # charts_x = extraction.dataClean(charts, [], 'charttime',  [211, 1332, 1341, 618, 3603, 8113, 676, 677, 678, 679, 455,818,1531,198], pts_sofa_hiv_id)

    charts_x = charts_x[['new_id', 'charttime', 'itemid', 'value1num']]
    charts_x2 = charts_x.rename(columns={'value1num': 'valuenum'})

    #=======================add CD4 counts======================================
    # add oxygen saturation
    labs = pd.read_csv("pts_bld_sepsis3_labs_hospitaladm.csv")
    labs_x = extraction.dataClean(labs, [], 'charttime', [50015], pts_sofa_hiv_id)
    o2 = labs_x[['new_id', 'charttime', 'itemid', 'valuenum']]

    #combine the lactates from both charts and tests=========
    Xs = pd.concat([charts_x2, o2])
    # Xs_2 = consistentItemid_2(Xs)


    # ===============select patients ==================================
    with open('pts_abx_bld.pickle', 'rb') as f:
        pts_abx_bld = pickle.load(f)
    pts_abx_bld_id = list(set(pts_abx_bld['new_id'].values))
    charts_sepsis_hiv_pts_ids = list(set(pts_abx_bld_id).intersection(pts_sofa_hiv_id))

    infos_hiv_sepsis_hospital_info = infos_sofa2[infos_sofa2['new_id'].isin(charts_sepsis_hiv_pts_ids)]
    infos_hiv_sepsis_hospital_info = infos_hiv_sepsis_hospital_info[['new_id', 'hospital_admit_dt', 'hospital_expire_flg']]
    sepsis_hiv_infos = pd.merge(infos_hiv_sepsis_hospital_info, pts_abx_bld, on='new_id', how='left')
    sepsis_hiv_infos['adm_sepsis_time'] = extraction.timeDiff(sepsis_hiv_infos['charttime'], sepsis_hiv_infos['hospital_admit_dt'])
    sepsis_hiv_infos['adm_sepsis_time'].describe()


    # charts_variables0 = variableExtract(Xs_2, pts_abx_bld, charts_sepsis_hiv_pts_ids)
    charts_variables0 = variableExtract(Xs, pts_abx_bld, charts_sepsis_hiv_pts_ids)
    charts_variables = pd.DataFrame(charts_variables0, columns=['new_id', 'charttime', 'itemid', 'valuenum'])
    charts_ids = set(charts_variables['new_id'].values)# 189 patients $ 143 patients

    charts_table = charts_variables.pivot_table(values='valuenum', index='new_id', columns='itemid', aggfunc=lambda d: d[-1:])
    charts_table.columns = ['GCS', 'HR', 'SBP', 'RR', 'Temp', 'Oxyg']# charts_table.columns = ['GCS', 'HR', 'SBP', 'RR', 'Temp', 'Lactate', 'CD4']
    charts_table['new_id'] = charts_table.index
    infos_sofa2_hospital_expire = infos_sofa2[['new_id', 'hospital_expire_flg']]
    # add abs CD4 counts
    cd4 = pd.read_csv("pts_bld_sepsis3_CD4_hospitaladm.csv")
    cd4_x = extraction.dataClean(cd4, [], [],  [50318], charts_ids)

    charts_table_cd = pd.merge(charts_table, cd4_x, on='new_id', how='left')
    charts_table_cd = pd.merge(charts_table_cd, infos_sofa2_hospital_expire, on='new_id', how='left')

    charts_table_cd = charts_table_cd[['new_id', 'GCS', 'HR', 'SBP', 'RR', 'Temp', 'Oxyg', 'valuenum', 'hospital_expire_flg']]
    charts_table_cd = charts_table_cd.rename(columns={'valuenum': 'CD4'})
    charts_table2 = imputeMedian(charts_table_cd, ['GCS', 'HR', 'SBP', 'RR', 'Temp', 'Oxyg', 'CD4'])

    charts_table2['MEWS_score'] = computeMEWS(charts_table2)
    charts_table2['NEWS_score'] = computeNEWS(charts_table2)

    charts_table3 = charts_table2.drop_duplicates()

    del charts_table3['new_id']
    charts_table3.to_csv('charts_sepsis_hiv_cd4.csv', header=True, index=False)


    #================================================
    charts_table3 = pd.read_csv('charts_sepsis_hiv_cd4_nomogram.csv')

    pos = charts_table3[charts_table3['hospital_expire_flg'] == 'Y']
    neg = charts_table3[charts_table3['hospital_expire_flg'] == 'N']

    # predictors = np.array(charts_table4[['GCS', 'HR', 'SBP', 'RR', 'Temp']])
    # predictors = np.array(charts_table3[['MEWS_score']])
    # predictors2 = np.array(charts_table3[['NEWS_score']])
    # predictors3 = np.array(charts_table3[['MEWS0_CD4_score0']])
    targets = np.array(charts_table3['hospital_expire_flg'].tolist())

    # model = linear_model.LogisticRegression(penalty='l1')
    # cv_scores = cross_validation.cross_val_predict(model, predictors, targets, cv=5)
    # metrics.accuracy_score(targets, cv_scores)
    # metrics.confusion_matrix(targets, cv_scores)
    # lr_report = metrics.classification_report(targets, cv_scores)
    # print(lr_report)

    targets_num = string2bin(targets)
    Mews0 = np.array(charts_table3['MEWS_score'].tolist())
    Mews = np.array(charts_table3['MEWS0_CD4_score0'].tolist())
    # cv_scores_num = string2bin(cv_scores)
    fpr, tpr, thresholds = metrics.roc_curve(targets_num, Mews0, pos_label=1)
    metrics.auc(fpr, tpr)

    res = otherMetrics(144, 22, tpr, fpr)
    MEWS_results = pd.DataFrame(np.array([thresholds, tpr, fpr]).T, columns=['Threshold', 'TPR', 'FPR'])
    res_pd = pd.DataFrame(np.array(res), columns=['tnr', 'precision', 'accuracy', 'F1', 'F2'])
    MEWS_results = pd.concat([MEWS_results, res_pd], axis=1)
    MEWS_results.to_csv('NEWS_prediction_pos_mortality2.csv', header=True)
    MEWS_results = pd.read_csv('NEWS_prediction_pos_mortality.csv')

