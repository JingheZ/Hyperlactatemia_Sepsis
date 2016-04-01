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
            if (diff_hours <= 24) and (diff_hours >= 0):
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

    # select the vitals of these patients

    # select vital signs for these patients
    charts = pd.read_csv("pts_bld_sepsis3_vitals_hospitaladm.csv")
    charts_x = extraction.dataClean(charts, [], 'charttime',  [211, 618, 678, 455, 198], pts_sofa_hiv_id)
    # # get the diastolic BP
    # charts_SDBP = charts_x[charts_x['itemid'] == 455]
    # charts_SDBP['itemid'] = 4552
    # charts_SDBP['value1num'] = charts_SDBP['value2num']
    # charts_x = pd.concat([charts_x, charts_SDBP])
    charts_x = charts_x[['new_id', 'charttime', 'itemid', 'value1num']]
    charts_x2 = charts_x.rename(columns={'value1num': 'valuenum'})

    with open('pts_abx_bld.pickle', 'rb') as f:
        pts_abx_bld = pickle.load(f)
    pts_abx_bld_id = list(set(pts_abx_bld['new_id'].values))
    charts_sepsis_hiv_pts_ids = list(set(pts_abx_bld_id).intersection(pts_sofa_hiv_id))

    charts_variables0 = variableExtract(charts_x2, pts_abx_bld, charts_sepsis_hiv_pts_ids) #214
    charts_variables = pd.DataFrame(charts_variables0, columns=['new_id', 'charttime', 'itemid', 'valuenum'])
    charts_ids = set(charts_variables['new_id'].values) #  57 patients

    charts_table = charts_variables.pivot_table(values='valuenum', index='new_id', columns='itemid', aggfunc=lambda d: d[-1:])

    imputer = preprocessing.Imputer(missing_values='NaN', strategy='median')
    charts_table2 = imputer.fit_transform(charts_table)
    charts_table2 = pd.DataFrame(charts_table2, columns=['GCS', 'HR', 'SBP', 'RR', 'Temp'])
    charts_table2['new_id'] = charts_table.index

    charts_table2['MEWS_score'] = computeMEWS(charts_table2)

    infos_sofa2_hospital_expire = infos_sofa2[['new_id', 'hospital_expire_flg']]
    charts_table3 = pd.merge(charts_table2, infos_sofa2_hospital_expire, on='new_id', how='left')
    charts_table4 = charts_table3.drop_duplicates()
    del charts_table4['new_id']
    charts_table4.to_csv('charts_sepsis_hiv.csv', header=True, index=False)

    pos = charts_table4[charts_table4['hospital_expire_flg'] == 'Y']
    neg = charts_table4[charts_table4['hospital_expire_flg'] == 'N']
    # # ========================== analyze patient mortalities ===================

    # ==================================
    #
    # # mortality analysis among these patients
    # infos_hiv_sepsis = infos[infos['new_id'].isin(hiv_sepsis_pts_id)]
    # infos_hiv_sepsis = pd.merge(infos_hiv_sepsis, hiv_sepsis_pts, how='inner', on='new_id')
    #
    # mortalities = mortalityComputation(infos_hiv_sepsis)
    #
    #
    # # ========================== analyze the prediction performance of MEWS =========================================
    # # ============== get the variables values in MEWS ========================================
    #
    # '''Respiratory rate (618); Heart rate (211); Systolic blood pressure (455 value1num ); Conscious level(GCS 198); Temperature (678); Hourly urine output (for previous 2 hours)'''
    # charts = pd.read_csv("pts_vitals1b.csv")
    # # get the hospital_seq info for sofa data
    # hospseqs = pd.read_csv('hospseqs_vitals.csv', header=None)
    # hospseqs.columns = ['hospital_seq']
    # charts['hospital_seq'] = hospseqs['hospital_seq']
    # items = [618, 678, 198, 455, 211]
    # charts2 = extraction.dataClean(charts, 'hospital_seq', 'charttime',  items, hiv_sepsis_pts_id)
    #
    # charts3 = charts2[['new_id', 'charttime', 'itemid', 'value1num']]
    # charts3 = charts3.sort(['new_id', 'charttime'], ascending=[1, 1], inplace=False)
    # ptids = list(set(hiv_sepsis_pts_id).intersection(set(charts3['new_id'].values)))
    # ptids.sort()
    #
    # # get the IV fluids of patients
    # variables = variableExtract(charts3, infos_hiv_sepsis, ptids)
    # x = infos_hiv_sepsis[['icustay_intime', 'icustay_outtime', 'charttime']]
