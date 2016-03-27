__author__ = 'jinghe'
'''
Analysis on patients with HIV and sepsis, where SOFA score is the max SOFA >= 2
'''


# import packages
import pandas as pd
import numpy as np
from datetime import datetime
import PopulationExtraction as extraction
import cPickle as pickle


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



if __name__ == '__main__':


    # ============================= get all HIV patients that meet the sepsis criteria ==============================
    # read hiv patient identified using ICD-9 codes and create patient ids

    hiv_pts = pd.read_csv('HIVicd9.csv')
    hiv_pts['new_id'] = extraction.createID(hiv_pts['subject_id'], hiv_pts['hospital_seq'], hiv_pts['icustay_seq'])
    hiv_pts_id = hiv_pts['new_id'].values
    # get the max sofa score of patients

    infos = pd.read_csv('pticu_infos.csv')
    infos['new_id'] = extraction.createID(infos['subject_id'], infos['hospital_seq'], infos['icustay_seq'])
    # select only patients with max SOFA >= 2
    infos_sofa2 = infos[infos['sofa_max'] >= 2]
    infos_sofa2_id = infos_sofa2['new_id'].values

    with open('pts_abx_bld.pickle', 'rb') as f:
        pts_abx_bld = pickle.load(f)

    # get the pt ids of both max sofa >=2 and HIV

    pts_sofa_hiv = list(set(infos_sofa2_id).intersection(set(hiv_pts_id)))

    hiv_sepsis_pts = pts_abx_bld[pts_abx_bld['new_id'].isin(pts_sofa_hiv)]
    hiv_sepsis_pts_id = hiv_sepsis_pts['new_id'].values

    # ========================== analyze patient mortalities =====================================================

    # mortality analysis among these patients
    infos_hiv_sepsis = infos[infos['new_id'].isin(hiv_sepsis_pts_id)]
    infos_hiv_sepsis = pd.merge(infos_hiv_sepsis, hiv_sepsis_pts, how='inner', on='new_id')

    mortalities = mortalityComputation(infos_hiv_sepsis)


    # ========================== analyze the prediction performance of MEWS =========================================
    # ============== get the variables values in MEWS ========================================

    '''Respiratory rate; Heart rate; Systolic blood pressure; Conscious level; Temperature; Hourly urine output (for previous 2 hours)'''
