__author__ = 'jinghe'

'''
This script is used to extract the sepsis patient population according to the Sepsis 3 definition
# working directory
cd /media/sf_Box_Sync/Hyperlactemia sepsis project_jinghe/Data/

'''

# import packages
import pandas as pd
import numpy as np
from datetime import datetime

# create the patient id column
def createID(subject_id, hospadm_seq, icustay_seq):
    f = lambda x, y, z: str(x) + '#%#' + str(y) + '#%#' + str(z)
    ids = map(f, subject_id, hospadm_seq, icustay_seq)
    return ids

# clean dataset, such as removing NaN and convert to correct data types
def dataClean(data, retrieved_id_name, charttime_name, column_itemid, selected_pt_ids):
    if len(retrieved_id_name) > 0:
        data = data[data[retrieved_id_name].notnull()]
    pd.to_datetime(data[charttime_name])
    # select columns according to itemid and remove NaN
    data = data[data.itemid.notnull()]
    if len(column_itemid) > 0:
        data = data[data.itemid.isin(column_itemid)]
    data = data[data[charttime_name].notnull()]
    # get the patient id for this dataset
    data['new_id'] = createID(data['subject_id'], data['hospital_seq'], data['icustay_seq'])
    if len(selected_pt_ids) > 0:
        data = data[data['new_id'].isin(selected_pt_ids)]
    return data


def abxCategory():
    category = {}
    category['Cephalosporins'] = ['CeftazIDIME', 'Ceftazidime', 'CeftAZIDime', 'Ceftriaxone', 'CefTRIAXone', 'CeftriaXONE', 'CefePIME', 'Cefepime']
    category['Penicillins'] = ['Piperacillin-Tazobactam Na', 'Piperacillin Sodium', 'Piperacillin', 'Ampicillin-Sulbactam',
                               'Ampicillin Sodium', 'Ampicillin Sodium/Sulbactam', 'Unasyn', 'Penicillin G Potassium', 'Nafcillin',
                               "Ampicillin Desensitization", 'Penicillin G K Desensitization', 'Ampicillin']
    category['Fluoroquinolones'] = ['Ciprofloxacin', 'Ciprofloxacin IV', 'Levofloxacin']
    category['Aminoglycosides'] = ['Tobramycin Sulfate', 'Gentamicin Sulfate', 'Gentamicin', 'Amikacin']
    category['Carbapenems'] = ['Imipenem-Cilastatin', 'Meropenem']
    category['Macrolide'] = ['Azithromycin ', 'Erythromycin Lactobionate', 'Erythromycin']
    category['Antiviral'] = ['Ribavirin *NF*', 'Acyclovir', 'Acyclovir Sodium']
    category['Antifungal'] = ['Amphotericin B', 'Fluconazole', "Amphotericin B Lipo. (Ambisome)"]
    category['Vancomycin'] = ['Vancomycin', 'Vancomycin HCl', "Vancomycin Enema"]
    category['Others'] = ['Aztreonam', 'Tigecycline', 'Doxycycline Hyclate', 'Colistin', 'Linezolid', 'MetRONIDAZOLE (FLagyl)',
                          'Metronidazole', 'Clindamycin', 'Clindamycin Phosphate', 'Sulfameth/Trimethoprim', 'Daptomycin',
                          'Rifampin', 'Quinupristin/Dalfopristin', 'Pentamidine Isethionate', 'Cefotetan']
    category_reverse = {}
    for key, value in category.iteritems():
        for v in value:
            if not category_reverse.__contains__(v):
                category_reverse[v] = key
    return category_reverse


def GeneralizeAbx(category_reverse, data):
    category_names = []
    for d in data:
        if not category_reverse.__contains__(d):
            category_reverse[d] = 'Others'
        category_names.append(category_reverse[d])
    return category_names, category_reverse


def selectSort(data, selected_ids):
    '''
    select data rows with the selected_ids and then sort the data frame according to the new_id and charttime
    :param data:
    :param selected_ids:
    :return:
    '''
    data2 = data[data['new_id'].isin(selected_ids)]
    data2.sort(['new_id', 'charttime'], ascending=[1, 1], inplace=True)
    return data2


def selectPatients(abx, bld, ptids):
    '''
    select patients that meet the abx and bld criteria of new sepsis definition 3
    :param abx:
    :param bld:
    :param ptids:
    :return:
    '''
    select_index = []
    earlier_time = []
    for i in range(len(ptids)):
        pid = ptids[i]
        abx_times = abx[abx['new_id'] == pid]['charttime'].values
        f = lambda x: datetime.strptime(x, '%Y-%m-%d %H:%M:%S')
        abx_times2 = map(f, abx_times)
        bld_times = bld[bld['new_id'] == pid]['charttime'].values
        bld_times2 = map(f, bld_times)
        select = 0
        for j in range(len(abx_times2)):
            for m in range(len(bld_times2)):
                abx_bld_hours = float((abx_times2[j] - bld_times2[m]).seconds) / 3600
                if (abx_bld_hours >= -24) and (abx_bld_hours <= 72):
                    select_index.append(i)
                    earlier_time.append(min(abx_times2[j], bld_times2[m]))
                    select = 1
                    break
            if select == 1:
                break
    ptids_selected = np.array(ptids)[select_index]
    pts_abx_bld = pd.DataFrame({'new_id': ptids_selected, 'earlier_time': earlier_time})
    return pts_abx_bld


def selectPatients2(abx_bld, sofa, ptids):
    '''
    select patients that meet the abx and bld, who also meet the sofa value and timing criteria of new sepsis definition 3
    '''
    select_index = []
    for i in range(len(ptids)):
        pid = ptids[i]
        abx_bld_time = abx_bld[abx_bld['new_id'] == pid]['earlier_time'].values[0]

        f = lambda x: datetime.strptime(x, '%Y-%m-%d %H:%M:%S')
        sofa_times = sofa[sofa['new_id'] == pid]['charttime'].values
        sofa_times2 = map(f, sofa_times)
        sofa_values = sofa[sofa['new_id'] == pid]['value1num'].values
        prior_sofa_values = []
        post_sofa_values = []
        for m in range(len(sofa_times2)):
            abxbld_sofa_hours = float((abx_bld_time - sofa_times2[m]).seconds) / 3600
            if (abxbld_sofa_hours > 0) and (abxbld_sofa_hours <= 24):
                prior_sofa_values.append(sofa_values[m])
            elif (abxbld_sofa_hours <= 0) and (abxbld_sofa_hours >= -24):
                post_sofa_values.append(sofa_values[m])
        if len(prior_sofa_values) == 0:
            prior_sofa_values.append(0)
        elif (len(post_sofa_values) > 0) and (min(post_sofa_values) - min(prior_sofa_values) >= 2):
            select_index.append(i)
    ptids_selected = np.array(ptids)[select_index]
    pts_abx_bld_sofa = abx_bld[abx_bld['new_id'].isin(ptids_selected)]
    return pts_abx_bld_sofa

if __name__ == '__main__':

    #===============Preprocess Blood Culture Data==============================
    # read the blood culture data
    bld = pd.read_csv("pts_bldcultures3_new.csv")

    # get the icuseq info for bld culture
    icuseqs = pd.read_csv('icuseqs_blds.csv', header=None)
    icuseqs.columns = ['icu_seq']
    bld['icustay_seq'] = icuseqs['icu_seq']
    bld = bld.rename(columns={'spec_itemid': 'itemid'})
    bld2 = dataClean(bld, 'icustay_seq', 'charttime',  [70011, 70012], [])

    # export the patient ids to select the antibiotics data from SQL database
    pts_bld_ids_sepsis3 = bld2[['subject_id', 'hospital_seq', 'icustay_seq']]
    # pts_bld_ids_sepsis3.to_csv('pts_bld_ids_sepsis3.csv', header=True)

    # extract the useful bld data columns
    pts_bld = bld2[['new_id', 'itemid', 'charttime']]
    # pts_bld.to_csv('pts_bld.csv', header=True)

    # sort bld culture data by id and charttime
    pts_bld.sort(['new_id', 'charttime'], ascending=[1, 1], inplace=True)

    # get the unique blood culture patient new_ids
    pts_bld_id = list(set(pts_bld['new_id'].values))
    pts_bld_id.sort()


     #===============Preprocess Antibiotics Data==============================
    # read the abx data
    abx = pd.read_csv("pts_bld_sepsis3_antibiotics.csv")
    abx = abx.rename(columns={'start_dt': 'charttime'})
    abx = abx.rename(columns={'drug': 'itemid'})
    abx2 = dataClean(abx, [], 'charttime', [], pts_bld_id)

    # extract the useful abx data columns
    pts_abx = abx2[['new_id', 'itemid', 'charttime', 'stop_dt']]
    # pts_abx.to_csv('pts_abx.csv', header=True)

    # get the unique abx patient new_ids
    pts_abx_id = list(set(pts_abx['new_id'].values))
    pts_abx_id.sort()

    # get the intersect patient ids from bld and abx data
    pts_abx_bld_ids0 = list(set(pts_abx_id).intersection(set(pts_bld_id)))
    pts_abx_bld_ids0.sort()
    # select patients according to the pts_abx_bld_ids0
    pts_abx2 = selectSort(pts_abx, pts_abx_bld_ids0)
    pts_bld2 = selectSort(pts_bld, pts_abx_bld_ids0)

    # select patients with abx and bld timinig meets the new sepsis definition
    pts_abx_bld = selectPatients(pts_abx2, pts_bld2, pts_abx_bld_ids0)
    pts_abx_bld_id = list(set(pts_abx_bld['new_id'].values))

    #===============Preprocess SOFA Data==============================
    # read the sofa data
    sofa = pd.read_csv("pts_vitals1b.csv")
    # get the hospital_seq info for sofa data
    hospseqs = pd.read_csv('hospseqs_vitals.csv', header=None)
    hospseqs.columns = ['hospital_seq']
    sofa['hospital_seq'] = hospseqs['hospital_seq']
    sofa2 = dataClean(sofa, 'hospital_seq', 'charttime',  [20009], pts_abx_bld_id)

    # extract the useful abx data columns
    pts_sofa = sofa2[['new_id', 'itemid', 'charttime', 'value1num']]
    pts_sofa = pts_sofa[pts_sofa['value1num'].notnull()]

    # get the unique sofa patient new_ids
    pts_sofa_id = list(set(pts_sofa['new_id'].values))
    pts_sofa_id.sort()

    # get the intersect patient ids from bld and abx data
    pts_abx_bld_sofa_ids0 = list(set(pts_abx_bld_id).intersection(set(pts_sofa_id)))
    pts_abx_bld_sofa_ids0.sort()

    # select patients according to the pts_abx_bld_sofa_ids0
    pts_abx2 = selectSort(pts_abx, pts_abx_bld_sofa_ids0)
    pts_bld2 = selectSort(pts_bld, pts_abx_bld_sofa_ids0)
    pts_sofa2 = selectSort(pts_sofa, pts_abx_bld_sofa_ids0)

    # select patients with abx,bld,sofa timinig and value meets the new sepsis definition
    pts_abx_bld_sofa = selectPatients2(pts_abx_bld, pts_sofa2, pts_abx_bld_sofa_ids0)
    pts_abx_bld_sofa_id = list(set(pts_abx_bld_sofa['new_id'].values))
    pts_abx_bld_sofa.to_csv('pts_abx_bld_sofa.csv', header=True)
