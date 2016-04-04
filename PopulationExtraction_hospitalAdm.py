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
import cPickle as pickle


# create the patient id column
def createID(subject_id, hospadm_seq):
    f = lambda x, y: str(x) + '#%#' + str(y)
    ids = map(f, subject_id, hospadm_seq)
    return ids


# clean dataset, such as removing NaN and convert to correct data types
def dataClean(data, retrieved_id_name, charttime_name, column_itemid, selected_pt_ids):
    # get the patient id for this dataset
    data['new_id'] = createID(data['subject_id'], data['hospital_seq'])
    if len(selected_pt_ids) > 0:
        data = data[data['new_id'].isin(selected_pt_ids)]
    if len(retrieved_id_name) > 0:
        data = data[data[retrieved_id_name].notnull()]
    if len(charttime_name) > 0:
        # pd.to_datetime(data[charttime_name])
        data = data[data[charttime_name].notnull()]
    # select columns according to itemid and remove NaN
    data = data[data.itemid.notnull()]
    if len(column_itemid) > 0:
        data = data[data.itemid.isin(column_itemid)]
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
    data2 = data2.sort(['new_id', 'charttime'], ascending=[1, 1], inplace=False)
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
                abx_bld_hours = (abx_times2[j] - bld_times2[m]).total_seconds() / 3600.
                if (abx_bld_hours >= -24) and (abx_bld_hours <= 72):
                    select_index.append(i)
                    earlier_time.append(min(abx_times2[j], bld_times2[m]))
                    select = 1
                    break
            if select == 1:
                break
    ptids_selected = np.array(ptids)[select_index]
    pts_abx_bld = pd.DataFrame({'new_id': ptids_selected, 'charttime': earlier_time})
    return pts_abx_bld


def selectPatients2(abx_bld, sofa, ptids):
    '''
    select patients that meet the abx and bld, who also meet the sofa value and timing criteria of new sepsis definition 3
    '''
    select_index = []
    for i in range(len(ptids)):
        pid = ptids[i]
        abx_bld_time = abx_bld[abx_bld['new_id'] == pid]['charttime'].values[0]

        f = lambda x: datetime.strptime(x, '%Y-%m-%d %H:%M:%S')
        sofa_times = sofa[sofa['new_id'] == pid]['charttime'].values
        sofa_times2 = map(f, sofa_times)
        sofa_values = sofa[sofa['new_id'] == pid]['value1num'].values
        prior_sofa_values = []
        post_sofa_values = []
        for m in range(len(sofa_times2)):
            abxbld_sofa_hours = (abx_bld_time - sofa_times2[m]).total_seconds() / 3600.
            if (abxbld_sofa_hours > 0) and (abxbld_sofa_hours <= 24):
                prior_sofa_values.append(sofa_values[m])
            elif (abxbld_sofa_hours <= 0) and (abxbld_sofa_hours >= -24):
                post_sofa_values.append(sofa_values[m])
        if len(prior_sofa_values) == 0:
            prior_sofa_values.append(0)
        if (len(post_sofa_values) > 0) and (max(post_sofa_values) - min(prior_sofa_values) >= 2):
            select_index.append(i)
    ptids_selected = np.array(ptids)[select_index]
    pts_abx_bld_sofa = abx_bld[abx_bld['new_id'].isin(ptids_selected)]
    return pts_abx_bld_sofa


def selectPatients3(abx_bld, sofa, ptids):
    '''
    select patients that meet the abx and bld, who also meet the sofa value and timing criteria of new sepsis definition 3
    '''
    select_index = []
    for i in range(len(ptids)):
        pid = ptids[i]
        abx_bld_time = abx_bld[abx_bld['new_id'] == pid]['charttime'].values[0]

        f = lambda x: datetime.strptime(x, '%Y-%m-%d %H:%M:%S')
        sofa_times = sofa[sofa['new_id'] == pid]['charttime'].values
        sofa_times2 = map(f, sofa_times)
        sofa_values = sofa[sofa['new_id'] == pid]['value1num'].values
        prior_post_sofa_values = []
        for m in range(len(sofa_times2)):
            abxbld_sofa_hours = (abx_bld_time - sofa_times2[m]).total_seconds() / 3600.
            if (abxbld_sofa_hours <= 24) and (abxbld_sofa_hours >= -24):
                prior_post_sofa_values.append(sofa_values[m])
        increase = 0
        if len(prior_post_sofa_values) == 1:
            increase = prior_post_sofa_values[0]
        elif len(prior_post_sofa_values) > 1:
            increase_all = []
            for n in range(len(prior_post_sofa_values)-1):
                increase_all.append(max(prior_post_sofa_values[n+1:]) - prior_post_sofa_values[n])
            increase = max(increase_all)
        if increase >= 2:
            select_index.append(i)
    ptids_selected = np.array(ptids)[select_index]
    pts_abx_bld_sofa = abx_bld[abx_bld['new_id'].isin(ptids_selected)]
    return pts_abx_bld_sofa


def lactateTimes(ptEvents, lactates, ptids):
    select_pt_hl = [] # ptids of those with high lactate > 2
    lactateVevent = []
    lactateTevent = []
    clear_times = []
    normal_times = []
    last_lactate_times = []
    sepsis_time = []
    for i in range(len(ptids)):
        pid = ptids[i]
        eventtime = ptEvents[ptEvents['new_id'] == pid]['charttime'].values[0]
        f = lambda x: datetime.strptime(x, '%Y-%m-%d %H:%M:%S')
        lactate_times = lactates[lactates['new_id'] == pid]['charttime'].values
        lactate_times2 = map(f, lactate_times)
        lactate_values = lactates[lactates['new_id'] == pid]['valuenum'].values
        lactateT = ''
        lactateV = ''
        t = 0
        # get the lactate closest to the time of event
        for m in range(len(lactate_times2)):
            diff_hours = (eventtime - lactate_times2[m]).total_seconds() / 3600.
            if (diff_hours <= 48) and (diff_hours >= 0) and (lactate_values[m] > 2):
                lactateT = lactate_times2[m]
                lactateV = lactate_values[m]
                t = m
        if lactateV != '':
            select_pt_hl.append(pid)
            lactateTevent.append(lactateT)
            lactateVevent.append(lactateV)
            # get the lactate clear time and normalization time
            clearT = ''
            normalT = ''
            for n in range(t+1, len(lactate_times2)):
                if lactate_values[n] <= 0.9 * lactateV:
                    clearT = lactate_times2[n]
                if lactate_values[n] <= 2:
                    normalT = lactate_times2[n]
            clear_times.append(clearT)
            normal_times.append(normalT)
            last_lactate_times.append(lactate_times2[-1])
            sepsis_time.append(eventtime)

    return [select_pt_hl, sepsis_time, lactateVevent, lactateTevent, clear_times, normal_times, last_lactate_times]


def timeDiff(times1, times2):
    f = lambda x, y: (x - y).total_seconds() / 3600. if x != '' else ''
    diffs = map(f, times1, times2)
    return diffs


def fillResponse(time):
    f = lambda x: 1 if x != '' else 0
    response = map(f, time)
    return response


def fillEmpty(time1, time2):
    f = lambda x, y: x if x != '' else y
    response = map(f, time1, time2)
    return response


def extractPredictor(data, events, ptids):
    data_new = []
    for i in range(len(ptids)):
        pid = ptids[i]
        eventtime = events[events['new_id'] == pid]['sepsis_lac_time'].values[0]
        f = lambda x: datetime.strptime(x, '%Y-%m-%d %H:%M:%S')
        times = data[data['new_id'] == pid]['charttime'].values
        times2 = map(f, times)
        variables = []
        data2 = data[data['new_id'] == pid]
        for m in range(len(times2)):
            diff_hours = (eventtime - times2[m]).total_seconds() / 3600.
            if (diff_hours <= 12) and (diff_hours >= -24):
                variables.append(data2.iloc[m])
        data_new += variables
    return data_new


def consistentItemid(labs_x):
    labs_x2 = labs_x
    labs_x2['itemid'] = labs_x2['itemid'].replace([50172], 50025)
    labs_x2['itemid'] = labs_x2['itemid'].replace([50014], 50013)
    labs_x2['itemid'] = labs_x2['itemid'].replace([50112], 50006)
    labs_x2['itemid'] = labs_x2['itemid'].replace([50386], 50007)
    labs_x2['itemid'] = labs_x2['itemid'].replace([50395], 50419)
    labs_x2['itemid'] = labs_x2['itemid'].replace([50460], 50440)
    labs_x2['itemid'] = labs_x2['itemid'].replace([50429], 50428)
    labs_x2['itemid'] = labs_x2['itemid'].replace([50149], 50009)
    labs_x2['itemid'] = labs_x2['itemid'].replace([50159], 50012)
    labs_x2['itemid'] = labs_x2['itemid'].replace([50468], 50316)
    return labs_x2

if __name__ == '__main__':

    #===============Preprocess Blood Culture Data==============================
    # read the blood culture data
    bld = pd.read_csv("pts_bldcultures3_new2.csv")

    # get the icuseq info for bld culture
    # icuseqs = pd.read_csv('icuseqs_blds.csv', header=None)
    # icuseqs.columns = ['icu_seq']
    # bld['icustay_seq'] = icuseqs['icu_seq']
    bld = bld.rename(columns={'spec_itemid': 'itemid'})
    bld2 = dataClean(bld, [], 'charttime',  [70011, 70012], [])

    # export the patient ids to select the antibiotics data from SQL database
    pts_bld_ids_sepsis3 = bld2[['subject_id', 'hospital_seq']]
    pts_bld_ids_sepsis3 = pts_bld_ids_sepsis3.drop_duplicates()
    # pts_bld_ids_sepsis3.to_csv('pts_bld_ids_sepsis3.csv', header=True)

    # extract the useful bld data columns
    pts_bld = bld2[['new_id', 'itemid', 'charttime']]
    # pts_bld.to_csv('pts_bld.csv', header=True)

    # sort bld culture data by id and charttime
    pts_bld = pts_bld.sort(['new_id', 'charttime'], ascending=[1, 1], inplace=False)

    # get the unique blood culture patient new_ids
    pts_bld_id = list(set(pts_bld['new_id'].values))
    pts_bld_id.sort() # 13,812


     #===============Preprocess Antibiotics Data==============================
    # read the abx data
    abx = pd.read_csv("pts_bld_sepsis3_antibiotics_hospitaladm.csv")
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
    pts_abx_bld_ids0.sort() # 9,775
    # select patients according to the pts_abx_bld_ids0
    pts_abx2 = selectSort(pts_abx, pts_abx_bld_ids0)
    pts_bld2 = selectSort(pts_bld, pts_abx_bld_ids0)

    # select patients with abx and bld timinig meets the new sepsis definition
    pts_abx_bld = selectPatients(pts_abx2, pts_bld2, pts_abx_bld_ids0)
    pts_abx_bld_id = list(set(pts_abx_bld['new_id'].values))  #9,117 patients
    #
    with open('pts_abx_bld.pickle', 'wb') as f:
        pickle.dump(pts_abx_bld, f)

    # with open('pts_abx_bld_id.pickle', 'wb') as f:
    #     pickle.dump(pts_abx_bld_id, f)

    #===============Preprocess SOFA Data==============================
    # read the sofa data
    charts = pd.read_csv("pts_bld_sepsis3_vitals_hospitaladm.csv")
    # get the hospital_seq info for sofa data
    # hospseqs = pd.read_csv('hospseqs_vitals.csv', header=None)
    # hospseqs.columns = ['hospital_seq']
    # sofa['hospital_seq'] = hospseqs['hospital_seq']
    sofa2 = dataClean(charts, [], 'charttime',  [20009], [])

    # extract the useful abx data columns
    pts_sofa = sofa2[['new_id', 'itemid', 'charttime', 'value1num']]
    pts_sofa = pts_sofa[pts_sofa['value1num'].notnull()]

    # get the unique sofa patient new_ids
    pts_sofa_id = list(set(pts_sofa['new_id'].values))
    pts_sofa_id.sort()

    # get the intersect patient ids from bld and abx data
    pts_abx_bld_sofa_ids0 = list(set(pts_abx_bld_id).intersection(set(pts_sofa_id)))
    pts_abx_bld_sofa_ids0.sort()  # 8,904

    # select patients according to the pts_abx_bld_sofa_ids0
    # pd.to_datetime(pts_abx_bld['charttime'])
    # pts_abx_bld2 = selectSort(pts_abx_bld, pts_abx_bld_sofa_ids0)
    pts_sofa2 = selectSort(pts_sofa, pts_abx_bld_sofa_ids0)

    # select patients with abx,bld,sofa timinig and value meets the new sepsis definition
    pts_abx_bld_sofa = selectPatients2(pts_abx_bld, pts_sofa2, pts_abx_bld_sofa_ids0)
    pts_abx_bld_sofa_id = set(pts_abx_bld_sofa['new_id'].values) # 2,836
    # pts_abx_bld_sofa.to_csv('pts_abx_bld_sofa.csv', header=True)

    pts_abx_bld_sofa2 = selectPatients3(pts_abx_bld, pts_sofa2, pts_abx_bld_sofa_ids0)
    pts_abx_bld_sofa2_id = set(pts_abx_bld_sofa2['new_id'].values) #2,915

    pts_abx_bld_sofa_id_all = list(pts_abx_bld_sofa_id.union(pts_abx_bld_sofa2_id)) # 2,947 pts

    with open('ptids_sepsis3def.pickle', 'wb') as f:
        pickle.dump(pts_abx_bld_sofa_id_all, f)

    pd.Series(pts_abx_bld_sofa_id_all).to_csv('ptids_sepsis3def.csv', header=True)

#=================================================================================================================

    with open('ptids_sepsis3def.pickle', 'rb') as f:
        pts_abx_bld_sofa_id_all = pickle.load(f)


    # get the patients with lactate values
    # ===lactate values from charts=========
    lactate_charts = dataClean(charts, [], 'charttime',  [1531, 818], pts_abx_bld_sofa_id_all)
    lactate_charts = lactate_charts[['new_id', 'charttime', 'itemid', 'value1num']]
    lactate_charts = lactate_charts.rename(columns={'value1num': 'valuenum'})

    # ===lactate values from lab tests=========
    labs = pd.read_csv("pts_bld_sepsis3_labs_hospitaladm.csv")
    lactate_labs = dataClean(labs, [], 'charttime',  [50010], pts_abx_bld_sofa_id_all)
    lactate_labs = lactate_labs[['new_id', 'charttime', 'itemid', 'valuenum']]

    #combine the lactates from both charts and tests=========
    lactates = pd.concat([lactate_charts, lactate_labs])
    lactates = lactates.sort(['new_id', 'charttime'], ascending=[1, 1], inplace=False)

    sepsis_lactate_ids0 = list(set(lactates['new_id'].values).intersection(set(pts_abx_bld_sofa_id_all)))
    sepsis_lactate_ids0.sort()  # 2,501

    # get the lactate clearance time and
    sepsis_lactate_infos = lactateTimes(pts_abx_bld, lactates, sepsis_lactate_ids0)
        # 795 patients when use -24 ~ 0 hours with at least one lactate > 2
        # 836 patients when use -48 ~ 0 hours with at least one lactate > 2 *

    sepsis_lactate_infos_pd = pd.DataFrame(sepsis_lactate_infos).transpose()
    sepsis_lactate_infos_pd.columns = ['new_id', 'sepsis_time', 'sepsis_lac_value', 'sepsis_lac_time', 'clear_time', 'normalize_time', 'last_lac_time']

    with open('sepsis_lactate_infos.pickle', 'wb') as f:
        pickle.dump(sepsis_lactate_infos_pd, f)

    sepsis_lactate_infos_pd.to_csv('sepsis_lactate_infos.csv', header=True)


    sepsis_lactate_id = sepsis_lactate_infos[0]
    with open('ptids_sepsis3def_lactate2.pickle', 'wb') as f:
        pickle.dump(sepsis_lactate_id, f)

    pd.Series(sepsis_lactate_id).to_csv('ptids_sepsis3def_lactate2.csv', header=True)


    with open('sepsis_lactate_infos.pickle', 'rb') as f:
        sepsis_lactate_infos_pd = pickle.load(f)


    #=====================extract the lab and vitals of the patients who meet the sepsis 3 definition and lactate > 2==============
    # y response: the patient is able to clear lactate in 12 hours from sepsis onset
    sepsis_lactate_infos_pd['clear_duration'] = timeDiff(sepsis_lactate_infos_pd['clear_time'], sepsis_lactate_infos_pd['sepsis_lac_time'])
    sepsis_lactate_infos_pd['normalize_duration'] = timeDiff(sepsis_lactate_infos_pd['normalize_time'], sepsis_lactate_infos_pd['sepsis_lac_time'])
    sepsis_lactate_infos_pd['last_duration'] = timeDiff(sepsis_lactate_infos_pd['last_lac_time'], sepsis_lactate_infos_pd['sepsis_lac_time'])

    sepsis_lactate_infos_pd['cleared?'] = fillResponse(sepsis_lactate_infos_pd['clear_time'])
    sepsis_lactate_infos_pd['normalized?'] = fillResponse(sepsis_lactate_infos_pd['normalize_time'])

    sepsis_lactate_infos_pd['clear_duration'] = fillEmpty(sepsis_lactate_infos_pd['clear_duration'], sepsis_lactate_infos_pd['last_duration'])
    sepsis_lactate_infos_pd['normalize_duration'] = fillEmpty(sepsis_lactate_infos_pd['normalize_duration'], sepsis_lactate_infos_pd['last_duration'])

    sepsis_lactate_infos_pd2 = sepsis_lactate_infos_pd[sepsis_lactate_infos_pd['last_duration'] > 0] # 741 patients
    sepsis_lactate_id2 = sepsis_lactate_infos_pd2['new_id'].values

    with open('sepsis_lactate_infos_times.pickle', 'wb') as f:
        pickle.dump(sepsis_lactate_infos_pd2, f)

    sepsis_lactate_infos_pd2.to_csv('sepsis_lactate_infos_times.csv', header=True)

    # select vital signs for these patients
    charts = pd.read_csv("pts_bld_sepsis3_vitals_hospitaladm.csv")
    charts_x = dataClean(charts, [], 'charttime',  [211, 618, 678, 455, 198, 52], sepsis_lactate_id2)
    # get the diastolic BP
    charts_SDBP = charts_x[charts_x['itemid'] == 455]
    charts_SDBP['itemid'] = 4552
    charts_SDBP['value1num'] = charts_SDBP['value2num']
    charts_x = pd.concat([charts_x, charts_SDBP])
    charts_x = charts_x[['new_id', 'charttime', 'itemid', 'value1num']]
    charts_x2 = charts_x.rename(columns={'value1num': 'valuenum'})

    charts_variables0 = extractPredictor(charts_x2, sepsis_lactate_infos_pd2, sepsis_lactate_id2)
    charts_variables = pd.DataFrame(charts_variables0, columns=['new_id', 'charttime', 'itemid', 'valuenum'])
    charts_ids = set(charts_variables['new_id'].values) # 723 patients

    with open('charts_variables_sepsis3.pickle', 'wb') as f:
        pickle.dump(charts_variables, f)

    # select lab tests for these patients
    labs = pd.read_csv("pts_bld_sepsis3_labs_hospitaladm.csv")
    item_labs = [50014, 50460, 50060, 50061, 50062, 50073, 50002, 50025, 50172, 50177, 50079, 50090, 50013, 50006, 50112, 50386, 50383,
                 50007, 50184, 50140, 50419, 50395, 50015, 50440, 50016, 50018, 50428, 50429, 50019, 50009, 50149, 50439,
                 50399, 50012, 50159, 50170, 50171, 50316, 50122, 50188, 50068, 50091, 50148, 50451, 50316, 50468]
    labs_x = dataClean(labs, [], 'charttime',  item_labs, sepsis_lactate_id2)

    labs_x = labs_x[['new_id', 'charttime', 'itemid', 'valuenum']]
    labs_x2 = consistentItemid(labs_x)

    labs_variables0 = extractPredictor(labs_x, sepsis_lactate_infos_pd2, sepsis_lactate_id2)
    labs_variables = pd.DataFrame(labs_variables0, columns=['new_id', 'charttime', 'itemid', 'valuenum'])
    labs_ids = set(labs_variables['new_id'].values)  # 741 patients

    with open('labs_variables_sepsis3.pickle', 'wb') as f:
        pickle.dump(labs_variables, f)

    charts_labs = pd.concat([charts_variables, labs_variables])
    charts_labs = charts_labs.sort(['new_id', 'itemid', 'charttime'], ascending=[1, 1, 1], inplace=False)
    with open('predictors_sepsis3.pickle', 'wb') as f:
        pickle.dump(charts_labs, f)
    charts_labs_ids = set(charts_labs['new_id'].values)  # 739 patients
    print 'patient numbers: %i' % len(charts_labs_ids)

