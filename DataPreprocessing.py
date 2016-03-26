__author__ = 'jinghe'
import datetime
import pandas as pd

def readmultiICUInfoData(filename):
    '''Read the patient icu stay info which contains patient id, hospital admit id, hospital seq, hospital admit time,
    hospital discharge time, icustay id, icustay seq, icustay intime, icustay outtime; read the data into dictionaries by patient'''
    f = open(filename, 'rb')
    pts = {}
    i = 0
    for line in f.xreadlines():
        if i > 0:
            line = line.replace('\r\n', '')
            line = line.split(',')
            if not pts.__contains__(line[0]):
                pts[line[0]] = {}
            if not pts[line[0]].__contains__(line[2]):
                pts[line[0]][line[2]] = {}
                # pts[line[0]][line[2]]['admtime'] = line[3]
                # pts[line[0]][line[2]]['distime'] = line[4]
            if not pts[line[0]][line[2]].__contains__(line[6]):
                pts[line[0]][line[2]][line[6]] = {}
                pts[line[0]][line[2]][line[6]]['intime'] = datetime.datetime.strptime(line[7], '"%Y-%m-%d %H:%M:%S"')
                pts[line[0]][line[2]][line[6]]['outtime'] = datetime.datetime.strptime(line[8], '"%Y-%m-%d %H:%M:%S"')
        i += 1
    print('The number of rows: %s', i)
    return pts


def write(data, filename):
    # convert list to pandas data frame and write to csv
    df = pd.DataFrame(data)
    df.to_csv(str(filename)+'.csv', sep=',', index=False, header=False)


def timeconversion(dt):
    # convert time to half hours
    # dt = datetime.datetime.strptime(t0, '"%Y-%m-%d %H:%M:%S"')
    mins = (dt.minute / 30) * 30
    dt1 = dt.replace(minute=mins, second=0)
    return dt1


def getHALFHOURS(filename):
    # convert time of character type to datetime format
    f = open(filename, 'rb')
    i = 0
    halfhours = []
    for line in f.xreadlines():
        if i > 0:
            line = line.split(',')
            if len(line[1]) > 2 and line[1] != 'NA':
                line[1] = line[1].split('\n')[0]
                t = datetime.datetime.strptime(line[1], '"%Y-%m-%d %H:%M:%S"')
                t_hr = timeconversion(t)
            else:
                t_hr = line[1]
            halfhours.append(t_hr)
        i += 1
    print('The number of rows: %s', i)
    return halfhours


def findICUSeq(filename, infos):
    '''the bldculture data does not have an icu stay id. hence, we use the chartime of the test and compare it with the
     intime and outtime of an icu stay to determine the icu stay id of the observation'''
    f = open(filename, 'rb')
    i = 0
    icuseqall = []
    unmatch = []
    halfhours = []
    for line in f.xreadlines():
        if i > 0:
            line = line.split(',')
            if len(line[3]) > 2:
                t = datetime.datetime.strptime(line[3], '"%Y-%m-%d %H:%M:%S"')
                t_hr = timeconversion(t)
            else:
                t = line[3]
                t_hr = line[3]
            halfhours.append(t_hr)
            if not infos.__contains__(line[0]):
                icuseq = 0
            elif not infos[line[0]].__contains__(line[2]):
                icuseq = 0
            else:
                icuseq = 0
                if len(line[3]) > 2:
                    for key, value in infos[line[0]][line[2]].iteritems():
                        try:
                            if t >= value['intime'] and t <= value['outtime']:
                                icuseq = int(key)
                                break
                        except ValueError:
                            icuseq = 0
                if icuseq == 0 and len(infos[line[0]][line[2]]) == 1:
                    icuseq = 1
            icuseqall.append(icuseq)
            if icuseq == 0:
                unmatch.append(line)
        i += 1
    print('The number of rows: %s', i)
    return icuseqall, unmatch, halfhours


def findICUSeq2(filename, infos):
    '''the bldculture data does not have an icu stay id. hence, we use the chartime of the test and compare it with the
     intime and outtime of an icu stay to determine the icu stay id of the observation; this function expands the time span
     from findICUSeq so that the item taken outside an icu but within 24 hours is also counted as for that icu stay'''
    f = open(filename, 'rb')
    i = 0
    icuseqall = []
    unmatch = []
    halfhours = []
    for line in f.xreadlines():
        if i > 0:
            line = line.split(',')
            if len(line[3]) > 2:
                t = datetime.datetime.strptime(line[3], '"%Y-%m-%d %H:%M:%S"')
                t_hr = timeconversion(t)
            else:
                t = line[3]
                t_hr = line[3]
            halfhours.append(t_hr)
            if not infos.__contains__(line[0]):
                icuseq = 0
            elif not infos[line[0]].__contains__(line[2]):
                icuseq = 0
            else:
                icuseq = 0
                if len(line[3]) > 2:
                    for key, value in infos[line[0]][line[2]].iteritems():
                        try:
                            if t >= value['intime'] and t <= value['outtime']:
                                icuseq = int(key)
                                break
                        except ValueError:
                            icuseq = 0
                if icuseq == 0:
                    if len(infos[line[0]][line[2]]) == 1:
                        icuseq = 1
                    elif len(line[3]) > 2:
                        gaptime = {}
                        for key, value in infos[line[0]][line[2]].iteritems():
                            ingap = max(t, value['intime']) - min(t, value['intime'])
                            outgap = max(t, value['outtime']) - min(t, value['outtime'])
                            gaptime[key] = min(ingap, outgap)
                        mingap = min(gaptime, key=gaptime.get)
                        if gaptime[mingap] <= datetime.timedelta(0, 86400):
                            icuseq = int(mingap)
            icuseqall.append(icuseq)
            if icuseq == 0:
                unmatch.append(line)
        i += 1
    print('The number of rows: %s', i)
    return icuseqall, unmatch, halfhours


pts_hadmicu = readmultiICUInfoData('ptsicustayinfo.csv')
# icuseqs, unmatch, halfhours = findICUSeq('pts_bldcultures3_new.csv', pts_hadmicu)
# halfhours = getHALFHOURS('severe_sepsis_all.csv')
# halfhours_df = pd.DataFrame(halfhours)
#
# write(halfhours, 'halfhours_severe_sepsis_all')
# write(unmatch, 'unmatch_blds')
# write(icuseqs, 'icuseqs_blds')
#
#
# a = 0
# for i in icuseqs:
#     if i == 0:
#         a += 1
# print a #7892
#
# ap = set() #unique patients in the unmatch data
# missingtime = 0 #number of lines with missing time stamp
# ap_missing = set()
# for i in unmatch:
#     ap.add(i[0]) #1023 patients
#     if i[3] == '""':
#         missingtime += 1
#         ap_missing.add(i[0]) #506 patients
# print missingtime #1064
#
# icuseqs_df = pd.DataFrame(icuseqs)
# icuseqs_df.to_csv('icuseqs.csv', sep=',', index=False, header=False)
#
# halfhours_df = pd.DataFrame(halfhours)
# halfhours_df.to_csv('halfhours_blds.csv', sep=',', index=False, header=False)

icuseqs, unmatch, halfhours = findICUSeq2('pts_bldcultures3_new.csv', pts_hadmicu)


a = 0
for i in icuseqs:
    if i == 0:
        a += 1
print a #5646

ap = set()
missingtime = 0
ap_missing = set()
for i in unmatch:
    ap.add(i[0]) #839 patients
    if i[3] == '""':
        missingtime += 1
        ap_missing.add(i[0]) #506 patients
print missingtime #1064


icuseqs_df = pd.DataFrame(icuseqs)
icuseqs_df.to_csv('icuseqs_blds.csv', sep=',', index=False, header=False)

halfhours_df = pd.DataFrame(halfhours)
halfhours_df.to_csv('halfhours_blds.csv', sep=',', index=False, header=False)


icuseqs, unmatch, halfhours = findICUSeq2('pts_catheters.csv', pts_hadmicu)
icuseqs_df = pd.DataFrame(icuseqs)
icuseqs_df.to_csv('icuseqs_catheters.csv', sep=',', index=False, header=False)

halfhours_df = pd.DataFrame(halfhours)
halfhours_df.to_csv('halfhours_catheters.csv', sep=',', index=False, header=False)

def readmultiHospInfoData(filename):
    '''Read the patient icu stay info which contains patient id, hospital admit id, hospital seq, hospital admit time,
    hospital discharge time, icustay id, icustay seq, icustay intime, icustay outtime; read the data into dictionaries by patient'''
    f = open(filename, 'rb')
    pts = {}
    i = 0
    for line in f.xreadlines():
        if i > 0:
            line = line.replace('\r\n', '')
            line = line.split(',')
            if not pts.__contains__(line[0]):
                pts[line[0]] = {}
            if not pts[line[0]].__contains__(line[2]):
                pts[line[0]][line[2]] = set()
            pts[line[0]][line[2]].add(line[5])
                # pts[line[0]][line[2]]['admtime'] = datetime.datetime.strptime(line[3], '"%Y-%m-%d %H:%M:%S"')
                # pts[line[0]][line[2]]['distime'] = datetime.datetime.strptime(line[4], '"%Y-%m-%d %H:%M:%S"')
        i += 1
    print('The number of rows: %s', i)
    return pts


def getcomorb(filename):
    f = open(filename, 'rb')
    comorb = {}
    i = 0
    issues = []
    for line in f.xreadlines():
        if i > 0:
            line = line.split(',')
            try:
                for j in range(len(line)):
                    line[j] = int(line[j])
                if not comorb.__contains__(line[0]):
                    comorb[line[0]] = {}
                if not comorb[line[0]].__contains__(line[2]):
                    comorb[line[0]][line[2]] = sum(line[3:])
            except ValueError:
                issues.append(line)
        i += 1
    print('The number of rows: %s', i)
    return comorb, issues


def findhospSeqComorb(filename, infos, comorbinfo):
    '''get the hospital seq for vitals signs data as well as their comorbidity score'''
    f = open(filename, 'rb')
    i = 0
    hospseqall = []
    comorbscore = []
    for line in f.xreadlines():
        if i > 0:
            line = line.split(',')
            for j in range(0, 3):
                line[j] = int(line[j])
            if not infos.__contains__(line[0]):
                seq = 1
            else:
                seq = 0
                for key, value in infos[line[0]].iteritems():
                    if line[1] in value:
                        seq = int(key)
                        break
            if seq == 0:
                print line
            hospseqall.append(seq)
            # seq = str(seq)
            if comorbinfo.__contains__(line[0]):
                if comorbinfo[line[0]].__contains__(seq):
                    comorbscore.append(comorbinfo[line[0]][seq])
                else:
                    comorbscore.append("NA")
            else:
                comorbscore.append("NA")
        i += 1
    print('The number of rows: %s', i)
    return hospseqall, comorbscore


def comorb(filename, comorbinfo, hseq):
    f = open(filename, 'rb')
    i = 0
    comorbscore = []
    for line in f.xreadlines():
        if i > 0:
            line = line.split(',')
            for j in range(0, hseq + 1):
                line[j] = int(line[j])
            if comorbinfo.__contains__(line[0]):
                if comorbinfo[line[0]].__contains__(line[hseq]):
                    comorbscore.append(comorbinfo[line[0]][line[hseq]])
                else:
                    comorbscore.append("NA")
            else:
                comorbscore.append("NA")
        i += 1
    print('The number of rows: %s', i)
    return comorbscore


def getinfecdysfunc(filename):
    f = open(filename, 'rb')
    icds = {}
    i = 0
    for line in f.xreadlines():
        if i > 0:
            line = line.split(',')
            if len(i[2]) > 0:
                if not icds.__contains__(line[0]):
                    icds[line[0]] = {}
                if not icds[line[0]].__contains__(line[1]):
                    icds[line[0]][line[1]] = 1
        i += 1
    print('The number of rows: %s', i)
    return icds



pts_multihospital = readmultiHospInfoData('multi_hospitaladmit.csv')

comorbinfo, issues = getcomorb('pts_comorbids_new.csv')

hospseqs, comorbs_vitals = findhospSeqComorb('pts_vitals1b_new.csv', pts_multihospital, comorbinfo)

hospseqs_bp, comorbs_vitals_bp = findhospSeqComorb('pts_bldpressures.csv', pts_multihospital, comorbinfo)

bldcomorbs = comorb('pts_bldcultures3_new.csv', comorbinfo, 2)
labcomorbs = comorb('pts_lab_new.csv', comorbinfo, 3)


write(labcomorbs, 'comorbscore_labs')
write(bldcomorbs, 'comorbscore_blds')
write(comorbs_vitals, 'comorbscore_vitals')
write(hospseqs, 'hospseqs_vitals')


halfhours1 = getHALFHOURS("control.group.charttime2.csv")
write(halfhours1, 'halfhours_control_group_charttime')


halfhours2 = getHALFHOURS("control.charttime2.csv")
write(halfhours2, 'halfhours_control_charttime')

halfhours3 = getHALFHOURS("target.charttime2.csv")
write(halfhours3, 'halfhours_target_charttime')

halfhours4 = getHALFHOURS("lab.target.charttime2.csv")
write(halfhours4, 'halfhours_lab_target_charttime')

# deal with the infectious codes
scode = "001, Cholera; 002, Typhoid/paratyphoid fever; 003, Other salmonella infection; 004, Shigellosis; " \
        "005, Other foodpoisoning; 008, Intestinal infection nototherwise classified; 009, Ill-defined intestinalinfection; \
010, Primary tuberculosisinfection; 011, Pulmonary tuberculosis; 012, Other respiratory tuberculosis; \
013, Central nervous system tuberculosis; 014, Intestinal tuberculosis; 015, Tuberculosisof bone and joint; \
016, Genitourinarytuberculosis; 017, Tuberculosisnot otherwise classified; 018, Miliary tuberculosis; 020, Plague; \
021, Tularemia; 022, Anthrax; 023, Brucellosis; 024, Glanders; 025, Melioidosis; 026, Rat-bite fever; 027, Other bacterial zoonoses; \
030, Leprosy; 031, Other mycobacterial disease; \
032, Diphtheria; 033, Whooping cough; \
034, Streptococcal throat/scarlet fever; \
035, Erysipelas; 036, Meningococcal infection; \
037, Tetanus; 038, Septicemia; \
039, Actinomycotic infections; 040, Other\
bacterial diseases; 041, Bacterial infection\
in other diseases not otherwise specified; \
090, Congenital syphilis; 091, Early\
symptomatic syphilis; 092, Early syphilis\
latent; 093, Cardiovascular syphilis; 094,\
Neurosyphilis; 095, Other late symptomatic\
syphilis; 096, Late syphilis latent; \
097, Other and unspecified syphilis; 098,\
Gonococcal infections; 100, Leptospirosis; \
101, Vincent’s angina; 102, Yaws; 103,\
Pinta; 104, Other spirochetal infection; \
110, Dermatophytosis; 111, Dermatomycosis\
not otherwise classified or specified; \
112, Candidiasis; 114, Coccidioidomycosis; \
115, Histoplasmosis; 116, Blastomycotic\
infection; 117, Other mycoses; 118,\
Opportunistic mycoses; 320, Bacterial\
meningitis; 322, Meningitis, unspecified; \
324, Central nervous system abscess; 325,\
Phlebitis of intracranial sinus; 420, Acute\
pericarditis; 421, Acute or subacute endocarditis; \
451, Thrombophlebitis; 461,\
Acute sinusitis; 462, Acute pharyngitis; \
463, Acute tonsillitis; 464, Acute laryngitis/\
tracheitis; 465, Acute upper respiratory\
infection of multiple sites/not otherwise\
specified; 481, Pneumococcal\
pneumonia; 482, Other bacterial pneumonia; \
485, Bronchopneumonia with organism\
not otherwise specified; 486,\
Pneumonia, organism not otherwise\
specified; 491.21, Acute exacerbation of\
obstructive chronic bronchitis; 494,\
Bronchiectasis; 510, Empyema; 513,\
Lung/mediastinum abscess; 540, Acute\
appendicitis; 541, Appendicitis not otherwise\
specified; 542, Other appendicitis; \
562.01, Diverticulitis of small intestine\
without hemorrhage; 562.03, Diverticulitis\
of small intestine with hemorrhage; \
562.11, Diverticulitis of colon without\
hemorrhage; 562.13, Diverticulitis of colon\
with hemorrhage; 566, Anal and rectal\
abscess; 567, Peritonitis; 569.5, Intestinal\
abscess; 569.83, Perforation of\
intestine; 572.0, Abscess of liver; 572.1,\
Portal pyemia; 575.0, Acute cholecystitis; \
590, Kidney infection; 597, Urethritis/\
urethral syndrome; 599.0, Urinary tract\
infection not otherwise specified; 601,\
Prostatic inflammation; 614, Female pelvic\
inflammation disease; 615, Uterine inflammatory\
disease; 616, Other female\
genital inflammation; 681, Cellulitis, finger/\
toe; 682, Other cellulitis or abscess; \
683, Acute lymphadenitis; 686, Other local\
skin infection; 711.0, Pyogenic arthritis; \
730, Osteomyelitis; 790.7, Bacteremia; \
996.6, Infection or inflammation of\
device/graft; 998.5, Postoperative infection; \
999.3, Infectious complication of\
medical care not otherwise classified"

s1 = scode.split('; ')
s0 = []
for i in s1:
    s0.append(i.split(',')[0])
saproxm = []
for i in s0:
    i0 = i.replace('.', '')
    if len(i0) < 5:
        i = i + '%'
    saproxm.append(i)

print saproxm


