setwd('/media/sf_Box_Sync/Hyperlactemia sepsis project_jinghe/Data/')
# setwd('C:/Users/jz4kg/Box Sync/Hyperlactemia sepsis project_jinghe/Data/')

# select patients with max SOFA >= 2

infos <- read.csv('pticu_infos.csv', header = T)
infos$id <- paste(infos$subject_id, infos$hospital_seq, infos$icustay_seq, sep='#%#')
str(infos)
# pts_sofa_all_ids <- infos[(infos$sofa_max >= 2), c(1,2,3,6,7)]
# # write.csv(pts_sofa_all_ids, file = 'pts_sofa_all_ids.csv')
# pts_sofa <- infos$id[(infos$sofa_max >= 2)]

# get the blood culture and abx info for the patients above

# blood culture
bld<- read.csv("pts_bldcultures3_new.csv", header=TRUE)
#get the icuseq info for bld culture
icuseqs <- read.csv('icuseqs_blds.csv',header = F)
bld$icuseq <- as.vector(unlist(icuseqs), mode = 'numeric')
bld <- bld[(bld$icuseq > 0) & (!is.na(bld$spec_itemid)) & (!is.na(bld$charttime)), ]
bld$id <- paste(bld$subject_id, bld$hospital_seq, bld$icuseq, sep='#%#')

pts_bld <- bld
pts_bld <- pts_bld[(!is.na(pts_bld$charttime)) & (pts_bld$charttime != ''),]
pts_bld_ids <- data.frame(pts_bld$subject_id, pts_bld$hospital_seq, pts_bld$icuseq)

write.csv(pts_bld_ids, file = 'pts_bld_ids_sepsis3.csv')

pts_bld <- data.frame(pts_bld$id, pts_bld$spec_itemid, pts_bld$charttime)
names(pts_bld) <- c('id', 'itemid', 'charttime')
write.csv(pts_bld, file = 'pts_bld.csv')

# abx

#get the antibiotics data of those patients
antibiotics <- read.csv('pts_bld_sepsis3_antibiotics.csv')

antibiotics$id <- paste(antibiotics$subject_id, antibiotics$hospital_seq, antibiotics$icustay_seq, sep='#%#')
antibiotics <- antibiotics[antibiotics$id %in% pts_bld_sofa$id,]
antibiotics2 <- data.frame(antibiotics$id, antibiotics$start_dt, antibiotics$stop_dt, antibiotics$drug)
names(antibiotics2) <- c('id', 'start_time', 'stop_time', 'itemid')
antibiotics2 <- antibiotics2[(antibiotics2$itemid != 'Meropenem Desensitization') & (antibiotics2$itemid != 'Hepatitis B Immune Globulin (Nabi-HB)'), ]

#create a dictionary to group antibiotics into more general clusters
# Additional abx in this new population
# [1] "Amphotericin B Lipo. (Ambisome)" "Vancomycin Enema"               
# [3] "Ampicillin Desensitization"  
Cephalosporins <- c('CeftazIDIME', 'Ceftazidime', 'CeftAZIDime', 'Ceftriaxone', 'CefTRIAXone', 'CeftriaXONE', 'CefePIME', 'Cefepime')

Penicillins <- c('Piperacillin-Tazobactam Na', 'Piperacillin Sodium', 'Piperacillin', 'Ampicillin-Sulbactam', 'Ampicillin Sodium', 
                 'Ampicillin Sodium/Sulbactam', 'Unasyn', 'Penicillin G Potassium', 'Nafcillin', "Ampicillin Desensitization", 'Penicillin G K Desensitization')

Fluoroquinolones <- c('Ciprofloxacin', 'Ciprofloxacin IV', 'Levofloxacin')

Aminoglycosides <- c('Tobramycin Sulfate', 'Gentamicin Sulfate', 'Gentamicin', 'Amikacin')

Carbapenems<- c('Imipenem-Cilastatin', 'Meropenem')

Macrolide <- c('Azithromycin ', 'Erythromycin Lactobionate', 'Erythromycin')

Antiviral <- c('Ribavirin *NF*', 'Acyclovir', 'Acyclovir Sodium')

Antifungal <- c('Amphotericin B', 'Fluconazole', "Amphotericin B Lipo. (Ambisome)")

Others <- c('Aztreonam', 'Tigecycline', 'Doxycycline Hyclate', 'Colistin', 'Linezolid', 'MetRONIDAZOLE (FLagyl)', 'Metronidazole', 'Clindamycin',
            'Clindamycin Phosphate', 'Sulfameth/Trimethoprim', 'Daptomycin', 'Rifampin', 'Quinupristin/Dalfopristin', 'Pentamidine Isethionate', 'Cefotetan')

Vancomycin <- c('Vancomycin', 'Vancomycin HCl', "Vancomycin Enema")

matrix_Cephalosporins <- data.frame(Cephalosporins, rep('Cephalosporins', length(Cephalosporins)))
names(matrix_Cephalosporins) <- c('itemid', 'itemidG')

matrix_Penicillins <- data.frame(Penicillins, rep('Penicillins', length(Penicillins)))
names(matrix_Penicillins) <- c('itemid', 'itemidG')

matrix_Fluoroquinolones <- data.frame(Fluoroquinolones, rep('Fluoroquinolones', length(Fluoroquinolones)))
names(matrix_Fluoroquinolones) <- c('itemid', 'itemidG')

matrix_Aminoglycosides <- data.frame(Aminoglycosides, rep('Aminoglycosides', length(Aminoglycosides)))
names(matrix_Aminoglycosides) <- c('itemid', 'itemidG')

matrix_Carbapenems <- data.frame(Carbapenems, rep('Carbapenems', length(Carbapenems)))
names(matrix_Carbapenems) <- c('itemid', 'itemidG')

matrix_Macrolide <- data.frame(Macrolide, rep('Macrolide', length(Macrolide)))
names(matrix_Macrolide) <- c('itemid', 'itemidG')

matrix_Antiviral <- data.frame(Antiviral, rep('Antiviral', length(Antiviral)))
names(matrix_Antiviral) <- c('itemid', 'itemidG')

matrix_Antifungal <- data.frame(Antifungal, rep('Antifungal', length(Antifungal)))
names(matrix_Antifungal) <- c('itemid', 'itemidG')

matrix_Others <- data.frame(Others, rep('Others', length(Others)))
names(matrix_Others) <- c('itemid', 'itemidG')

matrix_Vancomycin <- data.frame(Vancomycin, rep('Vancomycin', length(Vancomycin)))
names(matrix_Vancomycin) <- c('itemid', 'itemidG')

# names(antibiotics3)

antibiotics.group <- rbind(matrix_Cephalosporins, matrix_Penicillins)
antibiotics.group <- rbind(antibiotics.group, matrix_Fluoroquinolones)
antibiotics.group <- rbind(antibiotics.group, matrix_Aminoglycosides)
antibiotics.group <- rbind(antibiotics.group, matrix_Carbapenems)
antibiotics.group <- rbind(antibiotics.group, matrix_Macrolide)
antibiotics.group <- rbind(antibiotics.group, matrix_Antiviral)
antibiotics.group <- rbind(antibiotics.group, matrix_Antifungal)
antibiotics.group <- rbind(antibiotics.group, matrix_Others)
antibiotics.group <- rbind(antibiotics.group, matrix_Vancomycin)

antibiotics2$itemid <- as.character(antibiotics2$itemid)
str(antibiotics.group)
antibiotics.group$itemid <- as.character(antibiotics.group$itemid)

antibiotics3 <- merge(antibiotics2, antibiotics.group, by.x = 'itemid', by.y = 'itemid', all.x = T)
names(antibiotics3)

antibiotics4 <- data.frame(antibiotics3$id, antibiotics3$itemidG, antibiotics3$start_time, antibiotics3$stop_time)
names(antibiotics4) <- c('id', 'itemid', 'charttime', 'stop_time')



# pts_bld_sofa$stop_time <- pts_bld_sofa$charttime
# 
# names(pts_bld_sofa) <-  c( 'charttime', 'itemid', 'id', 'stop_time')
# 
# 
# sofa_bld_abx_pts <- rbind(pts_bld_sofa, antibiotics4)
# sofa_bld_abx_pts <- sofa_bld_abx_pts[(!is.na(sofa_bld_abx_pts$charttime)) & (sofa_bld_abx_pts$charttime != ''), ]
# sofa_bld_abx_pts <- sofa_bld_abx_pts[order(sofa_bld_abx_pts$id, sofa_bld_abx_pts$charttime),]
# write.csv(sofa_bld_abx_pts, file = 'sofa_bld_abx_pts.csv')

# select patients with blood culture taken and abx -- 
# (If the antibiotic was given first, the culture sampling must have been obtained within 24 hours. 
# If the culture sampling was first, the antibiotic must have been ordered within 72 hours)

library(plyr)
library(lubridate)
abx <- antibiotics4 
abx$id <- as.character(abx$id)
abx$charttime <- ymd_hms(abx$charttime)
str(abx)
bld <- pts_bld_sofa
bld$id <- as.character(bld$id)
bld$charttime <- ymd_hms(bld$charttime)
str(bld)
ptids <- as.character(unique(abx$id))

select <- rep(0, length(ptids))
bld_abx_early_time <- rep(NA, length(ptids))
for (i in 1:length(ptids)) {
  pid <- ptids[i]
  abx_times <- abx$charttime[abx$id == pid]
  bld_times <- bld$charttime[bld$id == pid]
  for (j in 1:length(abx_times)) {
    for (m in 1:length(bld_times)) {
      if (((difftime(abx_times[j], bld_times[m], units = 'hours') >= -24) & (difftime(abx_times[j], bld_times[m], units = 'hours') <=0)) | ((difftime(abx_times[j], bld_times[m], units = 'hours') <= 72)) & (difftime(abx_times[j], bld_times[m], units = 'hours') >= 0)) {
        select[i] <- 1
        bld_abx_early_time[i] <- min(abx_times[j], bld_times[m])
        break
      }
    }
    if (select[i] == 1) {
      break
    }
  }
}

table(select)

ptids_2 <- data.frame(ptids, select, bld_abx_early_time)
names(ptids_2) <- c('id', 'select', 'event.time')
ptids_2$id <- as.character(ptids_2$id)
str(ptids_2)
ptids_selected <- ptids_2$id[ptids_2$select == 1]
pt_selected <- ptids_2[ptids_2$select == 1, ]

sofa_first <- data.frame(infos$id, infos$sofa_first)
name(sofa_first) <- c('id', 'sofa0')
sofa_first_ptids_selected <- sofa_first[sofa_first$id %in% ptids_selected, ]


# chart data for sofa scores
charts <- read.csv("pts_vitals1b.csv", header=TRUE)
#get the icuseq info for bld culture
hospseqs <- read.csv('hospseqs_vitals.csv',header = F)
charts$hospseq <- as.vector(unlist(hospseqs), mode = 'numeric')
sofa_charts <- charts[(charts$hospseq > 0) & (!is.na(charts$itemid)) & (!is.na(charts$charttime)) & (charts$itemid == 20009), ]
sofa_charts$id <- paste(charts$subject_id, charts$hospseq, charts$icustay_seq, sep='#%#')
sofa_charts <- charts[charts$id %in% ptids_selected, ]

sofas <- merge(sofa_charts, sofa_first_ptids_selected, by.x = 'id', by.y = 'id')  
sofas$increase <- sofas$valuenum1 - sofas$sofa0
 # after selecting patients according to bld and abx, now select the patients SOFA increase of at least 2
select <- rep(0, length(pt_selected$id))
for (i in 1:length(pt_selected$id)) {
  pid <- ptids[i]
  sofa_times <- sofas$charttime[sofas$id == pid]
  sofa_increase <- sofas$increase[sofas$id == pid]
  for (j in 1:length(sofa_times)) {
    if (((difftime(sofa_times[j], pt_selected[i,3], units = 'hours') >= -24) & (difftime(sofa_times[j], pt_selected[i,3], units = 'hours') <=0)) & (sofa_increase[j] >= 2)) {
      select[i] <- 1
      
      break
    }
  }
    if (select[i] == 1) {
      break
    }
  }

ptids_3 <- data.frame(pt_selected$id, select)
names(ptids_3) <- c('id', 'select')

ptsids_3_selected <- ptids_3$id[ptids_3$select == 1]
pt_selected_sepsis3 <- ptids_2[ptids_2$id %in% ptsids_3_selected, c(1,3)]

write.csv(pt_selected_sepsis3, file = 'pt_selected_sepsis3.csv')
save(pt_selected_sepsis3, file = 'pt_selected_sepsis3.RData')
