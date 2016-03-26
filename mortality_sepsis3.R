#mortality 

library(stringr)
library(lubridate)
library(plyr)
mortalityinfos2 <- read.csv('pticu_infos.csv', header = T)
mortalityinfos2$id <- paste(mortalityinfos2$subject_id, mortalityinfos2$hospital_seq, mortalityinfos2$icustay_seq, sep='#%#')
mortalityinfos2$id2 <- paste(mortalityinfos2$subject_id, mortalityinfos2$hospital_seq, sep='#%#')

head(mortalityinfos2)
#patients expired in the ICU
ptid.icu_mortality <- mortalityinfos2$id[!is.na(mortalityinfos2$icustay_expire_flg) & mortalityinfos2$icustay_expire_flg == 'Y']
ptid.hosp_mortality <- mortalityinfos2$id2[!is.na(mortalityinfos2$hospital_expire_flg) & mortalityinfos2$hospital_expire_flg == 'Y']

# # load('df2_times_2_5_mmol.RData')
# load('df2_times.RData')
# df1.timeofevent <- df2_times
# str(mortalityinfos2)
# head(df1.timeofevent)
# mortalityinfos2$dod <- ymd_hms(mortalityinfos2$dod)
# df1.mortality30 <- merge(df1.timeofevent, mortalityinfos2, by.x = 'id', by.y = 'id', all.x = T) 
# df1.mortality30 <- df1.mortality30[!is.na(df1.mortality30$dod),]
# dim(df1.mortality30)
# df1.mortality30.id <- c()
# for (i in 1:nrow(df1.mortality30)) 
#   if ((!is.na(df1.mortality30$dod[i])) & (!is.na(df1.mortality30$event.time[i])) & difftime(df1.mortality30$dod[i], df1.mortality30$event.time[i], units="hours") < 30*24) {
#     df1.mortality30.id <- c(df1.mortality30.id, as.character(df1.mortality30$id[i]))
#   }

library(stringr)
hospids <- function(ids) {
  ids2 <- lapply(ids, function (x) str_sub(x, 1, -5))
  ids2 <- unique(unlist(ids2))
  return(ids2)
}

patientids <- function(ids) {
  ids2 <- lapply(ids, function (x) unlist(strsplit(x, "[#%#]"))[1])
  ids2 <- unique(unlist(ids2))
  return(ids2)
}

load(file = 'ptids_sepsis3.RData')

data <- ptids_selected
data_id <- data
data_id2 <- hospids(data_id)
data_id3 <- patientids(data_id)


data_icu_mortality <- intersect(data_id, ptid.icu_mortality)
data_hosp_mortality <- intersect(data_id2, ptid.hosp_mortality)
# data.30_mortality <- intersect(clear.id, df1.mortality30.id)

length(data_icu_mortality) / length(data_id)
length(data_hosp_mortality) / length(data_id2)
# length(data.30_mortality) / length(data_id3)



