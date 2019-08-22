SleepingBeautyITRs = Biostrings::readDNAStringSet(filepath = 'rs20190820_SleepingBeautyITRs.fasta', use.names = T)

save(SleepingBeautyITRs, file = 'data/SleepingBeautyITRs.RData')
load('data/SleepingBeautyITRs.RData')
