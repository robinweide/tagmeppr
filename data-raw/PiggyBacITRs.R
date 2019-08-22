PiggyBacITRs = Biostrings::readDNAStringSet(filepath = 'rs20190820_PiggyBacITRs.fasta', use.names = T)

save(PiggyBacITRs, file = 'data/PiggyBacITRs.RData')
load('data/PiggyBacITRs.RData')
