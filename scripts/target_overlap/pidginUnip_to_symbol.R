source('/scratch/az338/ucc-fileserver/ucc_az/AZ_challenge/scripts/uniprot_to_symbol.R')

pidginV1_targets = read.table('/scratch/az338/ucc-fileserver/AZ_challenge_data/pidgin_targets_v1.txt',h=T,quote='',comment.char='',sep='\t')
pidginV2_targets = read.table('/scratch/az338/ucc-fileserver/AZ_challenge_data/pidgin_targets_v2.txt',h=T,quote='',comment.char='',sep='\t')


pidginV1_symbol = uniprot_to_symbol(as.character(pidginV1_targets$UNIPROT))
writeLines(pidginV1_symbol,'/scratch/az338/ucc-fileserver/AZ_challenge_data/pidgin_symbol_v1.txt')
pidginV2_symbol = uniprot_to_symbol(as.character(subset(pidginV2_targets,Organism == 'Homo sapiens')$Uniprot))
writeLines(pidginV2_symbol,'/scratch/az338/ucc-fileserver/AZ_challenge_data/pidgin_symbol_v2.txt')
