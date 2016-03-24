source('~/ucc-fileserver/ucc_az/global_tools/uniprot_symbol_conv.R')

# load chembl20 human targets
chembl20_targets = read.table('~/ucc-fileserver/AZ_challenge_data/chembl20_targets.txt',h=F,sep='\t',quote='',comment.char='',stringsAsFactors=F)
colnames(chembl20_targets) = c('UNIPROT','TARGET')
# convert uniprot to symbols
#chembl20_map = uniprot_to_symbol(chembl20_targets$UNIPROT)
#chembl20_symbols = chembl20_map$SYMBOL

# load challenge targets
AZchg_targets = read.csv('~/ucc-fileserver/AZ_challenge_data/drug_synergy_data/Drug_info_release_curated.csv',stringsAsFactor=F)
AZchg_symbols = unlist(lapply(AZchg_targets$Target.Official.Symbol., function(x) strsplit(x,',')))
AZchg_symbols = gsub('^ ','',AZchg_symbols)
AZchg_symbols = gsub(' $','',AZchg_symbols)

AZchg_uniprot = symbol_to_uniprot(AZchg_symbols)
AZchg_uniprot[which(is.na(AZchg_uniprot$UNIPROT)),]
# overlapping analysis
overlap_targets = intersect(chembl20_targets$UNIPROT,AZchg_uniprot$UNIPROT)

# export table with correpsonding uniprot ids
write.table(AZchg_uniprot[match(overlap_targets,AZchg_uniprot$UNIPROT),],'~/ucc-fileserver/AZ_challenge_data/chembl20_challenge_overlaping.txt',row.names=F,col.names=T)
