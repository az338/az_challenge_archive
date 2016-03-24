library(plyr)
library(ade4)
library(reshape2)
library(ggplot2)
library(GGally)

DATA= '~/ucc-fileserver/AZ_challenge_data/'
CMP_DATA = file.path(DATA,'dan_descriptors/')


# combination descriptors
CC_DATA  = file.path(CMP_DATA,'combinations/') 

cc_t1_SE = read.csv(file.path(DATA,'az_combination_test_1_side_effects.csv'))
cc_t1_other = read.csv(file.path(CC_DATA,'dan_combination_test_1_additional_descriptors.csv'))
cc_t1_cell = read.csv(file.path(CC_DATA,'dan_combination_test_1_annotated_cell_line.csv'))
cc_t1_diseases = read.csv(file.path(CC_DATA,'dan_combination_test_1_annotated_disease.csv'))
cc_t1_annoTargets = read.csv(file.path(CC_DATA,'dan_combination_test_1_annotated_targets.csv'))
cc_t1_tissue = read.csv(file.path(CC_DATA,'dan_combination_test_1_annotated_tissue.csv'))
cc_t1_predTargets = read.csv(file.path(CC_DATA,'dan_combination_test_1_predicted_targets.csv'))
cc_t1_pathways = read.csv(file.path(CC_DATA,'dan_combination_test_1_signalink_pathways.csv'))
cc_t1_STN = read.csv(file.path(CC_DATA,'dan_combination_test_1_signalink_targets_and_neighbours.csv'))
cc_t1_fingerprints = read.csv(file.path(CC_DATA,'dan_combination_test_1_structural_fingerprints.csv'))
#cc_t1_syn = read.csv(file.path(CC_DATA,'dan_combination_test_1_synergies.csv'))


cc_t2_SE = read.csv(file.path(DATA,'az_combination_test_2_side_effects.csv'))
cc_t2_other = read.csv(file.path(CC_DATA,'dan_combination_test_2_additional_descriptors.csv'))
cc_t2_cell = read.csv(file.path(CC_DATA,'dan_combination_test_2_annotated_cell_line.csv'))
cc_t2_diseases = read.csv(file.path(CC_DATA,'dan_combination_test_2_annotated_disease.csv'))
cc_t2_annoTargets = read.csv(file.path(CC_DATA,'dan_combination_test_2_annotated_targets.csv'))
cc_t2_tissue = read.csv(file.path(CC_DATA,'dan_combination_test_2_annotated_tissue.csv'))
cc_t2_predTargets = read.csv(file.path(CC_DATA,'dan_combination_test_2_predicted_targets.csv'))
cc_t2_pathways = read.csv(file.path(CC_DATA,'dan_combination_test_2_signalink_pathways.csv'))
cc_t2_STN = read.csv(file.path(CC_DATA,'dan_combination_test_2_signalink_targets_and_neighbours.csv'))
cc_t2_fingerprints = read.csv(file.path(CC_DATA,'dan_combination_test_2_structural_fingerprints.csv'))
#cc_t2_syn = read.csv(file.path(CC_DATA,'dan_combination_test_2_synergies.csv'))


cc_tr_SE = read.csv(file.path(DATA,'az_combination_train_side_effects.csv'))
cc_tr_other = read.csv(file.path(CC_DATA,'dan_combination_train_additional_descriptors.csv'))
cc_tr_cell = read.csv(file.path(CC_DATA,'dan_combination_train_annotated_cell_line.csv'))
cc_tr_diseases = read.csv(file.path(CC_DATA,'dan_combination_train_annotated_disease.csv'))
cc_tr_annoTargets = read.csv(file.path(CC_DATA,'dan_combination_train_annotated_targets.csv'))
cc_tr_tissue = read.csv(file.path(CC_DATA,'dan_combination_train_annotated_tissue.csv'))
cc_tr_predTargets = read.csv(file.path(CC_DATA,'dan_combination_train_predicted_targets.csv'))
cc_tr_pathways = read.csv(file.path(CC_DATA,'dan_combination_train_signalink_pathways.csv'))
cc_tr_STN = read.csv(file.path(CC_DATA,'dan_combination_train_signalink_targets_and_neighbours.csv'))
cc_tr_fingerprints = read.csv(file.path(CC_DATA,'dan_combination_train_structural_fingerprints.csv'))
#cc_tr_syn = read.csv(file.path(CC_DATA,'dan_combination_train_synergies.csv'))


# concatenate eveything together

cc_SE = do.call('rbind',list(cc_t1_SE,cc_t2_SE,cc_tr_SE))
cc_SE = cc_SE[which(!duplicated(cc_SE)),]
cc_other = do.call('rbind',list(cc_t1_other,cc_t2_other,cc_tr_other))
cc_other = cc_other[which(!duplicated(cc_other)),]
cc_cell = do.call('rbind',list(cc_t1_cell,cc_t2_cell,cc_tr_cell))
cc_cell = cc_cell[which(!duplicated(cc_cell)),]
cc_diseases = do.call('rbind',list(cc_t1_diseases,cc_t2_diseases,cc_tr_diseases))
cc_diseases = cc_diseases[which(!duplicated(cc_diseases)),]
cc_annoTargets = do.call('rbind',list(cc_t1_annoTargets,cc_t2_annoTargets,cc_tr_annoTargets))
cc_annoTargets = cc_annoTargets[which(!duplicated(cc_annoTargets)),]
cc_tissue = do.call('rbind',list(cc_t1_tissue,cc_t2_tissue,cc_tr_tissue))      
cc_tissue = cc_tissue[which(!duplicated(cc_tissue)),]
cc_predTargets = do.call('rbind',list(cc_t1_predTargets,cc_t2_predTargets,cc_tr_predTargets))
cc_predTargets = cc_predTargets[which(!duplicated(cc_predTargets)),]
cc_pathways = do.call('rbind',list(cc_t1_pathways,cc_t2_pathways,cc_tr_pathways))
cc_pathways = cc_pathways[which(!duplicated(cc_pathways)),]
cc_STN = do.call('rbind',list(cc_t1_STN,cc_t2_STN,cc_tr_STN))
cc_STN = cc_STN[which(!duplicated(cc_STN)),]
colnames(cc_STN)[-1] = paste('STN_',colnames(cc_STN)[-1],sep='')
cc_fingerprints = do.call('rbind',list(cc_t1_fingerprints,cc_t2_fingerprints,cc_tr_fingerprints))
cc_fingerprints = cc_fingerprints[which(!duplicated(cc_fingerprints)),]

allCC = Reduce(function(...) merge(..., all=T, by='ID'), list(cc_SE,cc_other,cc_cell,cc_diseases,cc_annoTargets,
                                                              cc_tissue,cc_predTargets,cc_pathways,cc_STN,cc_fingerprints))

# allCC = Reduce(function(...) merge(..., all=T, by='ID'), list(cc_other,cc_cell,cc_diseases,cc_annoTargets,
#                                                                cc_tissue,cc_predTargets,cc_pathways,cc_STN, cc_fingerprints))
# 

# replace NA by non informative point (i.e mean in centered PCA)
allCC_noNA = apply(allCC[,-1],2,function(x) {
  m <- try(mean(as.numeric(x), na.rm = TRUE))
  x[is.na(x)] <- m
  return(x)
})

# PCA
pca = dudi.pca(allCC_noNA ,center = T,scale = T,scannf=F,nf=15)
cumsum(round(pca$eig/sum(pca$eig),3)*100)[1:15] #variance explained is very poor (13% 1st PC, 7% 2nd PC, cumsum first 15: 54%)
pca$li$ID = allCC$ID
pca$li = pca$li[,c('ID',colnames(pca$li)[1:15])]

id_train = match(cc_tr_fingerprints$ID,pca$li$ID)
write.csv(pca$li[id_train,],file.path(DATA,'az_combination_train_pc15.csv'),row.names=F,quote=F)

id_test1 = match(cc_t1_fingerprints$ID,pca$li$ID)
write.csv(pca$li[id_test1,],file.path(DATA,'az_combination_test_1_pc15.csv'),row.names=F,quote=F)

id_test2 = match(cc_t2_fingerprints$ID,pca$li$ID)
write.csv(pca$li[id_test2,],file.path(DATA,'az_combination_test_2_pc15.csv'),row.names=F,quote=F)

# Below not very useful given poor amount variance explained
pca$co$class = NA
pca$co$class[which(!is.na(match(rownames(pca$co),colnames(cc_SE))))] = 'Side effects'
pca$co$class[which(!is.na(match(rownames(pca$co),colnames(cc_other))))] = 'Other'
pca$co$class[which(!is.na(match(rownames(pca$co),colnames(cc_cell))))] = 'Cell'
pca$co$class[which(!is.na(match(rownames(pca$co),colnames(cc_diseases))))] = 'Diseases'
pca$co$class[which(!is.na(match(rownames(pca$co),colnames(cc_annoTargets))))] = 'Annotated targets'
pca$co$class[which(!is.na(match(rownames(pca$co),colnames(cc_tissue))))] = 'Tissue'
pca$co$class[which(!is.na(match(rownames(pca$co),colnames(cc_predTargets))))] = 'Predicted targets'
pca$co$class[which(!is.na(match(rownames(pca$co),colnames(cc_pathways))))] = 'Pathways'
pca$co$class[which(!is.na(match(rownames(pca$co),colnames(cc_STN))))] = 'STN'
pca$co$class[which(!is.na(match(rownames(pca$co),colnames(cc_fingerprints))))] = 'FP'
ggpairs(pca$co, columns=1:5, colour= 'class')
ggpairs(subset(pca$co,!(class %in% c('FP','Side effects'))), columns=1:5, colour= 'class')
