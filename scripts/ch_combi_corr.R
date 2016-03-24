library(plyr)
library(reshape2)
library(ggplot2)

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

# reformat into matrix for correlation computation
allCCmat= as.matrix(allCC[,-1])
mode(allCCmat) = 'numeric'
# remove 0 variance descriptors
allCCmat = allCCmat[,which(unlist(alply(allCCmat,2,sd)) > 0)]

#correlation
cormat = cor(allCCmat,method='spearman')

# hier. clustering
dist = dist(cormat,method = 'binary')
hclust = hclust(dist)
hc = rownames(cormat)[order(hclust$order)]

# correlation df
cordf = melt(cormat)
colnames(cordf) = c('ft1','ft2','cor')

#reorder according to clustering
cordf$ft1 = factor(cordf$ft1,levels=hc)
cordf$ft2 = factor(cordf$ft2,levels=hc)

# heatmap 
ggplot(cordf,aes(ft1,ft2,fill=cor)) + 
  geom_raster() + 
  scale_fill_gradient2(low='blue',mid='beige',high='red')+
  theme(axis.text.x=element_blank(),
        axis.text.y=element_blank())

# get total correlation (sum correlation)
idx = which(colnames(cormat) %in% colnames(cc_SE))
cormat_SE = cormat[idx,idx]
cors = colSums(cormat_SE)

# retain SE features with low total correlation
SE_keep = names(which(cors < 50 & cors > -50))
#SE_keep = sample(SE_keep,round(.15*length(SE_keep)),replace=F)
SE_remove = colnames(cc_SE)[which(!(colnames(cc_SE) %in% SE_keep))]

# heatmap with reduced correlations
ggplot(subset(cordf,!(ft1 %in% SE_remove) & !(ft2 %in% SE_remove)),aes(ft1,ft2,fill=cor)) + 
  geom_raster() + 
  scale_fill_gradient2(low='blue',mid='beige',high='red') +
  theme(axis.text.x=element_blank(),
        axis.text.y=element_blank())

# heatmap of side effects fetures that are kept in the analysis
ggplot(subset(cordf,ft1 %in% SE_keep & ft2 %in% SE_keep),aes(ft1,ft2,fill=cor)) + 
  geom_raster() + 
  scale_fill_gradient2(low='blue',mid='beige',high='red')

dim(subset(cordf,ft1 %in% SE_remove))

# reduce combinations to those uncorrelated features
allCCmat = allCCmat[,!(colnames(allCCmat) %in% SE_remove)]

#subset and write in relevant files
write.csv(allCCmat[id_train,],file.path(DATA,'az_combination_train_allFeatures_uncorrelated.csv'),row.names=F,quote=F)
write.csv(allCCmat[id_test1,],file.path(DATA,'az_combination_test_1_allFeatures_uncorrelated.csv'),row.names=F,quote=F)
write.csv(allCCmat[id_test2,],file.path(DATA,'az_combination_test_2_allFeatures_uncorrelated.csv'),row.names=F,quote=F)
