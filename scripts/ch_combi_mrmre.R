library(plyr)
library(reshape2)
library(ggplot2)
library(gridExtra)
library(mRMRe)

DATA= '~/ucc-fileserver/AZ_challenge_data/'
CMP_DATA = file.path(DATA,'dan_descriptors/')


# combination descriptors
CC_DATA  = file.path(CMP_DATA,'combinations/') 



# here load training combinations are loaded for which we actually have synergy scores 
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
cc_tr_syn = read.csv(file.path(CC_DATA,'dan_combination_train_synergies.csv'))
cc_tr_syn = cc_tr_syn[,-c(2,3,5)]

# concatenate eveything together

allCC = Reduce(function(...) merge(..., all=T, by='ID'), list(cc_tr_SE,cc_tr_other,cc_tr_cell,cc_tr_diseases,cc_tr_annoTargets,
                                                              cc_tr_tissue,cc_tr_predTargets,cc_tr_pathways,cc_tr_STN,cc_tr_fingerprints,
                                                              cc_tr_syn))

allCC2 = data.frame(lapply(allCC, function(x) ordered(x)))[,-1]
FS_data = mRMR.data(data = allCC2)
synLab_idx = grep('Synergy_Label_BINARY',colnames(allCC2))

# ge relevant ids to be 
id_train = match(cc_tr_fingerprints$ID,allCC)
cc_t1_fingerprints = read.csv(file.path(CC_DATA,'dan_combination_test_1_structural_fingerprints.csv'))
cc_t2_fingerprints = read.csv(file.path(CC_DATA,'dan_combination_test_2_structural_fingerprints.csv'))
id_test1 = match(cc_t1_fingerprints$ID,allCC)
id_test2 = match(cc_t2_fingerprints$ID,allCC)

# feature selection
l_ply(c(10,20,50,100,200,500,1000,1500,2000),function(nF) {
  cat(nF,'\n')
  FS = mRMR.classic(FS_data,target_indices=c(synLab_idx),feature_count = nF)
  selFeatures = unlist(FS@filters)
  df = data.frame(ID=allCC$ID,allCC2[,selFeatures])
  write.csv(df,file.path(DATA,paste('az_combination_train_top',nF,'_mrmrSelected.csv',sep='')),row.names=F,quote=F)
  write.csv(df,file.path(DATA,paste('az_combination_test_1_top',nF,'_mrmrSelected.csv',sep='')),row.names=F,quote=F)
  write.csv(df,file.path(DATA,paste('az_combination_test_2_top',nF,'_mrmrSelected.csv',sep='')),row.names=F,quote=F)
})


# some plotting
FS = mRMR.classic(FS_data,target_indices=c(synLab_idx),feature_count = 50)
selFeatures = c(1,unlist(FS@filters)+1)
colnames(allCC2)[selFeatures]
dat = data.frame(ID=allCC$ID,allCC2[,selFeatures])
df = melt(dat,id.vars = 'ID')
df$value[df$value > 1] = 1
df$Syn = allCC2$Synergy_Label_BINARY
df = df[order(df$Syn),]
df$ID = factor(df$ID,levels=unique(df$ID))
p1 = ggplot(df,aes(y=ID,x=variable,fill=as.numeric(value))) +
  geom_raster() +
  scale_fill_gradient(low='beige',high='red') +
  theme(axis.text.y=element_blank())
p2 = ggplot(df,aes(y=ID,x='',fill=factor(Syn)))+
  geom_raster() +
  theme(axis.text.y=element_blank(),
        axis.ticks.x = element_blank())
grid.arrange(p1,p2,nrow=1,widths=c(5/7,2/7))

