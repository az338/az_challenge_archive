# ------------------------------------------------------------------------------------
# Description: AZ-Sanger Challenge scoring functions
# Authors: Michael P Menden, Julio Saez-Rodriguez
# Modified by : Azedine Zoufir (CMI Cmabridge) on 17/11/15
# ------------------------------------------------------------------------------------


# Get observation format for Subchallenge 1
getObs_ch1 <- function(ls) {
  return(data.frame(CELL_LINE=as.character(ls$CELL_LINE),
                    COMBINATION_ID=as.character(ls$COMBINATION_ID),
                    OBSERVATION=ls$SYNERGY_SCORE))
}

# Get the mean correlation and SD of Subchallenge 1 for topX 
getDrugCombiScore_ch1 <- function(obs, pred, confidence=NA, topX=10) {
  R <- c()
  obs <- read.csv(obs,stringsAsFactors = F)
  obs <- getObs_ch1(obs)
  pred <- read.csv(pred,stringsAsFactors=F)
  pred <- pred[match(paste(obs$CELL_LINE,obs$COMBINATION_ID),paste(pred$CELL_LINE,pred$COMBINATION_ID)),]

  pred$COMBINATION_ID <- gsub(" ", "", pred$COMBINATION_ID)
  for (i in as.character(unique(obs$COMBINATION_ID))) {
      R <- c(R, cor(obs[obs$COMBINATION_ID == i, 'OBSERVATION'], 
                    pred[pred$COMBINATION_ID == i, 'PREDICTION']))
  }
  #Make NA's in R = 0
  R[is.na(R)] = 0
  names(R) <- as.character(unique(obs$COMBINATION_ID))
  
  if (!file.exists(confidence))
    return(round(c(mean=mean(R),
             ste=sd(R),
             n=sum(!is.na(R))),2))
  
  confidence <- read.csv(confidence,stringsAsFactors=F)
  confidence <- confidence[match(unique(obs$COMBINATION_ID),confidence$COMBINATION_ID),]
  idx <- order(confidence$CONFIDENCE,decreasing = T)[1:round(topX * (length(R) / 100))]

  return(round(c(mean=mean(R[idx]),
           ste=sd(R[idx]),
           n=sum(!is.na(R[idx]))),2))
}

# ------------------------------------------------------------------------------------
# Get the global score of Subchallenge 1
# ------------------------------------------------------------------------------------
getGlobalScore_ch1 <- function(obs, pred) {
  obs <- read.csv(obs, stringsAsFactors=F)
  obs <- getObs_ch1(obs)
  pred <- read.csv(pred,stringsAsFactors=F)
  pred <- pred[match(paste(obs$CELL_LINE,obs$COMBINATION_ID),paste(pred$CELL_LINE,pred$COMBINATION_ID)),]

  x = obs$OBSERVATION
  y = pred$PREDICTION
  
  agg <- aggregate(OBSERVATION ~ CELL_LINE, obs, median)
  z0 <- agg$OBSERVATION[match(obs$CELL_LINE, agg$CELL_LINE)]
  
  agg <- aggregate(OBSERVATION ~ COMBINATION_ID, obs, median)
  z1 <- agg$OBSERVATION[match(obs$COMBINATION_ID, agg$COMBINATION_ID)]
   
  parCor <- function(u,v,w) {
    numerator = cor(u,v) - cor(u,w) * cor(w,v)
    denumerator = sqrt(1-cor(u,w)^2) * sqrt(1-cor(w,v)^2)
    return(numerator/denumerator)
  }
  
  numerator=parCor(x,y,z1) - parCor(x,z0,z1) * parCor(z0,y,z1)
  denumerator= sqrt(1-parCor(x,z0,z1)^2) * sqrt(1-parCor(z0,y,z1)^2)
  
  # partial out the mean of synergy across cell lines and combinationations
  return(c(score=numerator/denumerator))
}



# ------------------------------------------------------------------------------------
#MAIN part of the program
#calls the above functions and compute the relevant scores for the challenge
# ------------------------------------------------------------------------------------
args = commandArgs(trailingOnly = TRUE)
obs_syn = args[3]
pred_syn = args[1]
pred_priority = args[2]

cat('Global Correlation:', getGlobalScore_ch1(obs_syn,pred_syn),'\n')
cor_all = getDrugCombiScore_ch1(obs_syn, pred_syn, confidence=pred_priority, topX=100)
cor_top30 =  getDrugCombiScore_ch1(obs_syn, pred_syn, confidence=pred_priority, topX=30)
cor_top20 =  getDrugCombiScore_ch1(obs_syn, pred_syn, confidence=pred_priority, topX=20)
cor_top10 =   getDrugCombiScore_ch1(obs_syn, pred_syn, confidence=pred_priority, topX=10)
cat('Mean Correlation(all):', cor_all[1],'\n')
cat('SD Correlation(all):', cor_all[2],'\n')
cat('Mean Correlation(top30):', cor_top30[1],'\n')
cat('SD Correlation(top30):', cor_top30[2],'\n')
cat('Mean Correlation(top20):', cor_top20[1],'\n')
cat('SD Correlation(top20):', cor_top20[2],'\n')
cat('Mean Correlation(top10):', cor_top10[1],'\n')
cat('SD Correlation(top10):', cor_top10[2],'\n')


