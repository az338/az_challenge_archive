
# ------------------------------------------------------------------------------------
# Description: AZ-Sanger Challenge scoring functions
# Authors: Michael P Menden, Julio Saez-Rodriguez
# Modified by Azedine Zoufir (CMI Cambridge) on 18/11/15
# ------------------------------------------------------------------------------------
if(!require('ROCR',quietly=T)) {
    install.packages('ROCR',repos="http://cran.rstudio.com/")
    library(ROCR,quietly=T)
}

# ------------------------------------------------------------------------------------
# Get observation format for Subchallenge 2
# ------------------------------------------------------------------------------------
getObs_ch2 <- function(ls,threshold = F) {
  combi <- unique(ls$COMBINATION_ID)
  cell <- unique(ls$CELL_LINE)
  mat <- matrix(NA, nrow=length(combi), ncol=length(cell),
                dimnames=list(combi, cell))
  for (i in 1:nrow(ls)) 
    mat[as.character(ls$COMBINATION_ID[i]), 
        as.character(ls$CELL_LINE[i])] <- ls$SYNERGY_SCORE[i]
  
  if (is.numeric(threshold)) {
    mat[mat <= threshold] = 0
    mat[mat > threshold ] = 1
  }
  return(mat)
}
# ------------------------------------------------------------------------------------
# Get prediction format for Subchallenge 2
# ------------------------------------------------------------------------------------
getPred_ch2 <- function(pred) {
  if (all(row.names(pred) == c(1:nrow(pred)))) {
    row.names(pred) = pred[,1]
    pred = pred[,-1]
  }
  pred <- as.matrix(pred) 
  return(pred)
}

# ------------------------------------------------------------------------------------
# Get unsigned score from one dimensional ANOVA
# ------------------------------------------------------------------------------------
getNegLog10pVal_ch2 <- function(fit, obs) {
  s <- 0
  if (!is.na(fit$coefficients[2]) & sum(!is.na(obs)) > 2)
      s <- -log10(anova(fit)['pred','Pr(>F)'])
  return(s)
}
# -----------------------------------------------------------------------
# Scores from confusion Matrix
# -----------------------------------------------------------------------
getPrecision_ch2 <- function(obs, pred, threshold=30) {
  obs <- read.csv(obs)
  obs <- getObs_ch2(obs,threshold)
  
  pred <- read.csv(pred,stringsAsFactors=F,check.names = F)
  pred <- getPred_ch2(pred)  
  pred <- pred[match(row.names(obs),row.names(pred)),]
  pred <- pred[,match(colnames(obs),colnames(pred))]

  #Remove all NA's
  pred <- as.numeric(pred)[!is.na(obs)]
  obs <- as.numeric(obs)[!is.na(obs)]

  preds <- prediction(pred,obs)
  prec <- performance(preds,"prec") #precision (Acc + )
  sens <- performance(preds,"sens") #True positive rate (Sensitivity) (Cov +)
  npv <- performance(preds,"npv") #Negative predictive value (Acc - )
  spec <- performance(preds,"spec") #True negative rate(specificity) (Cov -)
  auc <- performance(preds,"auc") #Area under curve (AUC)
  phi <- performance(preds,"phi") #phi correlation coefficient, (matthews)
  aupr <- performance(preds, "prec", "rec") #Area under precision recall (AUPR)
  
  prec_val <- unlist(prec@y.values)[2]
  sens_val <- unlist(sens@y.values)[2]
  npv_val <- unlist(npv@y.values)[2]
  spec_val <- unlist(spec@y.values)[2]
  auc_val <- unlist(auc@y.values)
  phi_val <- unlist(phi@y.values)[2]
  BAC <- (sens_val + spec_val)/2
  F1 <- 2*preds@tp[[1]][2]/(2*preds@tp[[1]][2] + preds@fn[[1]][2] + preds@fp[[1]][2])
  aupr_val <- unlist(aupr@y.values)[2]

  return(round(c(prec=prec_val,
                 sens = sens_val,
                 npv = npv_val,
                 spec=spec_val,
                 auc=auc_val,
                 phi=phi_val,
                 BAC=BAC,
                 F1=F1,                 
                 aupr=aupr_val),2))

  
}

# ------------------------------------------------------------------------------------
# Get the drug combinations score of Subchallenge 2
# ------------------------------------------------------------------------------------
getOneDimScore_ch2 <- function(obs, pred, confidence="none", topX=10, rows=T) {
  obs <- read.csv(obs)
  obs <- getObs_ch2(obs)
  
  pred <- read.csv(pred,stringsAsFactors=F,check.names = F)
  pred <- getPred_ch2(pred)
  pred <- pred[match(row.names(obs),row.names(pred)),]
  pred <- pred[,match(colnames(obs),colnames(pred))]
  
  n <- ncol(obs)
  if (rows)
    n <- nrow(obs)
  
  s <- c()
  for (i in 1:n) {
    if (rows) {
      fit <- aov(obs[i,] ~ pred[i,])
      nlp <- getNegLog10pVal_ch2(fit,obs[i,])
    } else {
      fit <- aov(obs[,i] ~ pred[,i]) 
      nlp <- getNegLog10pVal_ch2(fit,obs[,i])
    }
    
    sign <- 1
    if (mean(obs[pred==1], na.rm=T) < mean(obs[pred==0], na.rm=T))
      sign <- -1
    
    s <- c(s, sign * nlp)
  }
  
  if (!file.exists(confidence))
    return(round(c(mean=mean(s),
             ste=sd(s)),2))
  
  confidence <- read.csv(confidence,stringsAsFactors=F,check.names = F)
  confidence <- getPred_ch2(confidence)
  confidence <- confidence[match(row.names(obs),row.names(confidence)),]
  confidence <- confidence[,match(colnames(obs),colnames(confidence))]
  
  if (rows) {
    idx <- order(rowSums(confidence), decreasing = T)[1:round(topX * (nrow(confidence) / 100))]
  } else {
    idx <- order(colSums(confidence), decreasing = T)[1:round(topX * (ncol(confidence) / 100))]
  }
  
  return(round(c(mean=mean(s[idx]),
           ste=sd(s[idx])),2))
}

# ------------------------------------------------------------------------------------
# Get the performance score of Subchallenge 2
# ------------------------------------------------------------------------------------
getGlobalScore_ch2 <- function(obs, pred) { 
  obs <- read.csv(obs)
  obs <- getObs_ch2(obs)
  
  pred <- read.csv(pred,stringsAsFactors=F,check.names = F)
  pred <- getPred_ch2(pred)
  pred <- pred[match(row.names(obs),row.names(pred)),]
  pred <- pred[,match(colnames(obs),colnames(pred))]
  
  # regress out combination bias
  cov <- rep(rownames(obs), ncol(obs))
  
  c0 <- rep(rownames(obs), ncol(obs))
  c1 <- as.vector(matrix(colnames(obs), ncol=ncol(obs), nrow=nrow(obs), byrow=T))
  
  obs <- as.vector(obs)
  pred <- as.vector(pred)
  
  # run anove with combination label as covariate
  fit <- aov(obs ~ c0 + c1 + pred)
  pVal <- -log10(anova(fit)['pred','Pr(>F)'])
  
  sign <- 1
  if (mean(obs[pred==1], na.rm=T) < mean(obs[pred==0], na.rm=T))
    sign <- -1
  
  return(round(sign * pVal,2))
}


# ------------------------------------------------------------------------------------
#MAIN part of the program
#calls the above functions and compute the relevant scores for the challenge
# ------------------------------------------------------------------------------------
args = commandArgs(trailingOnly = TRUE)
obs_syn = args[3]
pred_syn = args[1]
conf_mat= args[2]


print('##ANOVA-based evaluation:')
print('')
cat('global score (3-way anova):',getGlobalScore_ch2(obs_syn,pred_syn),'\n')
print('')
print('#1-way anova by drug combination')
top10_cc = getOneDimScore_ch2(obs_syn, pred_syn, confidence=conf_mat, topX=10, rows=T)
top20_cc = getOneDimScore_ch2(obs_syn, pred_syn, confidence=conf_mat, topX=20, rows=T)
top30_cc = getOneDimScore_ch2(obs_syn, pred_syn, confidence=conf_mat, topX=30, rows=T)
all_cc = getOneDimScore_ch2(obs_syn, pred_syn, confidence=conf_mat, topX=100, rows=T)
cat('top10: Mean:',top10_cc[1],'SD:',top10_cc[2],'\n')
cat('top20: Mean:',top20_cc[1],'SD:',top20_cc[2],'\n')
cat('top30: Mean:',top30_cc[1],'SD:',top30_cc[2],'\n')
cat('ALL:  Mean:',all_cc[1],'SD:',all_cc[2],'\n')
print('')
print('#1-way anova by cell-lines')
top10_cl = getOneDimScore_ch2(obs_syn, pred_syn, confidence=conf_mat, topX=10, rows=F)
top20_cl = getOneDimScore_ch2(obs_syn, pred_syn, confidence=conf_mat, topX=20, rows=F)
top30_cl = getOneDimScore_ch2(obs_syn, pred_syn, confidence=conf_mat, topX=30, rows=F)
all_cl = getOneDimScore_ch2(obs_syn, pred_syn, confidence=conf_mat, topX=100, rows=F)
cat('top10: Mean:',top10_cl[2],'SD:',top10_cl[2],'\n')
cat('top20: Mean:',top20_cl[2],'SD:',top20_cl[2],'\n')
cat('top30: Mean:',top30_cl[2],'SD:',top30_cl[2],'\n')
cat('ALL:  Mean:',all_cl[2],'SD:',all_cl[2],'\n')
print('')
print('##Binary evaluation:')
cat('              ',paste(names(getPrecision_ch2(obs_syn,pred_syn,20)),sep=' '),'\n')
cat('Threshold>=20:',getPrecision_ch2(obs_syn,pred_syn,20),'\n')
cat('Threshold>=30:',getPrecision_ch2(obs_syn,pred_syn,30),'\n')
cat('Threshold>=40:',getPrecision_ch2(obs_syn,pred_syn,40),'\n')
