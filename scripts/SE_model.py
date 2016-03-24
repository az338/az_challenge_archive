
from sklearn import svm, preprocessing
from sklearn.metrics import r2_score,mean_squared_error
from scipy.stats import pearsonr
from sklearn.cross_validation import StratifiedShuffleSplit
from sklearn.feature_selection import VarianceThreshold
from sklearn.ensemble import RandomForestClassifier
import pandas as pd
import numpy as np


DATA_DIR = '/scratch/az338/ucc-fileserver/AZ_challenge_data/'

# load side effect matrix
se_matrix = pd.read_csv(DATA_DIR+'challenge_compounds_SE.csv',index_col=False)

# load syngery training set
syn_training = pd.read_csv(DATA_DIR+'compound_synergy_trainSet_by_disease_area.csv')


# get indexes for which compounds overlap between SE data and synergy data
idx = []
for i, r in syn_training.iterrows():
    mapA = r['COMPOUND_A'] in list(se_matrix['symbol'])
    mapB = r['COMPOUND_B'] in list(se_matrix['symbol'])
    if mapA and mapB : 
        idx.append(i)

# reduce training dataset to combinations for which
# both targets have a SE profile in SE matrix
syn_training = syn_training.iloc[idx]
labels = list(set(syn_training['DISEASE_AREA']))


## Function used to check whether the training data
## for a specific label has the right number of class
def synData_check(train,label) :
    df = train[train['DISEASE_AREA'] == label][['COMPOUND_A','COMPOUND_B','SYNERGY_SCORE']]
    df = df[np.isfinite(df['SYNERGY_SCORE'])]
    df['binary.y'] = df['SYNERGY_SCORE'] >= 10
    nClasses = len(df['binary.y'].value_counts())
    if nClasses >= 2: 
        return label

# get the labels that should be kept (i.e fulfill the above criterion)
keep_label = []
for l in labels:
    tmp = synData_check(syn_training,l)
    if tmp != None:
        keep_label.append(tmp)

# only those labels will be used for the modelling
labels = keep_label
 

# function that builds the training and testing set
# for prediction of target combinations given a specific synergy label
# ARGUMENTS: 
# se_mat : the SE matrix containing the SE profiles for each target
# train : the challenge training dataset containing based on the target (rather than compound) and containing only those combinations for which both targets have a SE profile in the SE matrix 
# label : synergy label to be predicted,the y variable in the ML algorithm  (e.g 'Breast.Synergy' or 'Lung.Synergy' )
# seed (optional) : random seed to be used for splitting the data into training and test set
# returns : predictor training set, predictor test set, y (synergy) training set, y test set
def make_train_test(se_mat,train,label,seed=5):
    df = train[train['DISEASE_AREA'] == label]
    # remove NaN predictor as they cannot be used for performance evaluation 
    df = df[np.isfinite(df['SYNERGY_SCORE'])]  

    # define a binarised version of y as required for the SSS below
    df['binary.y'] = df['SYNERGY_SCORE'] >= 10

    
    # create combined SE profiles matrix to be used as predictors  
    combined_se = [] 
    for i, r in df.iterrows():
        # SE profile for 1st target
        seA = se_mat[se_mat['symbol'] == r['COMPOUND_A']].drop('symbol',1) 
        # SE profile for 2nd target
        seB = se_mat[se_mat['symbol'] == r['COMPOUND_B']].drop('symbol',1)
        #combined_se.append(seA.values/2. + seB.values/2.) # ADD: combined SE profile = 1/2*(se vector target1) + 1/2*(se vector target2) 
        combined_se.append(seA.values*seB.values) # PROD : combined SE profile = se vector target 1 x se vector target 2
        # OR operation:
        #seC = seA.values+seB.values
        #seC[seC > 0] = 1
        #combined_se.append(seC)

    # combined all profiles into matrix of dimension (n targets x n SE) 
    combined_se = np.squeeze(np.array(combined_se),axis=(1,))

    
    # Stratified Shuffle Split with 10% of values to be used as y
    sss = StratifiedShuffleSplit(y=df['binary.y'].values,n_iter=10,test_size=.1,random_state=seed)
    df = df.drop('binary.y',axis=1)
    train_x_list = []
    test_x_list = []
    train_y_list = []
    test_y_list = []
    for train, test in sss:    
        train_x_list.append(combined_se[train,]) 
        test_x_list.append(combined_se[test,] )
        train_y_list.append(df.iloc[train]) 
        test_y_list.append(df.iloc[test])
    return train_x_list, test_x_list, train_y_list, test_y_list



# create train and test set for each label, and test classifier on this
# split for each label 
for l in labels:
    train_x,test_x,train_y,test_y = make_train_test(se_matrix,syn_training,l,1402)
    r = []
    r2 = []
    rmse = []

    # for each iteration of SSS :
    for i in xrange(len(train_x)) :
        # pre-processing : 
        # center-scale predictor variables (X)
    
        # feature selection :
        # remove zero variance predictors/features in training and test set
        fs = VarianceThreshold()
        sel_idx = fs.fit(train_x[i]).get_support(indices=True)
        train_x[i] = train_x[i][:,sel_idx]
        test_x[i] = test_x[i][:,sel_idx]

        # initialize SVM classifier (default params)
        #se_clf = svm.SVC(probability=True)
    
        # initialize Random Forest
        se_clf = RandomForestClassifier()

        # fits the classifier to the training data 
        se_fit = se_clf.fit(train_x[i],train_y[i]['SYNERGY_SCORE'].values)
    
        # predict synergy on the test data
        test_y[i]['predicted'] = se_fit.predict(test_x[i])

        # print test set, predictions, and measures of error
        r.append(round(pearsonr(test_y[i]['SYNERGY_SCORE'],test_y[i]['predicted'])[0],2))
        r2.append(round(r2_score(test_y[i]['SYNERGY_SCORE'],test_y[i]['predicted']),2))
        rmse.append(round(mean_squared_error(test_y[i]['SYNERGY_SCORE'],test_y[i]['predicted']),2))
    #print l,'//','R2:',round(r2_score(test_y['SYNERGY_SCORE'],test_y['predicted']),2),';','RMSE:',round(mean_squared_error(test_y['SYNERGY_SCORE'],test_y['predicted']),2),';','Pearson R:',round(r[0],2)
    #test_y.to_csv('obs_vs_pred_'+l+'.csv')
    #scatter(test_y['SYNERGY_SCORE'],test_y['predicted'])
    print l, 'R:', np.mean(r), 'R2:',np.mean(r2), 'RMSE:',np.mean(rmse)
    
