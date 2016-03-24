folder_data = '/scratch/az338/ucc-fileserver/AZ_challenge_data/'
default_location_synergies = folder_data+"drug_synergy_data/ch1_train_combination_and_monoTherapy.csv"
default_location_cell_info = folder_data+"sanger_molecular_data/cell_info.csv"

default_threshold_synergism = 10 # any synergy above this is labelled as synergistic

import pandas, numpy, warnings
from collections import Counter
from sklearn.cross_validation import StratifiedKFold
from load_synergies import Synergies

class StratifiedKFolds:
    def __init__(self,**args):
        # args: data_synergies/location_synergies, location_cell_info
        self.synergies = Synergies(**args) #unpacks args
        self.location_cell_info = args['location_cell_info'] if 'location_cell_info' in args.keys() else default_location_cell_info 
        self.load_cell_info()
        
    def process_data(self,**args):
        # args: remove_low_quality, select_cell_lines, select_drug_combinations, remove_cell_lines, remove_drug_combinations
        self.synergies.process_data(**args)        
        
    def split(self,**args):
        # args: folds, seed, threshold_synergism
        assert 'folds' in args.keys(), 'Missing argument: folds.'   
        self.folds = args['folds']
        self.seed = args['seed'] if 'seed' in args.keys() else None
        
        self.threshold_synergism = args['threshold_synergism'] if 'threshold_synergism' in args.keys() else default_threshold_synergism
        self.threshold_antagonism = -self.threshold_synergism
        
        self.label_data() # create labels for stratification
        self.split_data() # split data into <self.folds> folds of (training,data) pandas dataframes
        
        
    ''' Methods for loading data '''
    def load_cell_info(self):
        self.raw_cell_info = pandas.read_csv(self.location_cell_info)
        
        
    ''' Methods for splitting data '''
    def label_data(self):
        # Label = <disease>|<synergy_label> (=synergistic/antagonistic/additive)
        self.labels = []
        for (_,_,cell_line,synergy) in self.synergies.as_list():
            cell_info_row = self.raw_cell_info[self.raw_cell_info['Sanger.Name'] == cell_line]
            disease = cell_info_row['Disease.Area'].values[0]
            synergy_label = 'synergistic' if synergy >= self.threshold_synergism else 'antagonistic' if synergy <= self.threshold_antagonism else 'additive'
            label = disease+"|"+synergy_label        
            self.labels.append(label)
        
        print "Categories for stratification:"
        for name, count in Counter(self.labels).items():
            print name, count
        
    def split_data(self):
        # Use the labels and skikit-learn to make stratified folds
        self.training_dataframes = []
        self.test_dataframes = []
        skf = StratifiedKFold(y=self.labels, n_folds=self.folds, shuffle=True, random_state=self.seed)
        for train_indices,test_indices in skf:
            training,test = self.synergies.raw_synergies.iloc[train_indices], self.synergies.raw_synergies.iloc[test_indices]
            self.training_dataframes.append(training)
            self.test_dataframes.append(test)
        self.test_splits()
        
    # Test whether any training splits have no entries for a cell line or drug combination
    def test_splits(self):
        all_splits = self.return_all_splits()
        for f,(train,_) in enumerate(all_splits):
            sums_rows = train.mask.sum(axis=1)
            sums_cols = train.mask.sum(axis=0)
            for i,r in enumerate(sums_rows):
                drug_combination_name = self.synergies.unique_drug_combinations[i]
                if r == 0:
                    print "WARNING: Empty row in fold %s, row index %s. Drug combination name: %s." % (f+1,i,drug_combination_name)
            for j,c in enumerate(sums_cols):
                cell_line_name = self.synergies.unique_cell_lines[j]
                if c == 0:
                    print "WARNING: Empty column in fold %s, column index %s. Cell line name: %s." % (f+1,j,cell_line_name)
            
    
    ''' Methods for returning splits '''
    # Return the ith data split - so the ith (training,test) pandas dataframes
    def return_split_raw(self,index=0):
        return (self.training_dataframes[index],self.test_dataframes[index])
        
    # Return the split but as Synergies objects
    def return_split(self,index=0):
        (train_data,test_data) = self.return_split_raw(index=index)
        args = {'select_cell_lines':self.synergies.unique_cell_lines,'select_drug_combinations':self.synergies.unique_drug_combinations}
        train, test = Synergies(data_synergies=train_data), Synergies(data_synergies=test_data)
        train.process_data(**args), test.process_data(**args)
        return (train,test)
        
    # Return all (training,test) pairs of pandas dataframes
    def return_all_splits_raw(self):
        return zip(self.training_dataframes,self.test_dataframes)
        
    def return_all_splits(self):
        return [self.return_split(i) for i in range(0,self.folds)]
