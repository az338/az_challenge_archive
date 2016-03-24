import pandas as pd

class Features:
    
    def __init__(self,location):
        # Load in the csv file stored at location. Use pandas to get a DataFrame.
        # Make sure rows are sorted alphabetically before continuing
        #self.dataframe = self.dataframe.sort(['name_of_row_names_column'])
        #self.location = location
        self.dataframe = pd.read_csv(location)
        self.name = self.dataframe.columns[0]
        self.dataframe = self.dataframe.sort_index(by=self.name,ascending=1)
        return
        
    def process_data(self,filter_rows=None):
        # Transform the pandas DataFrame into:
        # - matrix (rows=datapoints, columns=features)
        # - dictionary of lists, mapping row name to a list of feature values
        # - dictionary of DataFrames/Series, mapping row name to a pandas Series (i.e. that row in self.dataframe)
        # Before that, if argument filter_rows is not None but instead a list 
        # of row names, select only the rows with those names.
        # This allows us to easily filter out datapoints.
        self.filter_rows = filter_rows

        if filter_rows != None:
            self.dataquery = self.dataframe[self.dataframe[self.name].isin(filter_rows)]
        else:
            self.dataquery = self.dataframe
           
        self.matrix = self.dataframe.values[:,1:]    
        self.matrixquery = self.dataquery.values[:,1:]
        return   
    
    def as_name(self):
        return self.name    
    
    def as_dataframe(self):
        return self.dataframe
        
    def as_dataquery(self):
        return self.dataquery
        
    def as_matrix(self):
       return self.matrix

    def as_matrixquery(self):
       return self.matrixquery
        
    def as_dictionary_list(self):
        # Return dictionary from row name to list of feature values
        return
    
    def as_dictionary_dataframe(self):
        # Return dictionary from row name to pandas DataFrame
        return