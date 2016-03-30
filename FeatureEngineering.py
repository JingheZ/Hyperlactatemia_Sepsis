__author__ = 'jinghe'

'''
This script is used to process the predictors for prediction on lactate clearance/normalization
cd /media/sf_Box_Sync/Hyperlactemia sepsis project_jinghe/Data/

'''

# import packages
import pandas as pd
import numpy as np
from datetime import datetime
import cPickle as pickle
import PopulationExtraction_hospitalAdm as populations

def lastsValues(d, n=3):
    f = lambda d: d[-1]
    return f


def renameCols(cols, s):
    f = lambda x: str(x) + '-' + str(s)
    cols2 = map(f, cols)
    return cols2



if __name__ == '__main__':

    with open('predictors_sepsis3.pickle', 'rb') as f:
        charts_labs = pickle.load(f)

    # in general, vitals are recorded every one hour; labs are every 1 hours
    # groupby to create a subset of data for each patient


    charts_labs_table_1 = charts_labs.pivot_table(values='valuenum', index='new_id', columns='itemid', aggfunc=lambda d: d[-1:])
    charts_labs_table_2 = charts_labs.pivot_table(values='valuenum', index='new_id', columns='itemid', aggfunc=lambda d: d[-2:-1])
    charts_labs_table_3 = charts_labs.pivot_table(values='valuenum', index='new_id', columns='itemid', aggfunc=lambda d: d[-3:-2])


    charts_labs_table_1.columns = renameCols(charts_labs_table_1.columns, 1)
    charts_labs_table_2.columns = renameCols(charts_labs_table_2.columns, 2)
    charts_labs_table_3.columns = renameCols(charts_labs_table_3.columns, 3)

    charts_labs_table = pd.concat([charts_labs_table_1, charts_labs_table_2], axis=1)
    charts_labs_table = pd.concat([charts_labs_table, charts_labs_table_3], axis=1)

    #merge the