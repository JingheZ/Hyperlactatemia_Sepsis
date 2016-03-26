__author__ = 'jinghe'


import pandas as pd

data = pd.read_csv('iv_fluids_volumes.csv')
data_saline = data.loc[data['label'].str.contains('.9% Normal Saline')]
saline_itemids = data_saline['itemid']
saline_itemids.to_csv('saline_itemids.csv', index=False)

# select the normal saline data from the fluids data
data2 = pd.read_csv('infectdysfunc_fluids.csv')
data2_saline = data2.loc[data2['itemid'].isin(saline_itemids)]
data2_saline.to_csv('infectdysfunc_normalsaline.csv', index=False)