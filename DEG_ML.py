# -*- coding: utf-8 -*-
"""
Created on Tue Oct  5 20:04:02 2021

@author: merag
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from collections import defaultdict
from tqdm import tqdm
import os
import re
import PIL
import scipy.stats as stats


#Man data analysis
# cells data

matrix = pd.read_csv('results/normalized_data_male.csv',  index_col=0)
matrix.index = [x.upper() for x in matrix.index]
matrix = matrix.transpose()
matrix['sample_submitter_id'] = matrix.index
meta_data = pd.read_csv('data/HTseq/Metadata/clinical.tsv', sep ='\t', index_col=0)
sample_sheet = pd.read_csv('data/HTseq/Metadata/sample.tsv', sep ='\t', index_col=0)
sample_sheet = sample_sheet[['sample_submitter_id', 'case_submitter_id']]
meta_data = meta_data[['case_submitter_id','treatment_or_therapy']]
meta_data = meta_data.merge(sample_sheet, on = 'case_submitter_id', how = 'left')
matrix = matrix.merge(meta_data, on = 'sample_submitter_id', how = 'left')
matrix = matrix.drop(['sample_submitter_id', 'case_submitter_id'], axis = 1)

import h2o
import psutil



h2o.init()
    
        
        
     
tmp =  h2o.H2OFrame(matrix)
tmp['treatment_or_therapy'] = tmp['treatment_or_therapy'].asfactor()
        
t = tmp[0].runif()
        
train = tmp[ t < 0.7 ]
valid = tmp[ 0.7 <= t ] 
        
        
        
        
        
        #Supervisiored
        
        #Grid-share
        #Gradient Boosting Machine (GBM) 
        
        
ntrees_opt = [2, 5, 10]
            
max_depth_opt = [2, 3, 4]
        
learn_rate_opt = [0.1, 0.2]
        
hyper_parameters = {"ntrees": ntrees_opt, "max_depth":max_depth_opt, "learn_rate":learn_rate_opt}
        
from h2o.grid.grid_search import H2OGridSearch
        
from h2o.estimators.gbm import H2OGradientBoostingEstimator
        
        
gbm = H2OGridSearch(H2OGradientBoostingEstimator(distribution="multinomial"), hyper_params=hyper_parameters)
        
gbm.train(x=train.names[1:len(train.names)-2], y=train.names[len(train.names)-1], training_frame=train, validation_frame= valid)
        
gbm.get_grid('logloss', decreasing=False)
        
mod_gbm = gbm.logloss(valid = True)
        
mod_gbm = pd.DataFrame.from_dict(mod_gbm.items())
        
        
        #Generalized Linear Models (GLM)
        
from h2o.estimators.glm import H2OGeneralizedLinearEstimator
        
        
hyper_parameters = {'alpha': [0.01], 'lambda': [1e-5, 1e-6]}
        
        
glm = H2OGridSearch(H2OGeneralizedLinearEstimator(family="binomial"), hyper_params=hyper_parameters)
        
glm.train(x=train.names[1:len(train.names)-2], y=train.names[len(train.names)-1], training_frame=train, validation_frame= valid)
        
glm.get_grid('logloss', decreasing=False)
        
mod_glm = glm.logloss(valid = True)
        
mod_glm = pd.DataFrame.from_dict(mod_glm.items())
        
        
        
        
        
        
    
if min(mod_glm[1]) > min(mod_gbm[1]):   
    best_model = gbm.models[0]
    del(glm, mod_glm, gbm, mod_gbm)
    genes = best_model.varimp(use_pandas=True)
    best_model.train(x=tmp.names[1:len(tmp.names)-2], y=tmp.names[len(tmp.names)-1], training_frame=tmp)
else:
    best_model = glm.models[0]
    del(glm, mod_glm, gbm, mod_gbm)
    genes = best_model.varimp(use_pandas=True)
    best_model.train(x=tmp.names[1:len(tmp.names)-2], y=tmp.names[len(tmp.names)-1], training_frame=tmp)
        
        
        
genes = genes[genes['scaled_importance'] != 0]
    # best_model.varimp_plot(num_of_features = len(genes))
    
    
    #Interaction selecting - p-val and FC
    
tmp3 = pd.DataFrame(tmp.as_data_frame())
test = []
for i in tmp3.columns[tmp3.columns != 'treatment_or_therapy']:
    f,p = stats.ttest_ind(tmp3[i][tmp3['treatment_or_therapy'] == 'yes'], tmp3[i][tmp3['treatment_or_therapy'] == 'no'])
    test.append(p)
      


tmp3 = tmp3.groupby(tmp3['treatment_or_therapy']).mean()
tmp3 = pd.DataFrame(tmp3.transpose())
tmp3['p-val'] = test
tmp3 = tmp3.sort_values(by = 'p-val', ascending = True)
tmp3['rank'] = range(1, len(tmp3['p-val']) + 1)
tmp3['FDR'] = (tmp3['rank'] / (len(tmp3['p-val']) + 1))
tmp3['-log(p-val)'] = -np.log10(tmp3['p-val'])
tmp3['FC'] = ((np.array(tmp3['yes']))/(np.array(tmp3['no'])))-1
tmp3['log(FC)'] = np.log(np.array(tmp3['yes'])) - np.log(np.array(tmp3['no']))
tmp3 = tmp3.sort_values(by = ['p-val'], ascending = True)
tmp3['top100'] = 'blue'
tmp3['top100'][0:100] = 'red'
tmp3['top100'][tmp3['FDR'] >= 0.05] = 'grey'
df_tmp = tmp3[['log(FC)', 'p-val', 'FDR']]
df_tmp['relation'] = df_tmp.index
df_tmp = df_tmp[df_tmp['p-val'] < 0.05]
df_tmp = df_tmp[0:1000]
    
    
import matplotlib.pyplot as plt
plt.style.use('seaborn-whitegrid')
    
    #Vulcano positive plot
plt.scatter(x=tmp3['log(FC)'], y=tmp3['-log(p-val)'], color = tmp3['top100'])
plt.xlabel('log(FC)')
plt.ylabel('-log(p-val)')
plt.title('Volcano plot')
plt.savefig('results/volcano_man.png')
plt.clf()
plt.close()
    
                 
    # clusters_df = pd.concat([clusters_df, cluster_tmp])
 

#reactome





import pandas
import reactome2py
from reactome2py import analysis, content, utils
from collections import defaultdict
from itertools import chain
from operator import methodcaller



#####
df_tmp['path'] = None
    


for l in tqdm(range(0, len(df_tmp.index))):
    idsr = ','.join(df_tmp.index[l])
    resultr = analysis.identifiers(ids=idsr)
    tokenr = resultr['summary']['token']
    
    
    url = 'https://reactome.org/PathwayBrowser/#/DTAB=AN&ANALYSIS=' + tokenr
    
    token_resultr = analysis.token(tokenr, species = ['Homo sapiens', 'Mus musculus'], page_size='-1', page='-1', sort_by='ENTITIES_FDR', 
                                  order='ASC', resource='TOTAL', p_value='0.05', include_disease=True, 
                                  min_entities=None, max_entities=None)
    
   
    pathwaysr = token_resultr['pathways']
    
    pathways = [p['name'] for p in pathwaysr]
    if len(pathways) == 0:
        pathways = ['None']
    df_tmp['path'][l] = pathways
  

import seaborn



df_tmp_up = df_tmp[df_tmp['log(FC)'] > 0]
df_tmp_up = df_tmp_up[df_tmp_up['FDR'] < 0.05]

ab = []
for i in tqdm(df_tmp_up['path'][df_tmp_up['path'] != str(['None'])]):
    ab = ab + eval(str(i))
    
values, counts = np.unique(ab, return_counts=True)      
count = pd.DataFrame({'names':values, 'n':counts})
count = count[count['names'] != 'None']
count = count.sort_values('n', ascending=False)
count = count[0:20]
seaborn.barplot(count['names'], count['n'])
plt.xticks(rotation=90)
plt.ylabel('Number of paths')
plt.title('Pathways')
plt.savefig('results/path_man_up.png',  bbox_inches='tight')
plt.clf()
plt.close()
        


df_tmp_down = df_tmp[df_tmp['log(FC)'] < 0]
df_tmp_down = df_tmp_down[df_tmp_down['FDR'] < 0.05]

ab = []
for i in tqdm(df_tmp_down['path'][df_tmp_down['path'] != str(['None'])]):
    ab = ab + eval(str(i))
    
values, counts = np.unique(ab, return_counts=True)      
count = pd.DataFrame({'names':values, 'n':counts})
count = count[count['names'] != 'None']
count = count.sort_values('n', ascending=False)
count = count[0:20]
seaborn.barplot(count['names'], count['n'])
plt.xticks(rotation=90)
plt.ylabel('Number of paths')
plt.title('Pathways')
plt.savefig('results/path_man_down.png',  bbox_inches='tight')
plt.clf()
plt.close()
     
pd.DataFrame.to_csv(df_tmp, 'results/summary_man_treatment.csv')





# Semi-Unsupervisiored data        
#PCs   
matrix = pd.read_csv('results/normalized_data_male.csv',  index_col=0)
matrix.index = [x.upper() for x in matrix.index]
matrix = matrix.transpose()

tmp = h2o.H2OFrame(matrix)

from h2o.transforms.decomposition import H2OPCA

pca_decomp = H2OPCA(k = 50, transform="NONE", pca_method="Power")

knee = pca_decomp.train(x=tmp.names[1:len(tmp.names)-1], training_frame=tmp)
knee = pd.DataFrame(knee.varimp())
knee = knee.transpose()
sd = list(knee.loc[0,:])
sd = sd[1:len(sd)+1]
sd = [round(float(num)) for num in sd]
plt.scatter(range(1,51), sd)
plt.ylabel('Standard deviation')
plt.xlabel('PCs')
plt.savefig('results/knee_pcs_man.png',  bbox_inches='tight')


pred = pca_decomp.predict(tmp)

plt.scatter(pred.as_data_frame()['PC1'], pred.as_data_frame()['PC2'])
plt.ylabel('PC2 - ' + str(round(knee.loc[2,1], 2)*100) + '%')
plt.xlabel('PC1 - ' + str(round(knee.loc[1,1], 2)*100) + '%')
plt.savefig('results/PC_man.png',  bbox_inches='tight')
plt.clf()
plt.close()

matrix2 = pred.as_data_frame()
matrix2 = matrix2.transpose()
matrix2.columns = matrix.index
matrix2 = matrix2[0:30]



#K-means

tmp = h2o.H2OFrame(matrix2.transpose())


r = tmp[0].runif()

train = tmp[ r < 0.7 ]
valid = tmp[ 0.7 <= r ] 



from h2o.estimators.kmeans import H2OKMeansEstimator


cluster_estimator = H2OKMeansEstimator(k=2, estimate_k=False, standardize=True)

cluster_estimator.train(x= train.names[1:], training_frame=train,  validation_frame= valid) 

pred = cluster_estimator.predict(tmp)

clusters = np.array(pred.as_data_frame())

#Lackof clusters!



#Woman data analysis
# cells data

matrix = pd.read_csv('results/normalized_data_female.csv',  index_col=0)
matrix.index = [x.upper() for x in matrix.index]
matrix = matrix.transpose()
matrix['sample_submitter_id'] = matrix.index
meta_data = pd.read_csv('data/HTseq/Metadata/clinical.tsv', sep ='\t', index_col=0)
sample_sheet = pd.read_csv('data/HTseq/Metadata/sample.tsv', sep ='\t', index_col=0)
sample_sheet = sample_sheet[['sample_submitter_id', 'case_submitter_id']]
meta_data = meta_data[['case_submitter_id','treatment_or_therapy']]
meta_data = meta_data.merge(sample_sheet, on = 'case_submitter_id', how = 'left')
matrix = matrix.merge(meta_data, on = 'sample_submitter_id', how = 'left')
matrix = matrix.drop(['sample_submitter_id', 'case_submitter_id'], axis = 1)

import h2o
import psutil



h2o.init()
    
        
        
     
tmp =  h2o.H2OFrame(matrix)
tmp['treatment_or_therapy'] = tmp['treatment_or_therapy'].asfactor()
        
t = tmp[0].runif()
        
train = tmp[ t < 0.7 ]
valid = tmp[ 0.7 <= t ] 
        
        
        
        
        
        #Supervisiored
        
        #Grid-share
        #Gradient Boosting Machine (GBM) 
        
        
ntrees_opt = [2, 5, 10]
            
max_depth_opt = [2, 3, 4]
        
learn_rate_opt = [0.1, 0.2]
        
hyper_parameters = {"ntrees": ntrees_opt, "max_depth":max_depth_opt, "learn_rate":learn_rate_opt}
        
from h2o.grid.grid_search import H2OGridSearch
        
from h2o.estimators.gbm import H2OGradientBoostingEstimator
        
        
gbm = H2OGridSearch(H2OGradientBoostingEstimator(distribution="multinomial"), hyper_params=hyper_parameters)
        
gbm.train(x=train.names[1:len(train.names)-2], y=train.names[len(train.names)-1], training_frame=train, validation_frame= valid)
        
gbm.get_grid('logloss', decreasing=False)
        
mod_gbm = gbm.logloss(valid = True)
        
mod_gbm = pd.DataFrame.from_dict(mod_gbm.items())
        
        
        #Generalized Linear Models (GLM)
        
from h2o.estimators.glm import H2OGeneralizedLinearEstimator
        
        
hyper_parameters = {'alpha': [0.01], 'lambda': [1e-5, 1e-6]}
        
        
glm = H2OGridSearch(H2OGeneralizedLinearEstimator(family="binomial"), hyper_params=hyper_parameters)
        
glm.train(x=train.names[1:len(train.names)-2], y=train.names[len(train.names)-1], training_frame=train, validation_frame= valid)
        
glm.get_grid('logloss', decreasing=False)
        
mod_glm = glm.logloss(valid = True)
        
mod_glm = pd.DataFrame.from_dict(mod_glm.items())
        
        
        
        
        
        
    
if min(mod_glm[1]) > min(mod_gbm[1]):   
    best_model = gbm.models[0]
    del(glm, mod_glm, gbm, mod_gbm)
    genes = best_model.varimp(use_pandas=True)
    best_model.train(x=tmp.names[1:len(tmp.names)-2], y=tmp.names[len(tmp.names)-1], training_frame=tmp)
else:
    best_model = glm.models[0]
    del(glm, mod_glm, gbm, mod_gbm)
    genes = best_model.varimp(use_pandas=True)
    best_model.train(x=tmp.names[1:len(tmp.names)-2], y=tmp.names[len(tmp.names)-1], training_frame=tmp)
        
        
        
genes = genes[genes['scaled_importance'] != 0]
    # best_model.varimp_plot(num_of_features = len(genes))
    
    
    #Interaction selecting - p-val and FC
    
tmp3 = pd.DataFrame(tmp.as_data_frame())
test = []
for i in tmp3.columns[tmp3.columns != 'treatment_or_therapy']:
    f,p = stats.ttest_ind(tmp3[i][tmp3['treatment_or_therapy'] == 'yes'], tmp3[i][tmp3['treatment_or_therapy'] == 'no'])
    test.append(p)
      


tmp3 = tmp3.groupby(tmp3['treatment_or_therapy']).mean()
tmp3 = pd.DataFrame(tmp3.transpose())
tmp3['p-val'] = test
tmp3 = tmp3.sort_values(by = 'p-val', ascending = True)
tmp3['rank'] = range(1, len(tmp3['p-val']) + 1)
tmp3['FDR'] = (tmp3['rank'] / (len(tmp3['p-val']) + 1))
tmp3['-log(p-val)'] = -np.log10(tmp3['p-val'])
tmp3['FC'] = ((np.array(tmp3['yes']))/(np.array(tmp3['no'])))-1
tmp3['log(FC)'] = np.log(np.array(tmp3['yes'])) - np.log(np.array(tmp3['no']))
tmp3 = tmp3.sort_values(by = ['p-val'], ascending = True)
tmp3['top100'] = 'blue'
tmp3['top100'][0:100] = 'red'
tmp3['top100'][tmp3['FDR'] >= 0.05] = 'grey'
df_tmp = tmp3[['log(FC)', 'p-val', 'FDR']]
df_tmp['relation'] = df_tmp.index
df_tmp = df_tmp[df_tmp['p-val'] < 0.05]
df_tmp = df_tmp[0:1000]
    
    
import matplotlib.pyplot as plt
plt.style.use('seaborn-whitegrid')
    
    #Vulcano positive plot
plt.scatter(x=tmp3['log(FC)'], y=tmp3['-log(p-val)'], color = tmp3['top100'])
plt.xlabel('log(FC)')
plt.ylabel('-log(p-val)')
plt.title('Volcano plot')
plt.savefig('results/volcano_woman.png')
plt.clf()
plt.close()
    
                 
    # clusters_df = pd.concat([clusters_df, cluster_tmp])
 

#reactome





import pandas
import reactome2py
from reactome2py import analysis, content, utils
from collections import defaultdict
from itertools import chain
from operator import methodcaller



#####
df_tmp['path'] = None
    


for l in tqdm(range(0, len(df_tmp.index))):
    idsr = ','.join(df_tmp.index[l])
    resultr = analysis.identifiers(ids=idsr)
    tokenr = resultr['summary']['token']
    
    
    url = 'https://reactome.org/PathwayBrowser/#/DTAB=AN&ANALYSIS=' + tokenr
    
    token_resultr = analysis.token(tokenr, species = ['Homo sapiens', 'Mus musculus'], page_size='-1', page='-1', sort_by='ENTITIES_FDR', 
                                  order='ASC', resource='TOTAL', p_value='0.05', include_disease=True, 
                                  min_entities=None, max_entities=None)
    
   
    pathwaysr = token_resultr['pathways']
    
    pathways = [p['name'] for p in pathwaysr]
    if len(pathways) == 0:
        pathways = ['None']
    df_tmp['path'][l] = pathways
  

import seaborn



df_tmp_up = df_tmp[df_tmp['log(FC)'] > 0]
df_tmp_up = df_tmp_up[df_tmp_up['FDR'] < 0.05]

ab = []
for i in tqdm(df_tmp_up['path'][df_tmp_up['path'] != str(['None'])]):
    ab = ab + eval(str(i))
    
values, counts = np.unique(ab, return_counts=True)      
count = pd.DataFrame({'names':values, 'n':counts})
count = count[count['names'] != 'None']
count = count.sort_values('n', ascending=False)
count = count[0:20]
seaborn.barplot(count['names'], count['n'])
plt.xticks(rotation=90)
plt.ylabel('Number of paths')
plt.title('Pathways')
plt.savefig('results/path_woman_up.png',  bbox_inches='tight')
plt.clf()
plt.close()
        


df_tmp_down = df_tmp[df_tmp['log(FC)'] < 0]
df_tmp_down = df_tmp_down[df_tmp_down['FDR'] < 0.05]

ab = []
for i in tqdm(df_tmp_down['path'][df_tmp_down['path'] != str(['None'])]):
    ab = ab + eval(str(i))
    
values, counts = np.unique(ab, return_counts=True)      
count = pd.DataFrame({'names':values, 'n':counts})
count = count[count['names'] != 'None']
count = count.sort_values('n', ascending=False)
count = count[0:20]
seaborn.barplot(count['names'], count['n'])
plt.xticks(rotation=90)
plt.ylabel('Number of paths')
plt.title('Pathways')
plt.savefig('results/path_woman_down.png',  bbox_inches='tight')
plt.clf()
plt.close()
     
pd.DataFrame.to_csv(df_tmp, 'results/summary_woman_treatment.csv')





# Semi-Unsupervisiored data        
#PCs   
matrix = pd.read_csv('results/normalized_data_female.csv',  index_col=0)
matrix.index = [x.upper() for x in matrix.index]
matrix = matrix.transpose()

tmp = h2o.H2OFrame(matrix)

from h2o.transforms.decomposition import H2OPCA

pca_decomp = H2OPCA(k = 50, transform="NONE", pca_method="Power")

knee = pca_decomp.train(x=tmp.names[1:len(tmp.names)-1], training_frame=tmp)
knee = pd.DataFrame(knee.varimp())
knee = knee.transpose()
sd = list(knee.loc[0,:])
sd = sd[1:len(sd)+1]
sd = [round(float(num)) for num in sd]
plt.scatter(range(1,51), sd)
plt.ylabel('Standard deviation')
plt.xlabel('PCs')
plt.savefig('results/knee_pcs_woman.png',  bbox_inches='tight')


pred = pca_decomp.predict(tmp)

plt.scatter(pred.as_data_frame()['PC1'], pred.as_data_frame()['PC2'])
plt.ylabel('PC2 - ' + str(round(knee.loc[2,1], 2)*100) + '%')
plt.xlabel('PC1 - ' + str(round(knee.loc[1,1], 2)*100) + '%')
plt.savefig('results/PC_woman.png',  bbox_inches='tight')
plt.clf()
plt.close()

matrix2 = pred.as_data_frame()
matrix2 = matrix2.transpose()
matrix2.columns = matrix.index
matrix2 = matrix2[0:30]



#K-means

tmp = h2o.H2OFrame(matrix2.transpose())


r = tmp[0].runif()

train = tmp[ r < 0.7 ]
valid = tmp[ 0.7 <= r ] 



from h2o.estimators.kmeans import H2OKMeansEstimator


cluster_estimator = H2OKMeansEstimator(k=2, estimate_k=False, standardize=True)

cluster_estimator.train(x= train.names[1:], training_frame=train,  validation_frame= valid) 

pred = cluster_estimator.predict(tmp)

clusters = np.array(pred.as_data_frame())




### DEseq reactoma

man_deseq = pd.read_csv('results/DESEQ_male.csv')
man_deseq = man_deseq[man_deseq['pvalue'] < 0.05]
man_deseq = man_deseq.sort_values(by = ['pvalue'], ascending = True)
man_deseq = man_deseq[0:1000]
man_deseq['path'] = None
man_deseq.index = man_deseq['labels']


for l in tqdm(range(0, len(man_deseq.index))):
    idsr = ','.join(man_deseq.index[l])
    resultr = analysis.identifiers(ids=idsr)
    tokenr = resultr['summary']['token']
    
    
    url = 'https://reactome.org/PathwayBrowser/#/DTAB=AN&ANALYSIS=' + tokenr
    
    token_resultr = analysis.token(tokenr, species = ['Homo sapiens', 'Mus musculus'], page_size='-1', page='-1', sort_by='ENTITIES_FDR', 
                                  order='ASC', resource='TOTAL', p_value='0.05', include_disease=True, 
                                  min_entities=None, max_entities=None)
    
   
    pathwaysr = token_resultr['pathways']
    
    pathways = [p['name'] for p in pathwaysr]
    if len(pathways) == 0:
        pathways = ['None']
    man_deseq['path'][l] = pathways
  

import seaborn



df_tmp_up = man_deseq[man_deseq['log2FoldChange'] > 0]
df_tmp_up = df_tmp_up[df_tmp_up['padj'] < 0.05]

ab = []
for i in tqdm(df_tmp_up['path'][df_tmp_up['path'] != str(['None'])]):
    ab = ab + eval(str(i))
    
values, counts = np.unique(ab, return_counts=True)      
count = pd.DataFrame({'names':values, 'n':counts})
count = count[count['names'] != 'None']
count = count.sort_values('n', ascending=False)
count = count[0:20]
seaborn.barplot(count['names'], count['n'])
plt.xticks(rotation=90)
plt.ylabel('Number of paths')
plt.title('Pathways')
plt.savefig('results/path_DESEQman_up.png',  bbox_inches='tight')
plt.clf()
plt.close()
        


df_tmp_down = man_deseq[man_deseq['log2FoldChange'] < 0]
df_tmp_down = df_tmp_down[df_tmp_down['padj'] < 0.05]

ab = []
for i in tqdm(df_tmp_down['path'][df_tmp_down['path'] != str(['None'])]):
    ab = ab + eval(str(i))
    
values, counts = np.unique(ab, return_counts=True)      
count = pd.DataFrame({'names':values, 'n':counts})
count = count[count['names'] != 'None']
count = count.sort_values('n', ascending=False)
count = count[0:20]
seaborn.barplot(count['names'], count['n'])
plt.xticks(rotation=90)
plt.ylabel('Number of paths')
plt.title('Pathways')
plt.savefig('results/path_DESEQman_down.png',  bbox_inches='tight')
plt.clf()
plt.close()




#woman


woman_deseq = pd.read_csv('results/DESEQ_female.csv')
woman_deseq = woman_deseq[woman_deseq['pvalue'] < 0.05]
woman_deseq = woman_deseq.sort_values(by = ['pvalue'], ascending = True)
woman_deseq = woman_deseq[0:1000]
woman_deseq['path'] = None
woman_deseq.index = woman_deseq['labels']


for l in tqdm(range(0, len(woman_deseq.index))):
    idsr = ','.join(woman_deseq.index[l])
    resultr = analysis.identifiers(ids=idsr)
    tokenr = resultr['summary']['token']
    
    
    url = 'https://reactome.org/PathwayBrowser/#/DTAB=AN&ANALYSIS=' + tokenr
    
    token_resultr = analysis.token(tokenr, species = ['Homo sapiens', 'Mus musculus'], page_size='-1', page='-1', sort_by='ENTITIES_FDR', 
                                  order='ASC', resource='TOTAL', p_value='0.05', include_disease=True, 
                                  min_entities=None, max_entities=None)
    
   
    pathwaysr = token_resultr['pathways']
    
    pathways = [p['name'] for p in pathwaysr]
    if len(pathways) == 0:
        pathways = ['None']
    woman_deseq['path'][l] = pathways
  

import seaborn



df_tmp_up = woman_deseq[woman_deseq['log2FoldChange'] > 0]
df_tmp_up = df_tmp_up[df_tmp_up['padj'] < 0.05]

ab = []
for i in tqdm(df_tmp_up['path'][df_tmp_up['path'] != str(['None'])]):
    ab = ab + eval(str(i))
    
values, counts = np.unique(ab, return_counts=True)      
count = pd.DataFrame({'names':values, 'n':counts})
count = count[count['names'] != 'None']
count = count.sort_values('n', ascending=False)
count = count[0:20]
seaborn.barplot(count['names'], count['n'])
plt.xticks(rotation=90)
plt.ylabel('Number of paths')
plt.title('Pathways')
plt.savefig('results/path_DESEQwoman_up.png',  bbox_inches='tight')
plt.clf()
plt.close()
        


df_tmp_down = woman_deseq[woman_deseq['log2FoldChange'] < 0]
df_tmp_down = df_tmp_down[df_tmp_down['padj'] < 0.05]

ab = []
for i in tqdm(df_tmp_down['path'][df_tmp_down['path'] != str(['None'])]):
    ab = ab + eval(str(i))
    
values, counts = np.unique(ab, return_counts=True)      
count = pd.DataFrame({'names':values, 'n':counts})
count = count[count['names'] != 'None']
count = count.sort_values('n', ascending=False)
count = count[0:20]
seaborn.barplot(count['names'], count['n'])
plt.xticks(rotation=90)
plt.ylabel('Number of paths')
plt.title('Pathways')
plt.savefig('results/path_DESEQwoman_down.png',  bbox_inches='tight')
plt.clf()
plt.close()     