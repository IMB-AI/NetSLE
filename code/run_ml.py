import pandas as pd
from sklearn.ensemble import *
from sklearn.model_selection import *
from sklearn.metrics import *
import numpy as np
import copy as cp
import os, math
from tqdm import *
import argparse


parser = argparse.ArgumentParser()
parser.add_argument('--pagerank_path', required = False, default = '../result/pagerank_score/SLE_prior_knowledge')
parser.add_argument('--DEG', required = False, default = '../data/DEG/SLE_DEGs.tsv')
args = parser.parse_args()

prefix = os.path.basename(args.pagerank_path)
if(prefix == ''):
	prefix = os.path.basename(args.pagerank_path[:-1])
os.makedirs(os.path.join('../result', 'disease_probability', prefix), exist_ok = True)

if(args.pagerank_path == '../result/pagerank_score/SLE_prior_knowledge' or args.pagerank_path == '../result/pagerank_score/SLE_prior_knowledge/'):
	### Load pagerank score
	baff = pd.read_csv(os.path.join(args.pagerank_path, 'Drug_target.tsv'), sep = '\t')
	gwas = pd.read_csv(os.path.join(args.pagerank_path, 'GWAS.tsv'), sep = '\t')
	gc = pd.read_csv(os.path.join(args.pagerank_path, 'GC.tsv'), sep = '\t')
	ix = pd.read_csv(os.path.join(args.pagerank_path, 'IX.tsv'), sep = '\t')

	all_pagerank_score = cp.deepcopy(baff)
	all_pagerank_score.rename(columns = {'pagerank_score' : 'BAFF'}, inplace = True)
	all_pagerank_score = pd.merge(all_pagerank_score, gwas, on = 'geneSymbol')
	all_pagerank_score.rename(columns = {'pagerank_score' : 'GWAS'}, inplace = True)
	all_pagerank_score = pd.merge(all_pagerank_score, gc, on = 'geneSymbol')
	all_pagerank_score.rename(columns = {'pagerank_score' : 'GC'}, inplace = True)
	all_pagerank_score = pd.merge(all_pagerank_score, ix, on = 'geneSymbol')
	all_pagerank_score.rename(columns = {'pagerank_score' : 'ImmuneXpresso'}, inplace = True)

	feature_list = ['BAFF','GWAS','GC','ImmuneXpresso']
	score_name = 'pSLE'
	biomarker_name = 'NetSLE biomarker'

	### Load differentially expressed genes (DEGs)
	degs = pd.read_csv(args.DEG, sep = '\t')
	deg_list = degs['geneSymbol'].tolist()
else:
	### Load pagerank score
	file_list = os.listdir(args.pagerank_path)
	print('Feature list: ', [x.split('.')[0] for x in file_list])
	for i,file in enumerate(file_list):
		temp = pd.read_csv(os.path.join(args.pagerank_path, file), sep = '\t')
		temp.rename(columns = {'pagerank_score' : file.split('.')[0]}, inplace = True)
		if(i == 0):
			all_pagerank_score = cp.deepcopy(temp)
		else:
			all_pagerank_score = pd.merge(all_pagerank_score, temp, on = 'geneSymbol')

	feature_list = [x.split('.')[0] for x in file_list]
	score_name = 'pDisease'
	biomarker_name = 'Disease biomarker'

	### Load differntially expressed genes (DEGs)
	degs = pd.read_csv(args.DEG, sep = '\t')
	deg_list = degs.iloc[:,0].tolist()

### Make table for machine learning 
all_pagerank_score.loc[all_pagerank_score['geneSymbol'].isin(deg_list) == True, 'DEG'] = 1
all_pagerank_score.loc[all_pagerank_score['DEG'].isnull() == True, 'DEG'] = 0
all_pagerank_score['DEG'] = all_pagerank_score['DEG'].astype(int)


### Construct machine learning classifier (random forest)
rf = RandomForestClassifier(n_estimators = 1000, max_depth = 3)

gene_list = []; proba_list = []; auroc_list = []; auprc_list = []
for i in trange(1000):
	X_train, X_test, y_train, y_test = train_test_split(all_pagerank_score.loc[:,['geneSymbol'] + feature_list], all_pagerank_score['DEG'], stratify = all_pagerank_score['DEG'], random_state = i, test_size = 0.1)
	rf.fit(X_train.loc[:,feature_list], y_train.values.ravel())
	y_proba = rf.predict_proba(X_test.loc[:,feature_list])
	y_pred = rf.predict(X_test.loc[:,feature_list])

	auroc_list.append(roc_auc_score(y_test, y_proba[:,1]))
	precision, recall, threshold = precision_recall_curve(y_test, y_proba[:,1])
	auprc_list.append(auc(recall, precision))

	gene_list.extend(X_test['geneSymbol'].tolist())
	proba_list.extend(list(y_proba[:,1]))

df_perf = pd.DataFrame(data = {'Metric' : [x for x in ['AUROC','AUPRC'] for y in range(1000)], 'Performance' : auroc_list + auprc_list})

print('Baseline of AUPRC: ', all_pagerank_score['DEG'].value_counts()[1] / len(all_pagerank_score))
print(df_perf.groupby('Metric').mean('Performance'))

df_proba = pd.DataFrame(data = {'geneSymbol' : gene_list, score_name : proba_list})
df_proba_mean = df_proba.groupby('geneSymbol').mean(score_name)
df_proba_mean.insert(0, 'geneSymbol', list(df_proba_mean.index))
df_proba_mean = df_proba_mean.reset_index(drop = True)

DF = pd.merge(df_proba_mean, all_pagerank_score.loc[:,['geneSymbol','DEG']], on = 'geneSymbol', how = 'left')
DF = DF.sort_values(by = score_name, ascending = False).reset_index(drop = True)
DF.to_excel(os.path.join('../result', 'disease_probability', prefix, biomarker_name + '.xlsx'), index = False)
