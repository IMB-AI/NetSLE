import networkx as nx
import pandas as pd
import numpy as np
import os
import argparse

### Construct network
graph = nx.Graph()
link = pd.read_csv(os.path.join('../data','network','string_0.txt'), sep = '\t', header = None)
graph.add_edges_from(np.array(link))

parser = argparse.ArgumentParser()
parser.add_argument('--knowledge_path', required = False, default = '../data/SLE_prior_knowledge')
args = parser.parse_args()

prefix = os.path.basename(args.knowledge_path)
if(prefix == ''):
	prefix = os.path.basename(args.knowledge_path[:-1])
os.makedirs(os.path.join('../result', 'pagerank_score', prefix), exist_ok = True)

if(args.knowledge_path == '../data/SLE_prior_knowledge' or args.knowledge_path == '../data/SLE_prior_knowledge/'):
	print('### BAFF (TNFSF13B) pagerank')
	baff = {'TNFSF13B' : 1.0}

	baff_score = nx.pagerank_scipy(graph, personalization = baff, max_iter = 10000, alpha = 0.5)
	baff_score = pd.DataFrame(data = {'geneSymbol' : baff_score.keys(), 'pagerank_score' : baff_score.values()})

	baff_score.sort_values(by = 'pagerank_score', ascending = False).to_csv(os.path.join('../result', 'pagerank_score', prefix, 'Drug_target.tsv'), sep = '\t', index = False)


	print('### GWAS pagerank')
	sle_gwas = pd.read_csv(os.path.join(args.knowledge_path, 'SLE_GWAS.txt'), sep = '\t', header = None)
	gwas = dict()
	for i,j in zip(sle_gwas[0].tolist(), sle_gwas[1].tolist()):
		gwas[i] = j

	gwas_score = nx.pagerank_scipy(graph, personalization = gwas, max_iter = 10000, alpha = 0.5)
	gwas_score = pd.DataFrame(data = {'geneSymbol' : gwas_score.keys(), 'pagerank_score' : gwas_score.values()})

	gwas_score.sort_values(by = 'pagerank_score', ascending = False).to_csv(os.path.join('../result', 'pagerank_score', prefix, 'GWAS.tsv'), sep = '\t', index = False)


	print('### GC (Germinal center) pagerank')
	gc_ = pd.read_csv(os.path.join(args.knowledge_path, 'GCall_genelist.txt'), sep = '\t', header = None)
	gc = dict()
	for i in gc_[0].tolist():
		gc[i] = 1.0

	gc_score = nx.pagerank_scipy(graph, personalization = gc, max_iter = 10000, alpha = 0.5)
	gc_score = pd.DataFrame(data = {'geneSymbol' : gc_score.keys(), 'pagerank_score' : gc_score.values()})

	gc_score.sort_values(by = 'pagerank_score', ascending = False).to_csv(os.path.join('../result', 'pagerank_score', prefix, 'GC.tsv'), sep = '\t', index = False)


	print('### Cytokines (ImmuneXpresso) pagerank')
	IX = pd.read_csv(os.path.join(args.knowledge_path, 'Systemic_lupus_erythematosus.csv'))
	ix = dict()
	gene_list = IX.drop_duplicates([' Cytokine Ontology Label'])[' Cytokine Ontology Label'].tolist()
	for gene in gene_list:
		ix[gene] = 1.0

	ix_score = nx.pagerank_scipy(graph, personalization = ix, max_iter = 10000, alpha = 0.5)
	ix_score = pd.DataFrame(data = {'geneSymbol' : ix_score.keys(), 'pagerank_score' : ix_score.values()})

	ix_score.sort_values(by = 'pagerank_score', ascending = False).to_csv(os.path.join('../result', 'pagerank_score', prefix, 'IX.tsv'), sep = '\t', index = False)

else:
	knowledge_list = os.listdir(args.knowledge_path)
	for file in knowledge_list:
		print('###', file, sep = ' ')
		temp = pd.read_csv(os.path.join(args.knowledge_path, file), sep = '\t',  header = None)
		temp_dict = dict()
		if(len(temp.columns) == 1):
			for i in temp[0].tolist():
				temp_dict[i] = 1.0

		else:
			for i,j in zip(temp[0].tolist(), temp[1].tolist()):
				temp_dict[i] = j

		temp_score = nx.pagerank_scipy(graph, personalization = temp_dict, max_iter = 10000, alpha = 0.5)
		temp_score = pd.DataFrame(data = {'geneSymbol' : temp_score.keys(), 'pagerank_score' : temp_score.values()})

		temp_score.sort_values(by = 'pagerank_score', ascending = False).to_csv(os.path.join('../result', 'pagerank_score', prefix, file.split('.')[0] + '.tsv'), sep = '\t', index = False)
