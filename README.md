# NetSLE
# README for NetSLE

> Before start, create a virtual environment for analysis using anaconda program and ___NetSLE.yml___ file.
```
conda env create -f NetworkModel.yml
```

---------------------------------------
## Table of Contents (Order of code run)
1. Run pagerank.py for assessing the association bewteen differentially expressed genes (DEGs) and prior knowledge related to disease: [___pagerank.py___](#pagerank.py)
2. Run run_ml.py for probability of association with disease by integrating pagerank score from prior knowledge: [___run_ml.py___](#run_ml.py)

* * *
## <a name="pagerank.py"></a> Network propagation
* ### Input:
    knowledge_path (Directory path of disease-related prior knowledge information)
  
```
python knowledge_path ../data/SLE_prior_knowledge
```

* ### Output:
    Disease association score (network propagation score) in ../result/pagerank_score/~
  

* * *
## <a name="run_ml.py"></a> Probability of disease association
* ### Input:
    pagerank_path (Directory path of pagerank score)

    DEG (DEG list)

  
```
python run_ml.py --pagerank_path ../result/pagerank_score/SLE_prior_knowledge --DEG ../data/DEG/SLE_DEGs.tsv
```

* ### Output:
    Probability of disease association in ../result/disease_probability/~
  
