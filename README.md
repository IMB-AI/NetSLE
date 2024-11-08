# NetSLE
# README for NetSLE

> Before start, create a virtual environment for analysis using anaconda program and ___NetSLE.yml___ file. Then, unzip the network file ("string_0.txt.zip") in ./data/network.
```
conda env create -f NetSLE.yml
conda activate NetSLE
cd ./data/network
unzip string_0.txt.zip
```

---------------------------------------
## Table of Contents (Order of code run)
1. Run pagerank.py for assessing the association bewteen differentially expressed genes (DEGs) and prior knowledge related to disease: [___pagerank.py___](#pagerank.py)
2. Run run_ml.py for probability of association with disease by integrating pagerank score from prior knowledge: [___run_ml.py___](#run_ml.py)
#### By preparing prior knowledge data (refer to the format of "SLE_GWAS.txt" or "GCall_genelist.txt") and DEG information (refer to the format of "SLE_DEGs.tsv") for a specific disease in the required file format and organizing them in the correct directory, you can prioritize the biomarker list for that disease. For SLE cases, use the data and command below.

* * *
## <a name="pagerank.py"></a> Network propagation

Runtime depending on the number of prior knowledge (for SLE, ~1 minute).
* ### Input:
    knowledge_path (Directory path of disease-related prior knowledge information)
  
```
cd ../../code
python pagerank.py --knowledge_path ../data/SLE_prior_knowledge
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
  
Cite the code: [![DOI](https://zenodo.org/badge/844365006.svg)](https://doi.org/10.5281/zenodo.13939175)
