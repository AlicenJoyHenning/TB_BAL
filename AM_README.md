# Alveolar Macrophages 


## Pre-processing 
Alveolar Macrophages (MARCO+) cells from the BAL samples were isolated and investigated further.   
   

All BAL samples were included in this analysis, despite the presence of batch effects in the BAL samples seen during clustering, due to the high quality of their alignment metrics ([Alignment QC](https://github.com/AlicenJoyHenning/TB_BAL/blob/master/data/alignment.csv)
). Significantly expressed (p_adj<0.05) genes coming from technical differences (pseudogenes, genes with unannotated IDs, sex-linked genes XIST, DDX3Y) were not included in the variable gene expression analysis and clustering in the analysis. 

**Table** AMs isolated from each BAL sample

| **Sample ID** | ID_1098 | ID_1376 | ID_1483 | ID_1523 | ID_1566 | ID_1676 |
|:---|:---:|:---:|:---:|:---:|:---:|:---:|
| **AM Cell Number** | 577 | 1422 | 713 | 443 | 702 | 241 |

  
## AM subclustering 
The isolated AMs were further analyzed through subclustering with `Seurat`. Three distinct populations were identified, representing different stages of macrophage maturation and activity: 
 
- **Early activated AMs** (_MARCO_+, _ITGAM+) showing high expression of $ITGAM$ ($CD11b$), a classic marker of recruited myeloid cells, alongside activation markers like $CYP1B1$ and $EMP1$. 

   
- **Mature AMs** (_MARCO_+, ITGAM-) with cells have downregulated recruitment markers in favour of genes associated with maintenance of the alveolar niche, including as $THBS1$ (Thrombospondin 1) and $SIGLEC5$. 

- **Proliferating AMs** (_MARCO_+, MKI67+) a small population emerging as a bridge between them with definite proliferative marker gene expression. 

<br>


![](https://github.com/AlicenJoyHenning/TB_BAL/raw/master/plots/AM_clusters_umap.png){width=80%}

![](https://github.com/AlicenJoyHenning/TB_BAL/raw/master/plots/AM_clusters_features.svg){width=80%}


![](https://github.com/AlicenJoyHenning/TB_BAL/raw/master/plots/markers_dot.svg){width=80%}

**Figure** Alevolar Macrophage (AM) subclusters with marker gene expression 
  
<br>


**Table**  Differential expression of genes for each AM subcluster 

| AM Subset | Gene | $p_{adj}$ | $avg\_log2FC$ |
| :--- | :--- | :--- | :--- |
| **Early Activated** | *CYP1B1* | $4.21 \times 10^{-225}$ | 7.64 |
| | *EMP1* | $6.58 \times 10^{-162}$ | 2.79 |
| | *ITGAM* | $2.43 \times 10^{-180}$ | 1.91 |
| **Mature Resident** | *STON1* | $5.04 \times 10^{-74}$ | 3.60 |
| | *RSPO3* | $1.10 \times 10^{-71}$ | 3.59 |
| | *SLC19A3* | $8.88 \times 10^{-153}$ | 3.13 |
| **Proliferating** | *MKI67* | $0.00$ | 9.70 |
| | *BUB1B* | $9.12 \times 10^{-317}$ | 9.34 |
| | *TOP2A* | $0.00$ | 8.94 |



<br> 

## Subcluster Proportions 
While each BAL sample has a heterogenous colleciton of AM subclusters, BAL_1376 has a definitively higher proportion of early AMs compared to mature AMS.

<br>

![](https://github.com/AlicenJoyHenning/TB_BAL/raw/master/plots/BAL_subtype_proportion.svg){width=80%}

**Figure** Proportion of AM subtypes across BAL samples


**Table** Proporitons of AM subtypes in each BAL sample

| Sample ID | Total AMs | Early AMs (%) | Mature AMs (%) | Proliferating AMs (%) |
|:----------|:----------|:--------------|:---------------|:----------------------|
| ID_1098 | 577       | 62 (10.7%)    | 505 (87.5%)    | 10 (1.7%)             |
| ID_1376 | 1422      | 1009 (71.0%)  | 358 (25.2%)    | 55 (3.9%)             |
| ID_1483 | 713       | 34 (4.8%)     | 655 (91.9%)    | 24 (3.4%)             |
| ID_1523 | 443       | 18 (4.1%)     | 412 (93.0%)    | 13 (2.9%)             |
| ID_1566 | 702       | 41 (5.8%)     | 644 (91.7%)    | 17 (2.4%)             |
| ID_1676 | 241       | 31 (12.9%)    | 197 (81.7%)    | 13 (5.4%)             |


<br>

## Differentiation markers
However, no significant differences in expression were detected across the samples or subclusters according to classical differentiation genes, with the plot plow showing uniform expressions.


![](https://github.com/AlicenJoyHenning/TB_BAL/raw/master/plots/differentiation_features.svg)