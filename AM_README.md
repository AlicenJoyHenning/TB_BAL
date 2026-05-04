# Alveolar Macrophages 


## Pre-processing 
Alveolar Macrophages (MARCO+) cells from the BAL samples were isolated and investigated further.   
   

All BAL samples were included in this analysis, despite the presence of batch effects in the BAL samples seen during clustering, due to the high quality of their alignment metrics ([Alignment QC](https://github.com/AlicenJoyHenning/TB_BAL/blob/master/data/alignment.csv)
). Significantly expressed (p_adj<0.05) genes coming from technical differences (pseudogenes, genes with unannotated IDs, sex-linked genes XIST, DDX3Y) were not included in the variable gene expression analysis and clustering in the analysis. 

**Table** AMs isolated from each BAL sample

| Sample ID | AM Cell Number | 
|:----------|:----------------|
| ID_1098   | 577             |
| ID_1376   | 1422         |
| ID_1483   | 713          |
| ID_1523   | 443             |
| ID_1566   | 702        |
| ID_1676   | 241         |

  
## AM subclustering 
The isolated AMs were further analyzed through subclustering with `Seurat`. Three distinct populations were identified, representing different stages of macrophage maturation and activity: 
 
- **Early activated AMs** (MARCO+, ITGAM+) showing high expression of $ITGAM$ ($CD11b$), a classic marker of recruited myeloid cells, alongside activation markers like $CYP1B1$ and $EMP1$. 

   
- **Mature AMs** (MARCO+, ITGAM-) with cells have downregulated recruitment markers in favour of genes associated with maintenance of the alveolar niche, including as $THBS1$ (Thrombospondin 1) and $SIGLEC5$. 


- **Proliferating AMs** (MARCO+, MKI67+) a small population emerging as a bridge between them with definite proliferative marker gene expression. 

<br>

![](https://github.com/AlicenJoyHenning/TB_BAL/raw/master/plots/AM_clusters_features.svg)

![](https://github.com/AlicenJoyHenning/TB_BAL/raw/master/plots/AM_clusters_umap.png)

![ ](https://github.com/AlicenJoyHenning/TB_BAL/raw/master/plots/markers_dot.svg)

**Figure** Alevolar Macrophage (AM) subclusters with marker gene expression 



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


