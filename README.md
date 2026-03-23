# ScSpTITH

Here, we introduce the ScSpTITH (Single-Cell and Spatial Intra-/Inter-tumoral Heterogeneity), a novel computational framework for quantifying ITH for both scRNA-seq and spatial transcriptomic data.
&nbsp;

# Description

ScSpTITH serves as a comprehensive computational framework for dissecting tumor heterogeneity through spatial (across distinct anatomical regions or microenvironments), temporal (during disease progression or therapeutic intervention), and population components (among defined cell types, clusters, or patient cohorts).

&nbsp;

# Details

The function `ScSpTITHscore()` is used to calculate ITH score. It supports **two modalities** , containing six parameters:

**Parameters**

+ **`data`** : The input dataset, gene expression matrix with rows as genes and columns as cells/samples.

+ **`meta`** : Metadata containing at minimum. 
  * **`First column`** : cell/spot IDs (must match column names in data). 

  * **`patient_col`** :Patient ID column (required).

  * **`condition_col`** : Condition/cluster column (optional, only required for cluster mode).

+ **`mode`** : Calculation mode
  * **`overall`** : calculate for entire sample. 
  * **`cluster`** : calculate separately for each cluster.

+ **`patient_col`** : Column name in metadata containing patient IDs, default "Patient".

+ **`condition_col`** : Column name in metadata containing cluster/group information, default "condition".

+ **`top_n_genes`** : Number of highly variable genes used for correlation calculation, default 5000.

+ **`sample_fraction`** : Sampling fraction, default 1.0 means use all cells; Recommended: 0.5 (50%) for optimal balance between computational efficiency and statistical power. Users can adjust based on their needs.

  

# Installation

- Users can install the released version of **ScSpTITH** with:
  &nbsp;

```R
if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")

devtools::install_github("WangX-Lab/ScSpTITH")
```

&nbsp;
&nbsp;

# Examples
**Install ScSpTITH**

```R
library(ScSpTITH)
example_file_path <- system.file("extdata", "example.RData", package = "ScSpTITH")
load(example_file_path)
#ls()
#"data"  "meta"
```

**data**

```R
data[1:5,1:5]
```

| row.names    | SS12pt.10x.P1_AAACCTGTCACCTTAT_1 | SS12pt.10x.P1_AAACCTGTCAGTCAGT_1 | SS12pt.10x.P1_AAACCTGTCCAAAGTC_1 | SS12pt.10x.P1_AAACCTGTCCGTCATC_1 | SS12pt.10x.P1_AAACCTGTCTATGTGG_1 |
| ------------ | -------------------------------- | -------------------------------- | -------------------------------- | -------------------------------- | -------------------------------- |
| **AP006222** | 0                                | 1                                | 1                                | 0                                | 0                                |
| **SAMD11**   | 0                                | 0                                | 0                                | 0                                | 0                                |
| **NOC2L**    | 0                                | 0                                | 2                                | 1                                | 0                                |
| **PLEKHN1**  | 0                                | 1                                | 0                                | 1                                | 0                                |
| **HES4**     | 2                                | 0                                | 0                                | 0                                | 3                                |

**meta**

```R
meta[1:5,]
```

| **cell_name**                        | **sample** | cell_type |
| :----------------------------------- | ---------- | :-------: |
| **SS12pt.10x.P1_AAACCTGTCACCTTAT_1** | SyS12pt    | Malignant |
| **SS12pt.10x.P1_AAACCTGTCAGTCAGT_1** | SyS12pt    | Malignant |
| **SS12pt.10x.P1_AAACCTGTCCAAAGTC_1** | SyS12pt    | Malignant |
| **SS12pt.10x.P1_AAACCTGTCCGTCATC_1** | SyS12pt    | Malignant |
| **SS12pt.10x.P1_AAACCTGTCTATGTGG_1** | SyS12pt    | Malignant |

&nbsp;


## Apply ScSpTITH to 'overall' mode

```R
ScSpTITH = ScSpTITHscore(data, meta, mode = "overall",patient_col = "sample",top_n_genes = 5000,sample_fraction = 1.0) # Sampling fraction, default 1.0 means use all cells
```

ScSpTITH

```R
ScSpTITH
```

| **Patient** | ScSpTITHscore |
| ----------- | ------------- |
| **SyS12pt** | 0.6507015     |
| **SyS13**   | 0.6242106     |
| **SyS14**   | 0.6595786     |

&nbsp;

## **Apply ScSpTITH to 'cluster' mode**

```R
ScSpTITH = ScSpTITHscore(data, meta, mode = "cluster",patient_col = "sample",condition_col = 'cell_type',top_n_genes = 5000,sample_fraction = 1.0) # Sampling fraction, default 1.0 means use all cells
```

**ScSpTITH**

```R
ScSpTITH[1:5,]
```

| **Patient** | ITH_Score | Cluster     |
| :---------- | :-------: | ----------- |
| **SyS12pt** | 0.5498507 | Endothelial |
| **SyS12pt** | 0.6216496 | Fibroblast  |
| **SyS12pt** | 0.5799930 | Macrophage  |
| **SyS12pt** | 0.6436497 | Malignant   |
| **SyS13**   | 0.5706511 | Endothelial |

# Contact

E-mail any questions to Xiaosheng Wang (xiaosheng.wang@hotmail.com)
