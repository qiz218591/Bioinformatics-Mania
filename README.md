# Bioinformatics-Mania
ParathyroidGeneSE analysis
---
title: "parathyroidgenesSE expression analysis"
author: "Divya"
date: "2022-11-03"
output:
  word_document: default
  html_document: default
  pdf_document: default
editor_options: 
  markdown: 
    wrap: 72
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

%When writing this homework assignment, I did not recall how to %insert
code in a nice looking way into LaTeX documents, %so I referred to this
page on stackoverflow for help:
%<https://stackoverflow.com/questions/3175105>

\usepackage{minted}

\\begin{minted}[mathescape, linenos]{python} Code To Insert in \LaTeX...

```{r}
require(parathyroidSE)
data(parathyroidGenesSE)
parathyroidGenesSE
```

class: RangedSummarizedExperiment dim: 63193 27 metadata(1): MIAME
assays(1): counts rownames(63193): ENSG00000000003 ENSG00000000005 ...
LRG_98 LRG_99 rowData names(0): colnames: NULL colData names(8): run
experiment ... study sample

```{r}
colData(parathyroidGenesSE)
```

DataFrame with 27 rows and 8 columns run experiment patient treatment
time submission <character> <factor> <factor> <factor> <factor> <factor>
1 SRR479052 SRX140503 1 Control 24h SRA051611 2 SRR479053 SRX140504 1
Control 48h SRA051611 3 SRR479054 SRX140505 1 DPN 24h SRA051611 4
SRR479055 SRX140506 1 DPN 48h SRA051611 5 SRR479056 SRX140507 1 OHT 24h
SRA051611 ... ... ... ... ... ... ... 23 SRR479074 SRX140523 4 DPN 48h
SRA051611 24 SRR479075 SRX140523 4 DPN 48h SRA051611 25 SRR479076
SRX140524 4 OHT 24h SRA051611 26 SRR479077 SRX140525 4 OHT 48h SRA051611
27 SRR479078 SRX140525 4 OHT 48h SRA051611 study sample <factor>
<factor> 1 SRP012167 SRS308865 2 SRP012167 SRS308866 3 SRP012167
SRS308867 4 SRP012167 SRS308868 5 SRP012167 SRS308869 ... ... ... 23
SRP012167 SRS308885 24 SRP012167 SRS308885 25 SRP012167 SRS308886 26
SRP012167 SRS308887 27 SRP012167 SRS308887

```{r}
parathyroidGenesSE[, parathyroidGenesSE$patient %in% c(1, 4)]
```

class: RangedSummarizedExperiment dim: 63193 13 metadata(1): MIAME
assays(1): counts rownames(63193): ENSG00000000003 ENSG00000000005 ...
LRG_98 LRG_99 rowData names(0): colnames: NULL colData names(8): run
experiment ... study sample

```{r}
a<- parathyroidGenesSE[, parathyroidGenesSE$patient %in% c(1, 4)]
comparison <- c("a","a","a","a","a","a","b","b","b","b","b","b","b")
colData(a) <- cbind(colData(a), comparison)
a@colData
```

DataFrame with 13 rows and 9 columns run experiment patient treatment
time submission <character> <factor> <factor> <factor> <factor> <factor>
1 SRR479052 SRX140503 1 Control 24h SRA051611 2 SRR479053 SRX140504 1
Control 48h SRA051611 3 SRR479054 SRX140505 1 DPN 24h SRA051611 4
SRR479055 SRX140506 1 DPN 48h SRA051611 5 SRR479056 SRX140507 1 OHT 24h
SRA051611 ... ... ... ... ... ... ... 9 SRR479074 SRX140523 4 DPN 48h
SRA051611 10 SRR479075 SRX140523 4 DPN 48h SRA051611 11 SRR479076
SRX140524 4 OHT 24h SRA051611 12 SRR479077 SRX140525 4 OHT 48h SRA051611
13 SRR479078 SRX140525 4 OHT 48h SRA051611 study sample comparison
<factor> <factor> <character> 1 SRP012167 SRS308865 a 2 SRP012167
SRS308866 a 3 SRP012167 SRS308867 a 4 SRP012167 SRS308868 a 5 SRP012167
SRS308869 a ... ... ... ... 9 SRP012167 SRS308885 b 10 SRP012167
SRS308885 b 11 SRP012167 SRS308886 b 12 SRP012167 SRS308887 b 13
SRP012167 SRS308887 b

```{r}
dds <-DESeqDataSet(a, design = ~patient + treatment)
dds
```

class: DESeqDataSet dim: 63193 13 metadata(2): MIAME version assays(1):
counts rownames(63193): ENSG00000000003 ENSG00000000005 ... LRG_98
LRG_99 rowData names(0): colnames: NULL colData names(9): run experiment
... sample compariso

```{r}
library(DESeq2)
dds <-DESeqDataSet(a, design = ~patient +treatment)
dds
keep <- rowSums(counts(dds)) >=10
dds <- dds[keep,]
```

class: DESeqDataSet dim: 21112 13 metadata(2): MIAME version assays(1):
counts rownames(21112): ENSG00000000003 ENSG00000000005 ...
ENSG00000271662 ENSG00000271707 rowData names(0): colnames: NULL colData
names(9): run experiment ... sample comparison

```{r}
rld <- rlogTransformation(dds)
head(assay(rld))
```

                        1           2         3        4         5

ENSG00000000003 9.2382241 9.380269147 9.2224735 9.351588 9.2933052
ENSG00000000005 0.1639762 0.002252169 0.1233632 0.106903 0.1817621
ENSG00000000419 7.8943268 7.715942075 7.8753845 7.726800 7.8785360
ENSG00000000457 7.0218770 7.033326498 7.0628297 6.896348 6.8088652
ENSG00000000460 8.0872929 7.337137647 8.0137027 7.431175 7.9823095
ENSG00000000938 2.3025199 2.478988367 2.3493355 2.380409 2.1891370 6 7 8
9 ENSG00000000003 9.47232053 8.77070534 8.75825709 8.86150922
ENSG00000000005 0.02053447 -0.04362329 -0.04237947 -0.03164322
ENSG00000000419 7.76126424 7.88231056 7.90575412 7.86074948
ENSG00000000457 6.94642902 7.19043538 7.13039639 7.04520236
ENSG00000000460 7.45551772 6.82786331 8.10499750 7.17106709
ENSG00000000938 2.38651637 2.56103354 2.53611718 2.57237320 10 11 12 13
ENSG00000000003 8.74913888 8.81253928 8.65200908 8.61974382
ENSG00000000005 -0.04730499 -0.03914103 -0.03638684 -0.04919825
ENSG00000000419 7.75035853 7.95750278 7.94088729 7.92511159
ENSG00000000457 7.10257939 7.19089710 7.17432859 7.06191087
ENSG00000000460 7.14932889 8.17863268 7.36924851 7.41912133
ENSG00000000938 2.67850848 2.63142783 2.62349971 2.47134689

```{r}
sampleDists <- dist( t( assay(rld) ) )
sampleDists
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste( colData(rld)$treatment, 
                                     colData(rld)$patient, sep="-" )
library( "gplots" )   
heatmap.2( sampleDistMatrix, trace="none" )
```

###heatmap ![alt text here](F:%5Cdata%5CRplot13.png)

```{r}
print( plotPCA( rld, intgroup = c( "patient", "treatment") ) )
```

##pca ![alt text here](F:%5Cdata%5CRplot14.png)

```{r}
dds <-DESeq(dds)
res <-results(dds)
res
```

log2 fold change (MLE): treatment OHT vs Control Wald test p-value:
treatment OHT vs Control DataFrame with 21112 rows and 6 columns
baseMean log2FoldChange lfcSE stat <numeric> <numeric> <numeric>
<numeric> ENSG00000000003 539.81743 -0.000518099 0.115916 -0.00446963
ENSG00000000005 1.08741 0.247633283 1.476883 0.16767293 ENSG00000000419
231.80566 0.059917287 0.114663 0.52255003 ENSG00000000457 133.44246
-0.150237665 0.130291 -1.15309214 ENSG00000000460 208.82918 0.505333435
0.540814 0.93439396 ... ... ... ... ... ENSG00000271643 24.78963
0.190424 0.364763 0.522048 ENSG00000271644 3.62177 -0.564756 0.741107
-0.762044 ENSG00000271646 9.65884 0.902594 0.533941 1.690436
ENSG00000271662 2.43725 -0.814716 0.899735 -0.905507 ENSG00000271707
6.34875 -0.794836 0.539010 -1.474620 pvalue padj <numeric> <numeric>
ENSG00000000003 0.996434 0.999946 ENSG00000000005 0.866841 0.999946
ENSG00000000419 0.601287 0.999946 ENSG00000000457 0.248873 0.999946
ENSG00000000460 0.350101 0.999946 ... ... ... ENSG00000271643 0.6016368
0.999946 ENSG00000271644 0.4460338 0.999946 ENSG00000271646 0.0909446
0.999946 ENSG00000271662 0.3651969 0.999946 ENSG00000271707 0.1403146
0.999946

```{r}
summary(res)
res <- results(dds, alpha = 0.05, lfcThreshold = 2)
summary(res)
resultsNames(dds)
```

out of 21112 with nonzero total read count adjusted p-value \< 0.05 LFC
\> 2.00 (up) : 0, 0% LFC \< -2.00 (down) : 0, 0% outliers [1] : 2,
0.0095% low counts [2] : 0, 0% (mean count \< 1) [1] see 'cooksCutoff'
argument of ?results [2] see 'independentFiltering' argument of ?results

[1] "Intercept" "patient_4\_vs_1"\
[3] "treatment_DPN_vs_Control" "treatment_OHT_vs_Control"

```{r}
library(EnhancedVolcano)
EnhancedVolcano::EnhancedVolcano(res, lab = rownames(res), x = 'log2FoldChange', y ='pvalue')
```

##volcano plot ![alt text here](F:%5Cdata%5CRplot16.png) ##GO pathways
are in doc file please refer it.
