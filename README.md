# **16S rRNA Microbial Community Analysis**

This repository contains comprehensive analysis of 16S rRNA sequencing data to investigate microbial community dynamics across different stages. The study focuses on comparing microbial composition and diversity across three groups:

* Early
* Late
* Mock (control)

The study integrates alpha diversity and beta diversity (PCoA) analyses using multiple distance metrics.

## **Repository Structure**

```bash
microbiome-analysis/
├── data/
├── scripts/
├── results/
└── README.md

```

## **Methods**

* **Alpha Diversity Metrics**

** Observed Features (Richness)
<br>
** Shannon Diversity
<br>
** Faith’s Phylogenetic Diversity
<br>
** Pielou’s Evenness

* **Beta Diversity Metrics**
  
** Bray–Curtis (abundance-based)
<br>
** Weighted UniFrac (abundance + phylogeny)
<br>
** Unweighted UniFrac (presence/absence + phylogeny)

* **Tools Used**

** QIIME2
<br>
** DADA2
<br>
** R (phyloseq, ggplot2)


## **Results**

* **Alpha Diversity**

Kruskal-Wallis test:

|            Metric        |      H- value    |       p- value      |
| ------------------------ | ---------------- | --------------------|
| Shannon                  |      4.98        |      0.08           |
| Observed Features        |      2.92        |      0.23           |
| Faith’s PD               |      3.62        |      0.16           |
| Pielou’s Evenness        |      6.89        |      0.03           |

## Key Insight:

* Species richness remains stable
* Community evenness significantly changes

  
* **Beta Diversity**
📌 Bray–Curtis: Clear separation between Early and Late

Early = high variability
Late = tightly clustered

📌 Weighted UniFrac: Strong separation based on phylogeny + abundance

Late samples show high stability

📌 Unweighted UniFrac: Differences in taxa presence/absence

Late samples show consistent membership























