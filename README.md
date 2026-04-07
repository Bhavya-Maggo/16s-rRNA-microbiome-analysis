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

Alpha diversity was assessed using four complementary metrics. Statistical comparisons used the Kruskal-Wallis test with FDR correction.

|            Metric        |      H- value    |       p- value      |    Mean   |
| ------------------------ | ---------------- | --------------------|-----------|
| Shannon                  |      4.98        |      0.08           |    5.12   |
| Observed Features        |      2.92        |      0.23           |    83.2   |
| Faith’s PD               |      3.62        |      0.16           |    6.02   |
| Pielou’s Evenness        |      6.89        |      0.03           |    0.81   |

Microbial communities across mock, early, and late groups consist of roughly similar types and numbers of organisms, with no significant differences in richness or phylogenetic diversity. However, their distribution varies significantly, as reflected by changes in evenness, indicating that the transition from early to late stages is associated with shifts in relative abundance—where some taxa become dominant while others decline—rather than the loss or gain of taxa.

  
* **Beta Diversity**
📌 Bray–Curtis: Clear separation between Early and Late

Early = high variability
Late = tightly clustered

📌 Weighted UniFrac: Strong separation based on phylogeny + abundance

Late samples show high stability

📌 Unweighted UniFrac: Differences in taxa presence/absence

Late samples show consistent membership























