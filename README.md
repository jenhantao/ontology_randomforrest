ontology_randomforrest
======================
Introduction
----------------------
A series of scripts used to produce a heatmap like the one below. The heatmap conveys the relationship between experimentally determined interacting terms in the gene ontology and interacting terms in a random forrest trained to predict interactions between gene pairs in yeast. Data was generated from a series of synthetic lethality experiments in yeast. The x-axis shows the likelihood score of a pair of terms. Given two GO terms and the two sets of genes associated with those two terms A and B, the likelihood score is the fraction of terms that were experimentally determined to have an interaction out of all possible gene pairs formed between A-B and B-A.The Y-axis gives the average distance between two terms in the random forrest. y-axis values and color values for the heat map are log values.
![text](https://github.com/jenhantao/ontology_randomforrest/blob/master/tree_1_log.png?raw=true)

Workflow
----------------------
