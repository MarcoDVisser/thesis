Trade-offs, enemies & dispersal: cross-scale comparisons on tropical tree populations
=====

TeX code, R scripts and C++ programs associated with my PhD Thesis. 

### Quicklinks

-   [Abstract](#abstract)
-   [Contents](#contents)
-   [Supplemental Materials](#supplemental-materials)
-   [Simulation code](#simulation-code)
  

## Abstract

Tropical forests are so rich in tree species that we don't even know how may tree species exist, and identifying the mechanisms structure this diversity has long challenged ecologists. Ample interest in the subject has generated a vast body of work on the topic of diversity maintenance. However, previous work has almost exclusively focused on the fine spatial scales and limited levels of biological organization - such as a single trophic level and life stage.  Yet do processes relevant at a single spatial, temporal or organizational scale translate across scales? The overarching theme of this thesis is to scale up processes across spatial and organizational scales, and to evaluate whether and how processes at fine scales translate to higher levels of organization.

Investigations focus on Barro Colorado Island in the Republic of Panama, one of the most intensely studied tropical forests in existence.  In general, 8 chapters integrate several large-scale datasets looking at ecological patterns, such as negative density dependence and trade-offs, over many hectares of forests and across multiple life stages and vital rates. Two general scientific approaches are used. The first is a model system of the palm Attalea butyracea, studied on large spatial scales and for all life-stages (Chapters 2-4).  The second approach is a cross-species comparison, where 7 datasets are integrated that included millions of observations, spanning three decades, from seeds to forest giants (Chapters 5-8).   All chapters taken together demonstrate that scale can have a profound influence on how we interpret results and understand links between processes operating at different levels (reviewed in chapter 10).  Clear patterns at one scale may crumble to noise at another, only a handful of lower-level and fine scale processes may be influential at larger scales, and that abstraction toward larger scales tends to increase predictably.

## Contents

#### Chapter 1
An exercise in cross-scale integration

#### Chapter 2
Tri-trophic interactions affect density dependence of seed fate in a tropical forest palm. 

#### Chapter 3
Negative density-dependence of seed dispersal and seedling recruitment in a
Neotropical palm

#### Chapter 4
Population-level density dependence in a tropical forest:  Regulation and limitation of a common palm

#### Chapter 5
Functional traits as predictors of vital rates across the life-cycle of tropical trees

#### Chapter 6
Surviving in a cosexual world: a cost-benefit analysis of dioecy in tropical trees

#### Chapter 7
Lianas differentially impact population growth rates of tropical tree species

#### Chapter 8
Explaining variation among tree species in liana infestation 

#### Chapter 9
Speeding up ecological and evolutionary computations in R; essentials of high performance computing for biologists.

#### Chapter 10
Increasing the scale of inquiry (General discussion).

## Supplemental Materials
Supplemental materials, organized per chapter, can be found in this [subfolder](https://github.com/MarcoDVisser/thesis/tree/master/SupplementalMaterial).

## Simulation code
An adaptation of the simple patch model by Pacala 1997 in R, with parts refactored in C++ using Rcpp, is given in this [subfolder](https://github.com/MarcoDVisser/thesis/tree/master/R). This (simple) model explores the effects of negative density dependent dispersal as documented by Jansen et al. 2014 on species richness. 

#### References

Jansen PA, Visser MD, Wright SJ, Rutten G, Muller-Landau H (2014). Negative density-dependence of seed dispersal and seedling recruitment in a Neotropical palm. Ecology Letters. 17:1111–1120. DOI

Pacala, S.W. (1997). Dynamics of plant communities. In: Plant Ecology. pp. 532-555.
