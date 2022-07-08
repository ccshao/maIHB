An R package to detect the maternal allele inheritance from next-generation sequencing data
via the hierarchical Bayesian model.

## Installation
First install JAGS. 
Then install the package
```
install.packages("remotes")
remotes::install_github("ccshao/maIHB")
```

## Examples
Run hierarchical Bayesian model on builtin data F1.TypeA
```
byesRes.A <- BernHierModel(F1.TypeA[, -1], saveName = "Family1.TypeA")
```

Summarize the results with graphs
```
summaryRes.A <- smryMCMC(byesRes.A, saveName = "Family1.TypeA")
```

## Reference
Inference of maternal allele inheritance via hierarchical Bayesian model in noninvasive prenatal diagnosis.
bioRxiv doi: http://dx.doi.org/10.1101/051995
