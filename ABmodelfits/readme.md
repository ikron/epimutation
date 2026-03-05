The file ABhaploid.R contains a modified AlphaBeta model that is modified for haploidy, the file BOOTmodel_Haploid.R contains bootstrap functions for the haploid epimutation model.

The neutral haploid model can be called with the function:
```
my.neutral.output <- ABneutralHaploid(pedigree.data = my.pedigree, p0uu = my.unmethylated.prop, eqp my.unmethylated.prop, eqp.weight = 0.001, Nstarts = 100, out.dir = "outputdir", out.name = "myresults")
```
Testing whether model of neutral epimutation accumulation fits better than a null model can be done:
```
FtestRSSHaploid(pedigree.select = "path2mymodel/myneutral.Rdata", pedigree.null = "path2mymodel/mynull.Rdata")
```
The haploid model has fewer parameters than the default neutral model in AlphaBeta, so a modified FtestRSS functions needs to be used.

To perform bootstrap
```
my.boot.results <- BOOTmodel_Haploid(model.fit = my.neutral.model.output, Nboot = 1000, out.dir = "outputdir", out.name = "my_bootresults")
```
