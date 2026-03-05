The file ABhaploid.R contains a modified AlphaBeta model that is modified for haploidy, the file BOOTmodel_Haploid.R contains bootstrap functions for the haploid epimutation model.

The neutral haploid model is the function:
```
ABneutralHaploid(pedigree.data = my.pedigree, p0uu = my.unmethylated.prop, eqp my.unmethylated.prop, eqp.weight = 0.001, Nstarts = 100, out.dir = "outputdir", out.name = "myresults")
```
