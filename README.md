# prfect
PRFect is a tool to predict programmed ribosomal frameshifting in eukaryotic, prokaryotic, and viral genomes


The published manuscript is available at:
https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-024-05701-0


<br>

PRFect takes as input the genome and its annotated CoDing Sequences (CDS) as a GenBank file.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;*&nbsp; *If you only have a fasta file we recommend our brand new gene caller [Genotate](https://github.com/deprekate/genotate) that is* <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; *the only gene caller that can call gene fragments*

PRFect searches through a GenBank file looking for 8 different slippery site motifs associated with
backwards (-1) frameshifts and two motifs associated with forward (+1) frameshifts. When
a motif is encountered, various cellular properties and factors are assessed and
a prediction is made whether the site is involved in programmed ribosomal frameshifting.

<br>

To install:
```
python3 -m pip install prfect
```
To run:
```
prfect.py input.gbk
```
An example genome for SARS-Cov2 is provided in the test folder. The SARS-Cov2 genome contains 12 genes the first of which happens to be a PRF gene and is denoted as such through the use of the `join` keyword.  Any genes already present that use the `join` keyword are split into their two parts and subsequently predicted anew and then tagged with the /label=1 feature tag to indicate a TruePositive.  When the genome is run through PRFect the known PRF gene is correctly predicted to utilize programmed ribosomal frameshifting.

```
$ prfect.py test/covid19.gbk 

     CDS             join(266..13468,13468..21555)
                     /ribosomal_slippage
                     /direction=-1
                     /motif=is_threethree
                     /slippery_sequence=tttaaac
                     /label=1
                     /locus=NC_045512
                     /product="ORF1ab polyprotein"
                     /product="ORF1ab polyprotein"

```

Another example is bacteriophage lambda, which has the *geneG* and *geneGT* tail assembly chaperone gene that is known to frameshift.  The current genbank annotation file (NC_001416) does not have the gene properly denoted with the `join` keyword and so both pieces are in two separate CDS features.  When the genome is run through PRFect the gene is correctly identified as being a single PRF gene with the /label=0 to indicate that it is an UnknownPositive.

```
$ prfect.py test/lambda.gbk

     CDS             join(9711..10115,10115..10549)
                     /ribosomal_slippage
                     /direction=-1
                     /motif=is_threethree
                     /bases=gggaaag
                     /label=0
                     /locus=NC_001416
                     /product="minor tail protein G"
                     /product="tail assembly protein T"
```


You can show all the slippery sites that PRFect checked to make sure it evaluated a given site and to see if there were any near hits.
Using the `--dump` flag will show the calculated cellular properites at each potential slippery site:
```
$ prfect.py test/lambda.gbk --dump | head
LOCUS      SLIPSITE   LOC  LABEL  N  DIR RBS1 RBS2  A0     A1     LF50    HK50    LF100   HK100  PRED  PROB  MOTIF
NC_001416  gcaaaacgc  4278   0  159   1   13   1.8  0.015  0.025  -0.24   -0.236  -0.523  -0.306   0    1.0  three
NC_001416  ggaaagtgt  10115  0   18  -1    2     0  0.004  0.024  -0.313  -0.287  -0.668  -0.404  -1   0.88  threethree  
NC_001416  gcgaaagca  31034  0   30   1    2   1.0  0.029  0.032  -0.282  -0.243  -0.477  -0.326   0    1.0  three
NC_001416  tggaaacgc  33370  0   72   1    1     0  0.015  0.028  -0.124  -0.118  -0.482  -0.36    0    1.0  three
NC_001416  cgtaaatta  33388  0   90   1    0     0  0.009  0.012  -0.15   -0.138  -0.291  -0.237   0    1.0  three
NC_001416  gcagggtgg  33442  0  144   1    0     0  0.017  0.021  -0.092  -0.039  -0.388  -0.274   0    1.0  three
NC_001416  gaaaaggag  42081  0   42  -1    0     0  0.027  0.013  -0.246  -0.149  -0.176  -0.105   0    1.0  twofour
NC_001416  aaaaccttc  42206  0   66  -1    0     0  0.015  0.014  -0.403  -0.266  -0.367  -0.249   0    1.0  fivetwo
NC_001416  cgaaaaaat  43240  0    6   1    2     0  0.019  0.023  -0.513  -0.245  -0.395  -0.294   0   0.98  four
```




The columns are:
```
LOCUS     id of the sequence
SLIPSITE  bases of the slippery site
LOC       location within the bases of the slippery site
LABEL     whether the slippery site is already annotated: 0 not a joined gene, 1 a joined gene, -1 a joined gene but is >10bp away 
N         distance of the slippery site from the in-frame stop codon
DIR       direction of the shift
RBS1      Prodigal like ribosomal binding site interference score
RBS2      RAST like ribosomal binding site interference score
A0        frequency of the A-site codon usage in all genes
A1        frequency of the +1 A-site codon usage in all genes
LF50      normalized LinearFold minimum free energy calculation of the downstream 50bp window
LF100     normalized LinearFold minimum free energy calculation of the downstream 100bp window
HK50      normalized HotKnots minimum free energy calculation of the downstream 50bp window
HK100     normalized HotKnots minimum free energy calculation of the downstream 100bp window
PRED      type of shift predicted by PRFect to occur: -1 backwards, 0 no shift, +1 forwards
PROB      how sure PRFect was for the predicted (PRED) type
MOTIF     slippery sequence motif
```


You can even use the flag `-s` to scale the MFE calculations to account for extreme GCcontent/temp/salinity:
```
$ prfect.py test/lambda.gbk -s 1.5 --dump | head -n 2
LOCUS      SLIPSITE   LOC  LABEL  N  DIR RBS1 RBS2  A0     A1     LF50    HK50    LF100   HK100  PRED  PROB  MOTIF
NC_001416  gcaaaacgc  4278   0  159   1   13   1.8  0.015  0.025  -0.36   -0.354  -0.785  -0.459   0    1.0  three
NC_001416  ggaaagtgt  10115  0   18  -1    2     0  0.004  0.024  -0.47   -0.431  -1.002  -0.606  -1  0.999  threethree  
```
you will notice that the MFE values were scaled by 50% when compared to the above dump, which also caused the trained model to be more confident in the backward -1 PREDiction at LOCation 10115




