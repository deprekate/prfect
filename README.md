# prfect

PRFect is a tool to predict programmed ribosomal frameshifting in eukaryotic, prokaryotic, and viral genomes

The manuscript is currently in review and has been preprinted on bioRxiv:

https://www.biorxiv.org/content/10.1101/2023.04.10.536325v1
<br>

It takes as input the genome and its annotated CoDing Sequence (CDS) as a GenBank file. It 
searches through the file looking for 8 different slippery site motifs associated with
backwards (-1) frameshifts and two motifs associated with forward (+1) frameshifts. When
a motif is encountered, various cellular properties and contributing factors are assessed and
a prediction is made whether the site is indeed involved in programmed ribosomal frameshifting.
<br>
<br>
To install:
```
pip3 install prfect
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
