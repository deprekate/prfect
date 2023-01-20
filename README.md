# prfect

PRFect is a tool to predict programmed ribosomal frameshifting in eukaryotic, prokaryotic, and viral genomes

It takes as input the genome and its annotated CoDing Sequence (CDS) as a GenBank file. It 
searches through the file looking for 8 different slippery site motifs associated with
backwards (-1) frameshifts and two motifs associated with forward (+1) frameshifts. When
a motif is encountered, various cellular properties and contributing factors are assed and
a prediction is made whether the site is indeed involved in programmed ribosomal frameshifting


An example genome for SARS-Cov2 is provided in the test folder, and only a single command is
needed to run PRFect:
```
python3 prfect.py test/covid19.gbk
```
and the output are any genes predicted to utilize programmed ribosomal frameshifting
(any already PRF genes that use the join keyword, are split into their two parts and
also predicted and labeled with the /label=1 feature tag)

The SARS-Cov2 genome contains 12 genes the first of which happens to be a PRF gene and 
is denoted as such through the `join` keyword.  PRFect correctly predicts that the two
gene pieces are indeed one PRF gene (after they are split into two separate genes as
mentioned above).
```
ribo frameshift detected in develop/prfect/test/covid19.gbk

     CDS             join(266..13468,13468..21555)
                     /ribosomal_slippage=-1
                     /slippery_sequence=tttaaacggg
                     /motif=is_threethree
                     /label=1
                     /product="ORF1ab polyprotein","ORF1ab polyprotein"

```
