# prfect


```
python3 prfect.py test/covid19.gbk
```
and the output are genes that are predicted to utilize programmed ribosomal frameshifing
(any already PRF genes that use the join keyword, are split into their two parts and also predicted
```
ribo frameshift detected in develop/prfect/test/covid19.gbk

     CDS             join(266..13468,13468..21555)
                     /ribosomal_slippage=-1
                     /slippery_sequence=tttaaacggg
                     /motif=is_threethree
                     /label=1
                     /product="ORF1ab polyprotein","ORF1ab polyprotein"

```
