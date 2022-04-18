
## Run
```
$ VCF_FILE=...
$ python3 evaluate_vcf.py ${VCF_FILE}
Analyzing...
AD Summary [All]
        Num=51
        Average=0.03362277341541204
        STD=0.046136712823273106
        Max=0.125
        Min=0.0
AD Summary [0 < ad < 0.5]
        Num=19
        Average=0.09025060232557969
        STD=0.02276278129820824
        Max=0.125
        Min=0.037037037037037035
GQ [All]
        Num=51
        Average=33.98039215686274
        STD=23.75835869421828
        Max=99
        Min=0
Done.
```

## Notes
### GQ
* Reasonable GQ values should have an average around 70-99 w/ a std less than 30

### AD
* Not sure, but when defining suspect AD scores w/ `[0 < ad < 0.5]`, an `AVG ~ 0.35-0.4` & `STD ~ 0.01` was an indicator of  mostly valid entries