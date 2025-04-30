
# `get_fasta_start_stop.py`

```
$ python get_fasta_start_stop.py
Parsing >NM_006904.7...

FORWARD
START   11
Found 27 stop codons, printing first 3
CDS/CDS Complement      11..12397       # CORRECT ONE - see https://www.ncbi.nlm.nih.gov/nuccore/NM_006904.7
CDS/CDS Complement      11..12439
CDS/CDS Complement      11..12475

REVERSE
START   166
Found 184 stop codons, printing first 3
CDS/CDS Complement      166..195
CDS/CDS Complement      166..378
CDS/CDS Complement      166..447
```

# `get_seq_idx.py`
Description: Returns the nucloetide at the position provided
```
$ python3 get_seq_idx.py -f simple.fa -i 6
file=simple.fa
idx=6
[INFO]
        sequence_length=11
        target_nucl=(6,T)
        sequence=(1,A)... (3,A) (4,C) (5,A) (6,T) (7,A) (8,G) (9,A) ...(11,A)
$ cat simple.fa 
>simple
AGACATAGACA
```

Return the nucleotides around transcript for BRCA1 gene, specifically pathogenic variant [NM_007294.4(BRCA1):c.5558A>G (p.Tyr1853Cys)](https://www.ncbi.nlm.nih.gov/clinvar/variation/55627/?new_evidence=true)
```
$ python3 get_seq_idx.py -f NM_007294.4.fa -i 5558 -o NM_007294.4.cut.fa
file=NM_007294.4.fa
idx=5558
idx_range=3
Checking >NM_007294.4 Homo sapiens BRCA1 DNA repair associated (BRCA1), transcript variant 1, mRNA
[INFO]
     sequence_length=7088
     target_nucl=(5558,G)
     sequence=(1,G)... (5555,C) (5556,T) (5557,G) (5558,G) (5559,A) (5560,C) (5561,A) ...(7088,A)
Done.

$ cat NM_007294.4.cut.fa
>NM_007294.4 Homo sapiens BRCA1 DNA repair associated (BRCA1), transcript variant 1, mRNA (5555, 5561)
CTGGACA
```
