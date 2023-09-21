# Creates mock FASTQ files
* `./get_fastas.sh` - Returns FASTA of name/id in local BLAST DB
* `./make_fastqs.sh` - Returns FASTQ of input fasta file

## DEPENDENCIES
* Assumes BLAST is setup locally with the `DB` assigned in `get_fastas.sh`

## RUN
```
tax_id="11137"
name="Humancoronavirus 229E"
fasta=$(./get_fastas.sh ${tax_id} "${name}")
fastq=$(./make_fastqs.sh ${fasta})

printf "${tax_id}\t${name}\t${fasta}\t${fastq}\n"
```

```
$ head -3 taxid_name.tsv
taxid	name	subtype	type
10298	Human alphaherpesvirus 1	Herpes Simplex Virus (HSV)	virus
10366	Murid betaherpesvirus 1	Cytomegalovirus (CMV)
$ ./create_blast_FAQs.sh
```


