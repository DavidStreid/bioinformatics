# Haplotype Map Converter
## INPUTS
`INPUT_MAP`: Input Haplotype Map that needs to be converted
`CHAIN_FILE`: Chain file translating alignments of one genome to another. See http://hgdownload.cse.ucsc.edu/goldenpath/${GENOME}}/liftOver/ for any `GENOME`

## Process
* Converts `INPUT_MAP` to a temporary BED file
* Uses liftOver (See liftOver in [UCSC Utilities](http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/)) to conver the BED file to the target genome using `CHAIN_FILE`
* Converts temporary BED file back to the final `.map` file

## Run
```
CHAIN_FILE=hg38ToHg19.over.chain
./convert_haplotype_map.sh ${INPUT_MAP} ${CHAIN_FILE}
```
