import sys

sf_tsv=sys.argv[1]
bed_base=sf_tsv.split('/')[-2]
out_bed=f'{bed_base}___interesting.bed'

print(f'reading: {sf_tsv}')
print(f'writing: {out_bed}')

with open(sf_tsv, 'r') as f, open(out_bed, 'w') as out:
  headers = f.readline().strip('\n').split('\t')
  for line in f:
    variant = dict(zip(headers, line.strip('\n').split('\t')))
    start = int(variant['variant_start_position__c'])
    chrom = variant['chromosome__c']
    if chrom == 'chrM':
      continue
    out.write(f'{chrom}\t{max(start - 100, 0)}\t{start + 100}\n')


