import sys

RUN_STR = "python3 vcf_to_bed.py ${VCF}"

def read_vcf(vcf_f, bed_f):
  with open(vcf_f, 'r') as in_f, open(bed_f, 'w') as out_f:
    for raw_line in in_f:
      line = raw_line.strip()
      if line[0] == '#':
        continue
      
      cols = line.split('\t')
      chrom = cols[0]
      start = cols[1]
      info = cols[7]

      info_fields = info.split(';')
      end, svlen = None, None
      for field in info_fields:
        [k,v] = field.split('=')
        if k == 'END':
          end = v
        if k == 'SVLEN':
          svlen = v
       
      if end is not None:
         line = f"{chrom}\t{start}\t{end}\n"
      elif svlen is not None:
         # TODO - handle deletions where SVLEN is usually a negative number
         end = int(start) + int(svlen)
         line = f"{chrom}\t{start}\t{end}\n"
      else:
         print(f"skipping CHROM={chrom} POS={start} - no svlen or send in INFO")
         continue
      
      out_f.write(line)
      

if __name__ == '__main__':
  if len(sys.argv) == 1:
    print("[ERROR] No VCF file passed")
    print(RUN_STR)
    sys.exit(1)

  print("[WARN] I only do VCFs w/ structural variant insertions")

  vcf_f = sys.argv[1]
  bed_f = vcf_f.replace(".vcf", ".bed")
  print(f"INPUT={vcf_f}")
  read_vcf(vcf_f, bed_f)
  print(f"OUTPUT={bed_f}")
  print("Done.")
