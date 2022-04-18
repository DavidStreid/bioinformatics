from asyncore import write
import sys, statistics 

TRUNCATED_FORMAT_ENTRIES = []
SKIPPED_FORMAT_ENTRIES = []

def process_variant_entry(format_col, sample_list):
  entries = []

  format_fields = format_col.split(':')
  for sample in sample_list:
    sample_vals = sample.split(':')
    if len(sample_vals) != len(format_fields):
      if len(sample_vals) < len(format_fields):
        format_fields = format_fields[0:len(sample_vals)]
        TRUNCATED_FORMAT_ENTRIES.append([format_fields, sample_vals])
      else:
        SKIPPED_FORMAT_ENTRIES.append([format_fields, sample_vals])
        continue

    entries.append(dict(zip(format_fields, sample_vals)))
    
  return entries

def get_alt_af(variant_entry):
  if 'AD' not in variant_entry:
    return None

  ad = variant_entry['AD']
  ad_vals = ad.split(',')

  ref = int(ad_vals[0])
  alt = int(ad_vals[1])

  if (ref + alt) == 0:
    return 0

  return alt / float(ref + alt)


def process_gq(variant_entry):
  if 'GQ' not in variant_entry:
    return None
  try:
    return float(variant_entry['GQ'])
  except ValueError:
    return None


def summarize(title, nums):
  num_ad = len(nums)
  max_ad = max(nums)
  min_ad = min(nums)
  mean_ad = statistics.mean(nums)
  std_ad = statistics.stdev(nums)

  print(title)
  print(f"\tNum={num_ad}")
  print(f"\tAverage={mean_ad}")
  print(f"\tSTD={std_ad}")
  print(f"\tMax={max_ad}")
  print(f"\tMin={min_ad}")

def evaluate_vcf(vcf_file):
  vcf_entries = []
  ct = 0
  with open(vcf_file, 'r') as in_f:
    for line in in_f:
      if line[0] == '#':
        continue

      ct += 1
      if ct % 1000 == 0:
        print(f"\t{ct}")
      stripped = line.strip()
      if stripped == '':
        continue
      
      cols = stripped.split('\t')
      
      chrom = cols[0]
      pos = int(cols[1])
      id = cols[2]
      ref = cols[3]
      alt = cols[4]
      qual = cols[5]
      filter = cols[6]
      info = cols[7]
      format = cols[8]
      samples = cols[9:]

      variant_entries = process_variant_entry(format, samples)
      if len(variant_entries) > 0:
        vcf_entries.append(variant_entries)


  # Allele-Depth: If this is not around 0, 0.5, or 1 it is suspect
  gq_list = []
  ad_list = []
  for variant_entries in vcf_entries:
    for variant_entry in variant_entries:
      gq = process_gq(variant_entry)
      ad = get_alt_af(variant_entry)
      
      if gq is not None:
        gq_list.append(gq)

      if ad is not None:
        ad_list.append(ad)

  # Look at any suspect alt-allele frequencies
  summarize("AD Summary [All]", ad_list)

  suspect_ad_min_filter = 0
  suspect_ad_max_filter = 0.5
  title = f"AD Summary [{suspect_ad_min_filter} < ad < {suspect_ad_max_filter}]"
  suspect_ads = [ ad for ad in ad_list if (ad > suspect_ad_min_filter and ad < suspect_ad_max_filter) ]

  summarize(title, suspect_ads)

  summarize("GQ [All]", gq_list)


def write_warn_files(fname, format_entries):
  num_entries = len(format_entries)
  if len(format_entries) == 0:
    return

  print(f"Warning - {fname}: {num_entries}")
  with open(fname, "w") as o:
    for entry in format_entries:
      o.write(f"{':'.join(entry[0])}\t{':'.join(entry[1])}\n")

if __name__ == '__main__':
  in_f = sys.argv[1]

  print("Analyzing...")
  vcf = evaluate_vcf(in_f)

  write_warn_files('truncated.tsv', TRUNCATED_FORMAT_ENTRIES)
  write_warn_files('skipped.tsv', SKIPPED_FORMAT_ENTRIES)

  print("Done.")