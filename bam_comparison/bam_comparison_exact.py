'''
$ python bam_comparison_exact.py sample.1.bam sample.2.bam
Processing chrM
	chrM	0
	chrM	250000
	chrM	500000
	chrM	750000
	chrM	1000000
	chrM	1250000
READ_COMPARISON	CHROM=chrM	SHARED=637935	R1_ONLY=0	R2_ONLY=0
Processing chr1
	chr1	1500000
    ...
$ cut -f2,3 bam_comparison.chrM.columns_only.tsv  | grep -v "-" | sort | uniq -c
     12 1	flag
 597917 1	query_qualities
    249 1	tags__MD
    243 1	tags__NM
    206 2	flag,query_qualities
   8981 2	query_qualities
   8981 2	tags__MD
   8903 2	tags__NM
'''

import sys
import pysam

CORE_FIELDS = ['query_name', 'flag', 'reference_id', 'reference_start',
          'mapping_quality', 'cigarstring', 'next_reference_id',
          'next_reference_start', 'template_length', 'query_sequence', 'query_qualities']


def normalize_tags(tags):
  """Convert list of (tag, value) to a sorted dict for easy comparison."""
  return dict(tags)


def get_read_info(read1):
  """Compare core fields and tags of two reads."""
  # Compare core fields
  core_vals = [getattr(read1, field) for field in CORE_FIELDS]
  tags1 = normalize_tags(read1.get_tags())
  return (core_vals, tags1)


def compare_reads(r1_dict, r2_dict, current_ref, out_f):
  r1_names = set(r1_dict.keys())
  r2_names = set(r2_dict.keys())

  shared = r1_names.intersection(r2_names)
  r1_only = r1_names.difference(r2_names)
  r2_only = r2_names.difference(r1_names)

  print(f'READ_COMPARISON\tCHROM={current_ref}\tSHARED={len(shared)}\tR1_ONLY={len(r1_only)}\tR2_ONLY={len(r2_only)}')
  if len(r1_only) > 0 or len(r2_only) > 0:
    print(f'WARNING\tREF={current_ref}\tSHARED={len(shared)}\tR1_ONLY={len(r1_only)}\tR2_ONLY={len(r2_only)}')

  keys_only_cols = ['read_name', 'num_columns', 'columns']
  keys_only_out = f'{out_f}.{current_ref}.columns_only.tsv'
  all_diffs_cols = ['read_name', 'column', 'v1', 'v2']
  all_differences = f'{out_f}.{current_ref}.values.tsv'
  different, same = 0, 0
  with open(keys_only_out, 'w') as c_f, open(all_differences, 'w') as a_f:
    c_f.write('\t'.join(keys_only_cols) + '\n')
    a_f.write('\t'.join(all_diffs_cols) + '\n')
    for rname in shared:
      differences = {}
      if len(r1_dict[rname]) != len(r2_dict[rname]):
        print(f'WARNING\t{rname}\t{len(r1_dict[rname])} != {len(r2_dict[rname])}')
      for r1_read, r2_read in zip(r1_dict[rname], r2_dict[rname]):
        cv_1, tags1 = r1_read[0], r1_read[1]
        cv_2, tags2 = r2_read[0], r2_read[1]

        if tags1 != tags2:
          # Show only differing tags
          tag_diff = {k: (tags1.get(k), tags2.get(k)) for k in set(tags1.keys()).union(tags2.keys()) if tags1.get(k) != tags2.get(k)}
          differences['tags'] = tag_diff

        for idx, field in enumerate(CORE_FIELDS):
          val1, val2 = cv_1[idx], cv_2[idx]
          if val1 != val2:
            differences[field] = (val1, val2)
      n_diffs = len(differences)
      if n_diffs > 0:
        different += 1
        for t in differences.get('tags', []):
          label = f'tags__{t}'
          c_f.write('\t'.join([rname, str(n_diffs), label]) + '\n')
          a_f.write('\t'.join([rname, label] + [str(v) for v in differences['tags'][t]]) + '\n')
        if 'tags' in differences:
          del differences['tags']
        c_f.write('\t'.join([rname, str(n_diffs), ','.join(list(differences.keys()))]) + '\n')
        for d, v in differences.items():
          if len(v) != 2:
            print(f'WARNING\t{rname}\t{v}')
          a_f.write('\t'.join([rname, d] + [str(val) for val in v]) + '\n')
      else:
        same += 1
        c_f.write('\t'.join([rname, '0', '-']) + '\n')


def compare_bams(bam_path1, bam_path2, out_f):
  bam1 = pysam.AlignmentFile(bam_path1, "rb")
  bam2 = pysam.AlignmentFile(bam_path2, "rb")

  current_ref = None # reference_id
  r1_dict, r2_dict = {}, {}

  for idx, (read1, read2) in enumerate(zip(bam1.fetch(until_eof=True), bam2.fetch(until_eof=True))):
    if read1.reference_id != read2.reference_id:
      raise ValueError(f'Out of order\tB1={read1.query_name}\tB2={read2.query_name}')
    ref = bam1.get_reference_name(read1.reference_id)
    if ref != current_ref:
      if current_ref is not None:
        compare_reads(r1_dict, r2_dict, current_ref, out_f)
      r1_dict, r2_dict = {}, {}
      current_ref = ref
      print(f'Processing {current_ref}')

    if read1.query_name not in r1_dict:
      r1_dict[read1.query_name] = []
    r1_dict[read1.query_name].append(get_read_info(read1))
    if read2.query_name not in r2_dict:
      r2_dict[read2.query_name] = []
    r2_dict[read2.query_name].append(get_read_info(read2))

    if idx % 250_000 == 0:
      print(f'\t{current_ref}\t{idx}')

  bam1.close()
  bam2.close()

# Example usage:
if __name__ == "__main__":
  bam_file_1 = sys.argv[1]
  bam_file_2 = sys.argv[2]
  out_f = 'bam_comparison'
  differences = compare_bams(bam_file_1, bam_file_2, out_f)