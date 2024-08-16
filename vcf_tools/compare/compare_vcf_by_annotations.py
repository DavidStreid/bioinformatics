''' Compares VCF annotations as well as positions
TODO
 - Mismatch in samples (e.g. #CHROM POS ALT...S1 & #CHROM POS ALT...S2 S3)
'''
import sys, os
CHROM_ORDER = ['chrM', 'chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY']
INFO_ORDER_NOT_IMPORTANT = []

KEYWORD_MISSING = 'MISSING:MISSING'

WARNINGS = set([])

def run():
  f1 = sys.argv[1]
  f2 = sys.argv[2]
  process(f1, f2)


def process(fname_1, fname_2):
  print(f'vcf_1={fname_1}')
  print(f'vcf_1={fname_2}')
  diffs = []
  with open(fname_1, 'r') as f1, open(fname_2, 'r') as f2:
    get_past_headers(f1)
    get_past_headers(f2)

    v1 = get_nxt_var(f1)
    v2 = get_nxt_var(f2)
    while v1.valid and v2.valid:
      if are_same_pos(v1, v2):
        # Get all the variants at the same position (in case they weren't ordered)
        v1_list, v2_list = [v1], [v2]
        nxt_v1 = get_nxt_var(f1)
        nxt_v2 = get_nxt_var(f2)
        while are_same_pos(nxt_v1, v2) or are_same_pos(v1, nxt_v2):
          print(f"MULTIPLE_VARIANTS\t{v1.chrom}\t{v1.pos}")
          if are_same_pos(nxt_v1, v2):
            v1_list.append(nxt_v1)
            nxt_v1 = get_nxt_var(f1)
          elif are_same_pos(v1, nxt_v2):
            v2_list.append(nxt_v2)
            nxt_v2 = get_nxt_var(f2)
        
        # Check if any variants are the same and only compare variants that aren't
        # TODO - if for some reason there are duplicate variants in a single file, won't catch that
        same = []
        for _v1 in v1_list:
          for _v2 in v2_list:
            if same_variant(_v1, _v2):
              same.append((_v1, _v2))
        update_v1_list = [v for v in v1_list if v not in [p[0] for p in same]]
        update_v2_list = [v for v in v2_list if v not in [p[1] for p in same]]
        for _v1, _v2 in zip(update_v1_list, update_v2_list):
          diffs.extend(compare(_v1, _v2))

        v1_ct = len(update_v1_list)
        v2_ct = len(update_v2_list)
        if v1_ct > v2_ct:
          for unmatched_v1 in update_v1_list[v2_ct:]:
            diffs.extend([unmatched_v1.to_string(), KEYWORD_MISSING, unmatched_v1.to_string(), '-'])
        if v1_ct < v2_ct:
          for unmatched_v2 in update_v2_list[v1_ct:]:
            diffs.extend([unmatched_v2.to_string(), KEYWORD_MISSING, '-', unmatched_v2.to_string()])

        # Updated v1 & v2 to the next variant found
        v1 = nxt_v1
        v2 = nxt_v2
      elif is_later_variant(v1, v2):
        v1 = get_nxt_var(f1)
        diffs.extend([v1.to_string(), KEYWORD_MISSING, v1.to_string(), '-'])
      else:
        v2 = get_nxt_var(f2)
        diffs.extend([v2.to_string(), KEYWORD_MISSING, '-', v2.to_string()])
  
  if len(diffs) == 0:
    print('VCF files are the same!')
  else:
    b1 = os.path.basename(fname_1).replace('.vcf', '')
    b2 = os.path.basename(fname_2).replace('.vcf', '')
    output_name = f'diffs_{b1}.tsv' if b1 == b2 else f'diffs_{b1}_{b2}.tsv'
    process_diffs(diffs, output_name)
    # print('\n'.join(['\t'.join(vals) for vals in diffs]))


def process_diffs(diff_list, out_fname):
  diff_map = {
    v: {} for v in {d[0] for d in diff_list}
  }
  d_type_set = set([])
  d_subtype_set = set([])
  for diff in diff_list:
    variant = diff[0]
    d_subtype = diff[1]
    d_type, _ = d_subtype.split(':')
    v1 = diff[2]
    v2 = diff[3]
    d_type_set.add(d_type)
    d_subtype_set.add(diff[1])
    if d_type not in diff_map[variant]:
      diff_map[variant][d_type] = {}
    diff_map[variant][d_type][d_subtype] = [v1, v2]
  print(f'TOTAL_VARIANTS={len(diff_map)}')
  print(f'TOTAL_DIFFS={len(diff_list)}')
  print(f'\tTYPES={len(d_type_set)}\t[{", ".join(list(d_type_set)[:5])}...]')
  print(f'\tSUB_TYPES={len(d_subtype_set)}')
  
  print(f'Writing diffs to {out_fname}')
  output_order = sorted(list(d_subtype_set))
  headers = ['variant'] + output_order
  with open(out_fname, 'w') as out_f:
    out_f.write('\t'.join(headers) + '\n')
    for variant, v_map in diff_map.items():
      line_1, line_2 = [variant], [variant]
      for d_type, d_type_map in v_map.items():
        for st in output_order:
          line_1.append(d_type_map.get(st, ['-', '-'])[0])
          line_2.append(d_type_map.get(st, ['-', '-'])[1])
      out_f.write('\t'.join(line_1) + '\n')
      out_f.write('\t'.join(line_2) + '\n')



def compare(v1, v2):
  diffs = []
  if v1.chrom != v2.chrom:
    diffs.append(['chr', v1.chrom, v2.chrom])
  if v1.pos != v2.pos:
    diffs.append(['pos', v1.pos, v2.pos])
  if v1.id != v2.id:
    diffs.append(['id', v1.id, v2.id])
  if v1.ref != v2.ref:
    diffs.append(['ref', v1.ref, v2.ref])
  if v1.alt != v2.alt:
    diffs.append(['alt', v1.alt, v2.alt])
  if v1.qual != v2.qual:
    diffs.append(['qual', v1.qual, v2.qual])
  if v1.filter != v2.filter:
    diffs.append(['id', v1.filter, v2.filter])

  i1, i2 = get_info_dict(v1.info), get_info_dict(v2.info)
  diffs.extend(compare_dicts(i1, i2, 'info', INFO_ORDER_NOT_IMPORTANT))

  diffs.extend(compare_sets(v1.format, v2.format, 'format'))

  v1_smp_dicts = [get_fmt_smp_dict(v1.format, s) for s in v1.samples]
  v2_smp_dicts = [get_fmt_smp_dict(v2.format, s) for s in v2.samples]

  v1_ct, v2_ct = len(v1_smp_dicts), len(v2_smp_dicts)
  # TODO - FMT -> SAMPLE comparisons
  if v1_ct == v2_ct:
    if len(v1_smp_dicts) == 1:
      for (v1_gt_dict, v2_gt_dict) in zip(v1_smp_dicts, v2_smp_dicts):
        diffs.extend(compare_dicts(v1_gt_dict, v2_gt_dict, 'sample', []))
    else:
      WARNINGS.add('Cannot compare sample columns - can only do one sample column at this time')  
  else:
    WARNINGS.add('Cannot compare sample columns - there are different numbers of samples')
  return [[v1.to_string()] + [str(v) for v in vals] for vals in diffs]


def compare_dicts(i1, i2, label, order_not_important):
  diffs = []
  is1, is2 = set(i1.keys()), set(i2.keys())
  for k in is1.difference(is2):
    diffs.append([f'{label}:{k}', i1[k], None])
  for k in is2.difference(is1):
    diffs.append([f'{label}:{k}', None, i2[k]])
  for k in is1.intersection(is2):
    if i1[k] != i2[k]:
      if k in order_not_important:
        if set(i1[k]) == set(i2[k]):
          continue
      diffs.append([f'{label}:{k}', i1[k], i2[k]])
  return diffs


def compare_sets(l1, l2, label):
  diffs = []
  s1, s2 = set(l1), set(l2)
  for k in s1.difference(s2):
    diffs.append([f'{label}:{k}', k, None])
  for k in s2.difference(s1):
    diffs.append([f'{label}:{k}', None, k])
  return diffs


def same_variant(v1, v2):
  return len(compare(v1, v2)) == 0


def parse(raw_line):
  return raw_line.strip().split('\t')


def get_nxt_var(f):
  v = Variant()
  v.process(f.readline())
  return v


def are_same_pos(v1, v2):
  if not (v1.valid and v2.valid):
    return False
  if v1.chrom == v2.chrom:
    return v1.pos == v2.pos
  return False


def is_later_variant(v1, v2):
  if CHROM_ORDER.index(v1.chrom) > CHROM_ORDER.index(v2.chrom):
    return True
  elif CHROM_ORDER.index(v1.chrom) < CHROM_ORDER.index(v2.chrom):
    return False
  # Same chromosome
  if v1.pos > v2.pos:
    return True
  return False
  

def get_info_dict(info):
  info_pairs = info.split(';')
  info_dict = {}
  for x in info_pairs:
    if '=' in x:
      [k, v] = x.split('=')
      info_dict[k] = v
    else:
      # Parse flag as True
      info_dict[x] = True
  return info_dict


def get_fmt_smp_dict(fmt, smp):
  return dict(zip(fmt.split(':'), smp.split(':')))


def get_past_headers(f):
  line = f.readline()
  while not line.startswith('#CHROM'):
    line = f.readline()
  samples = parse(line)[9:]
  print(f'SAMPLES={samples}')
  return f


class Variant:
  def __init__(self):
    self.chrom = None
    self.pos = None
    self.id = None
    self.ref = None
    self.alt = None
    self.qual = None
    self.filter = None
    self.info = None
    self.format = None
    self.samples = None


  def process(self, raw_line):
    self.valid = True
    if raw_line in ['', '\n']:
      self.valid = False
      return
    vals = parse(raw_line)
    if len(vals) < 10:
      self.valid = False
      return
    self.chrom = vals[0]
    self.pos = int(vals[1])
    self.id = vals[2]
    self.ref = vals[3]
    self.alt = vals[4]
    self.qual = vals[5]
    self.filter = vals[6]
    self.info = vals[7]
    self.format = vals[8]
    self.samples = vals[9:]
  
  def to_string(self):
    return f'{self.chrom}#{self.pos}#{self.ref}#{self.alt}'


if __name__ == '__main__':
  run()