''' Compares VCF annotations as well as positions
TODO
 - Mismatch in samples
 - Difference in ordering of multiple variants at the same position
'''
import sys
CHROM_ORDER = ['chrM', 'chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY']
INFO_ORDER_NOT_IMPORTANT = []


def run():
  f1 = sys.argv[1]
  f2 = sys.argv[2]
  process(f1, f2)


def process(fname_1, fname_2):
  diffs = []
  with open(fname_1, 'r') as f1, open(fname_2, 'r') as f2:
    get_past_headers(f1)
    get_past_headers(f2)

    v1 = get_nxt_var(f1)
    v2 = get_nxt_var(f2)
    while v1.valid and v2.valid:
      if are_same_pos(v1, v2):
        diffs.extend(compare(v1, v2))
        v1 = get_nxt_var(f1)
        v2 = get_nxt_var(f2)
      elif is_later_variant(v1, v2):
        v1 = get_nxt_var(f1)
        diffs.extend([v1.to_string(), 'MISSING', v1.to_string(), '-'])
      else:
        v2 = get_nxt_var(f2)
        diffs.extend([v2.to_string(), 'MISSING', '-', v2.to_string()])
  
  if len(diffs) == 0:
    print('VCF files are the same!')
  else:
    print(f'VCF file diffs={len(diffs)}')
    print('\n'.join(['\t'.join(vals) for vals in diffs]))


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

  # TODO - FMT -> SAMPLE comparisons
  for s in v1.samples:
    s1, s2 = get_fmt_smp_dict(v1.format, s), get_fmt_smp_dict(v1.format, s)
    diffs.extend(compare_dicts(s1, s2, 'sample', []))
  
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


def parse(raw_line):
  return raw_line.strip().split('\t')


def get_nxt_var(f):
  v = Variant()
  v.process(f.readline())
  return v


def are_same_pos(v1, v2):
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