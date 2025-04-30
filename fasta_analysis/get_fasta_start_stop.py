'''
Analyzes a given transcript for start/stop positions
If an NCBI transcript, can confirm via https://www.ncbi.nlm.nih.gov/nuccore/<TRANSCRIPT_ID> and looking for "CDS" lines

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
'''

MAX_TO_PRINT = 3
stop_codons_rna = ['UAA', 'UAG', 'UGA']
stop_codons_dna = [c.replace('U', 'T') for c in stop_codons_rna]
stop_codons = set(stop_codons_rna + stop_codons_dna)


def run():
  seq = read_fa('refseq.fa')
  print('FORWARD')
  read_seq(seq)

  print('\nREVERSE')
  read_seq(reverse_complement(seq))


def read_fa(fname):
  seq = ''
  with open(fname, 'r') as in_f:
    for raw_line in in_f:
      line = raw_line.strip('\n')
      if line.startswith('>'):
        print(f"Parsing {line.split(' ')[0]}...\n")
        continue
      seq += line
  return seq


def read_seq(seq):
  for i in range(len(seq)):
    if seq[i:i+3] == 'ATG':
      start = i
      print(f'START\t{ncbi_pos(start)}')
      break
  stops = []
  i = start
  while i < len(seq):
    codon = seq[i:i+3]
    if codon in stop_codons:
      # print(f'STOP\t{ncbi_pos(i)}:{ncbi_pos(i+2)}\t{codon}')
      stops.append((i, i+2))
    i += 3

  output(start, stops, max=MAX_TO_PRINT)

  return start, stops


def reverse_complement(seq):
  complement = str.maketrans('ATCGNatcgn', 'TAGCNtagcn')
  return seq.translate(complement)[::-1]


def output(start, stops, max=3):
  print(f'Found {len(stops)} stop codons, printing first {max}')
  for i in range(max):
    stop = stops[i]
    print(f'CDS/CDS Complement\t{ncbi_pos(start)}..{ncbi_pos(stop[1])}')


def ncbi_pos(p):
  return p + 1


if __name__ == '__main__':
  run()
