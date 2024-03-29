import sys
import argparse

def main():
  parser = argparse.ArgumentParser(description='Returns nucleotide & indices of sequences')
  parser.add_argument('-f', dest='fasta', help='fasta file', required=True)
  parser.add_argument('-i', dest='index', help='nucleotide index of fasta (int)', required=True)
  parser.add_argument('-r', dest='range', help='range around index to reterieve nucelotides (int)', required=False)
  parser.add_argument('-o', dest='output', help='file to write sequence to', required=False)

  args = parser.parse_args()

  fa = args.fasta
  idx = int(args.index)
  idx_range = args.range
  if idx_range is None:
    idx_range = 3
  else:
    idx_range = int(idx_range)

  print(f"file={fa}")
  print(f"idx={idx}")
  print(f"idx_range={idx_range}")

  seq = ''
  with open(fa, 'r') as in_f:
    for raw_line in in_f:
      line = raw_line.strip()
      if line[0] == ">":
        if len(seq) > 0:  
          print(f"Checking {header}")
          analyze_seq(seq, idx)
          seq = ''
        header = line
      else:
        seq += line
    if header is not None:
      print(f"Checking {header}")
    sub_seq = analyze_seq(seq, idx, idx_range)

    if args.output is not None:
      with open(args.output, 'w') as out_f:
        header = header if header is not None else '>input'
        range = f"{idx - idx_range, idx + idx_range}"
        out_header = f"{header} {range}\n"
        out_f.write(out_header)
        out_f.write(sub_seq + "\n")

  print("Done.")


def analyze_seq(seq, idx, idx_range):
  seq_len = len(seq)

  if seq_len < idx + 1:
    print("[ERROR]")
    print(f"\tsequence_length={seq_len} is less than index={idx}")
    print(f"\tlast 10 nucleotides - \"seq[-10:]\"")
    sys.exit(1)

  index_nucl = seq[idx-1]

  print("[INFO]")
  print(f"\tsequence_length={seq_len}")
  print(f"\ttarget_nucl=({idx},{index_nucl})")
  region = ""

  min_idx = idx - idx_range
  max_idx = min(idx+idx_range+1, seq_len)

  subsection_sequence = []
  for i in range(min_idx, max_idx):
    nucl = seq[i-1]
    region += f" ({i},{nucl})"
    subsection_sequence.append(nucl)

  print(f"\tsequence=(1,{seq[0]})...{region} ...({seq_len},{seq[seq_len-1]})")
  return ''.join(subsection_sequence)

if __name__ == '__main__':
  main()
