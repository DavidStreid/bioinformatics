import sys

def main():
  if len(sys.argv) != 3:
    print("Need fasta & index")
    sys.exit(1)

  fa = sys.argv[1]
  idx = int(sys.argv[2])

  print(f"file={fa}")
  print(f"idx={idx}")

  seq = ''
  with open(fa, 'r') as in_f:
    header = in_f.readline().strip()
    if ">" in header:
      pass
    else:
      seq += header
    for line in in_f:
      seq += line.strip()
  
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
  for i in range(idx-3, idx+4):
    region += f" ({i},{seq[i-1]})"
  print(f"\tsequence=(1,{seq[0]})...{region} ...({seq_len},{seq[seq_len-1]})")
  

if __name__ == '__main__':
  main()
