import argparse
import collections
import matplotlib.pyplot as plt

EXPECTED_BLASTNAMES = set()   # BLAST `sblastnames` that are ignored, e.g. if looking for contamination in human saliva, add "primates"
OUTPUT_TSV = 'output.tsv'
IDENTITY_TSV = 'query_identities.tsv'

VALID_ID_PROPORTION = 0.5

MAX_E_VALUE = 0.1                     # E-Values greater than this are filtered out, NOTE - "E=1" means 1 result expected by chance
MAX_MAGNITUDE_DIFFERENCE_ALLOWED = 3  # E-Values more than 10^x the best hit are filtered out, e.g. "3" will filter out anything with E-value > (1000 * lowest E) 

class blast_result:
  def __init__(self, qaccver, pident, evalue, bitscore, scomnames, sblastnames, sskingdoms):
    self.qaccver = qaccver
    self.pident = pident
    self.evalue = evalue
    self.bitscore = bitscore
    self.scomnames = scomnames
    self.sblastnames = sblastnames
    self.sskingdoms = sskingdoms

  def to_identity_string(self):
    return '\t'.join([self.qaccver, self.scomnames, self.sblastnames])

  def to_string(self):
    return '\t'.join([self.qaccver, str(self.pident), str(self.evalue), str(self.bitscore), self.scomnames, self.sblastnames, self.sskingdoms])

  def is_the_best_match(self):
    statistically_significant = self.evalue < 1e-50
    evolutionarily_close = self.pident > 0.9


def graph_pie(labels, sizes, fname):
  # https://matplotlib.org/stable/gallery/pie_and_polar_charts/pie_features.html
  plt.title('BLAST IDENTIFICATIONS')

  fig1, ax1 = plt.subplots()
  ax1.pie(sizes, labels=labels, autopct='%1.1f%%', startangle=90)
  ax1.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle.
  plt.savefig("%s.pdf" % fname)
  plt.close()


def output_relevant_results(query_dic):
  total = 0
  reads_with_valid_blast_results = len(query_dic)
  with open(OUTPUT_TSV, 'w') as out, open(IDENTITY_TSV, 'w') as identity_out:
    header = '\t'.join(['qaccver', 'pident', 'evalue', 'bitscore', 'scomnames', 'sblastnames', 'sskingdoms'])
    out.write(f"{header}\n")
    for qaccver, blast_result_list in query_dic.items():
      highest_hit = blast_result_list[0]
      highest_hit_evalue = highest_hit.evalue

      total_hits = len(blast_result_list)

      total += total_hits
      if highest_hit_evalue < 1e-50:
        # If the first hit is this certain, skip everything else
        out.write(f'{highest_hit.to_string()}\n')
        identity_out.write(f'{highest_hit.to_identity_string()}\n')
        continue
      else:
        for br in blast_result_list:
          out.write(f'{br.to_string()}\n')

      dic_of_commonname = {}
      list_of_sblastnames = [ br.to_identity_string() for br in blast_result_list ]
      
      for identity_string in list_of_sblastnames:
        if identity_string not in dic_of_commonname:
          dic_of_commonname[identity_string] = 1
        else:
          dic_of_commonname[identity_string] += 1
      for identity_string, ct in dic_of_commonname.items():
        if ct / float(total_hits) > 0.5:
          identity_out.write(f'{identity_string}\n')
          # print(f"{qaccver} - {sbn}")

      # print(f"{qaccver} - {list(set([ br.sblastnames for br in blast_result_list ]))}")

  print("[OUT] BLAST RESULT SUMMARY")
  print(f"\treads_with_valid_blast_results={reads_with_valid_blast_results}")
  print(f"\ttotal={total}")

def get_query_dic(blast_results_tsv):
  ''' Returns a dictionary of read IDs to their blast results
  :return { 'read_id': blast_result[], ... }
  '''
  print("[IN] BLAST RESULT SUMMARY")
  print(f"\tINPUT = {blast_results_tsv}")
  
  removed_reads = set()

  total_reads = 0
  filtered_reads_dic = {}
  with open(blast_results_tsv, 'r') as in_f:
    for line in in_f:
      # qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids sscinames scomnames sblastnames sskingdoms
      cols = line.strip('\n').split('\t')

      qaccver = cols[0]
      pident = cols[2]          # % Identity
      evalue = float(cols[10])
      bitscore = float(cols[11])
      scomnames = cols[14]
      sblastnames = cols[15]
      sskingdoms = cols[16]

      br = blast_result(qaccver, pident, evalue, bitscore, scomnames, sblastnames, sskingdoms)

      if qaccver not in filtered_reads_dic:
        total_reads += 1  
        highest_value = evalue

      low_evalue = evalue < MAX_E_VALUE                                                     # Less than one result expected by chance
      similar_evalue = (evalue / highest_value) < (10 * MAX_MAGNITUDE_DIFFERENCE_ALLOWED)   # All results are within three orders of magnitude of the best result  

      if low_evalue and similar_evalue:
        # Read ID has a valid E-Value and will be analyzed
        if qaccver not in filtered_reads_dic:
          filtered_reads_dic[qaccver] = []
        filtered_reads_dic[qaccver].append(br)
      else:
        if qaccver not in filtered_reads_dic:
          # Read ID has no valid E-Values and won't analyzed
          removed_reads.add(qaccver)
        else:
          # Read ID has valid E-Values, just not this one
          pass
  
  print(f"\tTotal Reads = {total_reads}")
  print(f"\tTotal Reads with Valid Blast Results = {len(filtered_reads_dic)}")
  num_removed_reads = len(removed_reads)
  if num_removed_reads > 0:
    print(f"\t[WARNING] Removed Reads = {len(num_removed_reads)}")

  return filtered_reads_dic

def get_proportions(blast_result_list):
  blast_names = [ blast_result.sblastnames for blast_result in blast_result_list ]
  freq_dic = collections.Counter(blast_names)
  total = len(blast_names)

  for id, id_freq in freq_dic.items():
    if (id_freq/total) > VALID_ID_PROPORTION:
      return freq_dic, id

  # Nothing met the minimum to be ID'd
  return freq_dic, None


def aggregate_blast_results(query_dic):
  AMBIGUOUS_KEY = 'ambiguous'
  read_ids = {
    AMBIGUOUS_KEY: 0
  }
  
  for qaccver, blast_result_list in query_dic.items():
    qaccver_freq_dic, id = get_proportions(blast_result_list)
    if id is not None:
      if id in read_ids:
        read_ids[id] += 1
      else:
        read_ids[id] = 1
    else:
      read_ids[AMBIGUOUS_KEY] += 1
      # TODO - what to do with this?
      # for k,v in qaccver_freq_dic.items():
      #  print(f"\t{k}: {v}")

  print("IDENTIFICATIONS")
  identifications = []
  identification_counts = []
  for identification, id_count in read_ids.items():
    identifications.append(identification)
    identification_counts.append(id_count)
    print(f"\t{identification}: {id_count}")    

  graph_pie(identifications, identification_counts, 'blast_results')


def log_settings():
  print("SETTINGS")
  print(f"\tMAX_E_VALUE={MAX_E_VALUE}")
  print(f"\tMAX_MAGNITUDE_DIFFERENCE_ALLOWED={MAX_MAGNITUDE_DIFFERENCE_ALLOWED}")

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='Identify likelihood of contamination')
  parser.add_argument('-f', dest='blast_results_tsv', help='TSV of blast results', required=True)
  args = parser.parse_args()

  blast_results_tsv = args.blast_results_tsv
  log_settings()

  query_dic = get_query_dic(blast_results_tsv)
  blast_results = aggregate_blast_results(query_dic)
  output_relevant_results(query_dic)