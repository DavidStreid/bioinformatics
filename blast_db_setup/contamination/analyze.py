import argparse
import collections
import matplotlib.pyplot as plt

EXPECTED_BLASTNAMES = set()   # BLAST `sblastnames` that are ignored, e.g. if looking for contamination in human saliva, add "primates"

OUTPUT_TSV = 'blast_result_analysis.tsv'
IDENTITY_TSV = 'query_identities.tsv'
IDENTITY_SUMMARY_TSV = 'identity_summary.tsv'

MAX_E_VALUE = 0.1                     # E-Values greater than this are filtered out, NOTE - "E=1" means 1 result expected by chance
MAX_MAGNITUDE_DIFFERENCE_ALLOWED = 3  # E-Values more than 10^x the best hit are filtered out, e.g. "3" will filter out anything with E-value > (1000 * lowest E) 

class blast_result:
  @staticmethod
  def get_header():
    return ['qaccver', 'pident', 'evalue', 'magnitude_of_best_result', 'proportion_valid_blast_results', 'total_results', 'bitscore', 'scomnames', 'sblastnames', 'sskingdoms']
  
  @staticmethod
  def get_id_header():
    return ['qaccver', 'scomnames', 'sblastnames']

  def __init__(self, qaccver, pident, evalue, bitscore, scomnames, sblastnames, sskingdoms):
    self.qaccver = qaccver
    self.pident = pident
    self.evalue = evalue
    self.bitscore = bitscore
    self.scomnames = scomnames
    self.sblastnames = sblastnames
    self.sskingdoms = sskingdoms
    self.magnitude = None             # Order of magnitude higher this result is from the next BETTER result
    self.next_magnitude = 0           # Order of magnitude lower this result is from the next WORSE result
    self.prop_valid_blast_results = 0
    self.total_valid_blast_results = 0
    self.total_blast_results = 0
    self.best_evalue_result_sblastnames_proportion = 0

  def clone(self):
    br = blast_result(self.qaccver, self.pident, self.evalue, self.bitscore, self.scomnames, self.sblastnames, self.sskingdoms)
    br.next_magnitude = self.next_magnitude
    br.prop_valid_blast_results = self.prop_valid_blast_results
    br.best_evalue_result_sblastnames_proportion = self.best_evalue_result_sblastnames_proportion

    return br

  def set_best_evalue_result_sblastnames_proportion(self, best_evalue_result_sblastnames_proportion):
    self.best_evalue_result_sblastnames_proportion = best_evalue_result_sblastnames_proportion

  def set_total_blast_results(self, total_blast_results):
    self.total_blast_results = total_blast_results

  def set_total_valid_blast_results(self, total_valid_blast_results):
    self.total_valid_blast_results = total_valid_blast_results

  def set_prop_valid_blast_results(self, prop_valid_blast_results):
    self.prop_valid_blast_results = prop_valid_blast_results

  def set_next_magnitude(self, blast_result):
    self.next_magnitude = ("%.2f" % blast_result.magnitude)

  def set_magnitude(self, next_best_evalue):
    self.magnitude = self.evalue / float(next_best_evalue)

  def to_identity_string(self):
    return '\t'.join([self.qaccver, self.scomnames, self.sblastnames])

  def to_string(self):
    return '\t'.join([self.qaccver, str(self.pident), str(self.evalue), str(self.next_magnitude), str(self.prop_valid_blast_results), str(self.total_blast_results), str(self.bitscore), self.scomnames, self.sblastnames, self.sskingdoms])

  def is_the_best_match(self):
    statistically_significant = self.evalue < 1e-50
    evolutionarily_close = self.pident > 0.9


def output_relevant_results(query_dic):
  with open(OUTPUT_TSV, 'w') as out:
    header = '\t'.join(blast_result.get_header())
    out.write(f"{header}\n")
    for highest_hit in query_dic.values():
      out.write(f'{highest_hit.to_string()}\n')

  print("[OUT] BLAST RESULT SUMMARY")
  print(f"\tIDS & STATISTICS={OUTPUT_TSV}")
  print(f"\tIDS ONLY={IDENTITY_TSV}")

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
    last_evalue = None

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
        last_evalue = br.evalue
        filtered_reads_dic[qaccver] = []

      br.set_magnitude(last_evalue)       # Gets order of magnitude worse that the result is from the last value
      filtered_reads_dic[qaccver].append(br)

      last_evalue = br.evalue
  
  print(f"\tTotal Reads = {len(filtered_reads_dic)}")

  return filtered_reads_dic

def get_best_result_sblastnames_proportion(blast_result_list, best_result_blastname):
  blast_names = [ blast_result.sblastnames for blast_result in blast_result_list ]
  freq_dic = collections.Counter(blast_names)

  return freq_dic[best_result_blastname]

def is_valid_blast_result(blast_result, lowest_evalue_for_read_id):
  is_valid_evalue = blast_result.evalue < MAX_E_VALUE
  similar_evalue = (blast_result.evalue / lowest_evalue_for_read_id) < (10 * MAX_MAGNITUDE_DIFFERENCE_ALLOWED)

  return is_valid_evalue and similar_evalue

def aggregate_blast_results(query_dic):
  aggregated_query_dic = {}       # Clone of aggregated query dictionary with only the best blast result
  
  invalid_blast_ids = set()

  total_blast_results = 0
  for qaccver, blast_result_list in query_dic.items():
    total_blast_results_for_qaccver = len(blast_result_list)
    total_blast_results += total_blast_results_for_qaccver


    blast_result_list.sort(key=lambda br: br.evalue)
    best_evalue_result = blast_result_list[0]

    best_evalue_for_read_id = best_evalue_result.evalue

    if len(blast_result_list) > 1:
      best_evalue_result.set_next_magnitude(blast_result_list[1])
    
    num_valid_blast_results = len([ br for br in blast_result_list if is_valid_blast_result(br, best_evalue_for_read_id)])
    if num_valid_blast_results == 0:
      invalid_blast_ids.add(qaccver)
      continue

    best_evalue_result.set_total_blast_results(total_blast_results_for_qaccver)
    best_evalue_result.set_total_valid_blast_results(num_valid_blast_results)
    best_evalue_result.set_prop_valid_blast_results(num_valid_blast_results / float(total_blast_results_for_qaccver))

    identification = best_evalue_result.sblastnames

    best_evalue_result_sblastnames_proportion = get_best_result_sblastnames_proportion(blast_result_list, identification)
    best_evalue_result.set_best_evalue_result_sblastnames_proportion(best_evalue_result_sblastnames_proportion)

    aggregated_query_dic[qaccver] = best_evalue_result.clone()

  print("AGGREGATION")
  print(f"\ttotal_blast_results={total_blast_results}")
  num_removed_reads = len(invalid_blast_ids)
  if num_removed_reads > 0:
    print(f"\t[WARNING] Filtered Out Reads Reads = {len(num_removed_reads)}")
  print(f"\tTotal Reads with Valid Blast Results = {len(aggregated_query_dic)}")

  return aggregated_query_dic


class Graph:
  def __init__(self, x_values, y_values, title):
    self.x_values = x_values
    self.y_values = y_values
    self.title = title

  def graph(self):
    assert len(self.x_values) == len(self.y_values)

    # https://matplotlib.org/stable/gallery/pie_and_polar_charts/pie_features.html
    graph_title = self.title.upper() 
    fname = self.title.replace(" ", "_") + '.pdf'

    print("[OUT] GRAPHING")
    print(f"\tfile={fname} (n={len(self.x_values)})")

    plt.title(graph_title)
    _, ax1 = plt.subplots()
    ax1.pie(self.y_values, labels=self.x_values, autopct='%1.1f%%', startangle=90)
    ax1.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle.
    plt.savefig(fname)
    plt.close()


def log_settings():
  print("SETTINGS")
  print(f"\tMAX_E_VALUE={MAX_E_VALUE}")
  print(f"\tMAX_MAGNITUDE_DIFFERENCE_ALLOWED={MAX_MAGNITUDE_DIFFERENCE_ALLOWED}")


def parse_identification_summary(aggregated_query_dic):
  read_ids = {}

  for blast_result in aggregated_query_dic.values():
    identification = blast_result.sblastnames
    if identification in read_ids:
      read_ids[identification] += 1
    else:
      read_ids[identification] = 1

  print("IDENTIFICATIONS")
  identifications = []
  identification_counts = []

  highest_count_identification = ''
  highest_count = 0
  total_ct = 0

  for identification, id_count in read_ids.items():
    identifications.append(identification)
    identification_counts.append(id_count)

    total_ct += id_count

    if id_count > highest_count:
      highest_count = id_count
      highest_count_identification = identification

  proportion_of_highest_result_ct = "%.2f" % (highest_count / float(len(aggregated_query_dic)))
  print(f"MOST_IDENTIFIED={highest_count_identification} ({proportion_of_highest_result_ct})")

  id_tuples = list(zip(identifications, identification_counts))
  id_tuples.sort(key=lambda it: it[1], reverse=True)

  with open(IDENTITY_TSV, 'w') as identity_out, open(IDENTITY_SUMMARY_TSV, 'w') as identity_summary_out:
    header = '\t'.join(blast_result.get_id_header())
    identity_out.write(f"{header}\n")
    for highest_hit in aggregated_query_dic.values():
      identity_out.write(f'{highest_hit.to_identity_string()}\n')

    identity_summary_out.write("ID\tCount\tProportion\n")
    for it in id_tuples:
      ct = it[1]
      proportion = ct/float(total_ct)
      line = f"{it[0]}\t{ct}\t{proportion}\n"
      identity_summary_out.write(line)
      

  return Graph(identifications, identification_counts, 'blast identifications per read id')


if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='Identify likelihood of contamination')
  parser.add_argument('-f', dest='blast_results_tsv', help='TSV of blast results', required=True)
  args = parser.parse_args()

  blast_results_tsv = args.blast_results_tsv
  log_settings()

  query_dic = get_query_dic(blast_results_tsv)

  aggregated_blast_results = aggregate_blast_results(query_dic)

  blast_graph = parse_identification_summary(aggregated_blast_results)
  blast_graph.graph()

  output_relevant_results(aggregated_blast_results)