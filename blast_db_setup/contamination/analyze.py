import os
import argparse

from pandas import qcut

EXPECTED_BLASTNAMES = set([]) # Add list of blastnames expected to be in results, e.g. "primates" if seaching human reads
OUTPUT_TSV = 'output.tsv'
IDENTITY_TSV = 'query_identities.tsv'

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


def output_relevant_results(query_dic):
  total = 0
  num_filtered_queries = len(query_dic)
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

  
  print(f"num_filtered_queries={num_filtered_queries}")
  print(f"\ttotal={total}")

def get_query_dic(blast_results_tsv):
  query_dic = {}
  queries_with_expected_top_hit = set()

  total_queries = 0
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

      if qaccver in queries_with_expected_top_hit:
        continue

      if qaccver not in query_dic:
        total_queries += 1  
        if sblastnames in EXPECTED_BLASTNAMES:
          # If the first hit for the qaccver is expected, skip all other blast results for the query
          queries_with_expected_top_hit.add(qaccver)
          # print(f"{qaccver} - {sblastnames}")
        highest_value = evalue

      if sblastnames not in EXPECTED_BLASTNAMES and evalue < 1e-4 and (evalue / highest_value) < 10:
        # Ignore anything expected
        # Ignore anything that's an order of magnitude lower than the first BLAST hit
        if qaccver not in query_dic:
          query_dic[qaccver] = []
        query_dic[qaccver].append(br)
      else:
        pass
        # print(f"{qaccver} - {sblastnames} - {evalue}")

  print(f"Total Queries = {float(total_queries)}")
  print(f"\tQueries with expected blast hits = {len(queries_with_expected_top_hit)}")
  
  return query_dic

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='Identify likelihood of contamination')
  parser.add_argument('-f', dest='blast_results_tsv', help='TSV of blast results', required=True)

  args = parser.parse_args()

  blast_results_tsv = args.blast_results_tsv

  query_dic = get_query_dic(blast_results_tsv)
  output_relevant_results(query_dic)
