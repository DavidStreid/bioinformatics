#!/usr/bin/env python

import sys
import os
import pysam

USAGE="python bam_util_to_csv.py /PATH/TO/bam_util_output"

ERRORS = []

FW_FLAG = 'fw'
RC_FLAG = 'rc'
DUP_FW_FLAG = 'dup_fw'
DUP_RC_FLAG = 'dup_rc'

V1 = 'v1'
V2 = 'v2'

CONTROL_BAM = 'CONTROL_BAM'
TARGET_BAM = 'TARGET_BAM'

total_missing_reads = 0

def combine_flag_vals(dic, flag_list, val_type):
    flag_vals = []
    for f in flag_list:
         for val in dic[f][val_type]:
            entries = [ [f,v] for v in dic[f][val_type] ]
            flag_vals.extend(entries)
    return flag_vals

class Entry:
    read = None
    flag_to_val_dic = {}

    def __init__(self, read):
        self.read = read
        self.flag_to_val_dic = {}

    def add_flag_paths(self, flag):
        if flag not in self.flag_to_val_dic:
             self.flag_to_val_dic[flag] = {}
        val_dic = self.flag_to_val_dic[flag]
        if V1 not in val_dic:
            val_dic[V1] = []
        if V2 not in val_dic:
            val_dic[V2] = []

    def add_v1(self, flag, v1):
        int_v1 = int(v1)
        self.add_flag_paths(flag)
        self.flag_to_val_dic[flag][V1].append(int_v1)

    def add_v2(self, flag, v2):
        int_v2 = int(v2)
        self.add_flag_paths(flag)
        self.flag_to_val_dic[flag][V2].append(int_v2)

    def return_flag_val_pairs(self):
        global total_missing_reads

        # Note - this permanently modifies the flag_to_val_dic field
        pairs = []

        paired_up_flags = []
        # Add pairs

        for flag, val_dic in self.flag_to_val_dic.items():
            v1_vals = val_dic[V1]
            v2_vals = val_dic[V2]

            pairing_dic = {}
            while len(v1_vals) > 0:
                v1_v = v1_vals.pop()
                if v1_v in pairing_dic:
                    pairing_dic[v1_v] += 1
                else:
                    pairing_dic[v1_v] = 1

            unpaired_v2s = []
            while len(v2_vals) > 0:
                v2_v = v2_vals.pop()
                if v2_v in pairing_dic:
                    pairing_dic[v2_v] -= 1
                    # print("ADDING: {} {} {} {}".format(self.read, flag, v2_v, v2_v))
                    pairs.append([ flag, flag, v2_v, v2_v ])        # Same value
                    if pairing_dic[v2_v] == 0:          # Remove if all values have been paired up
                        del pairing_dic[v2_v]
                else:
                    unpaired_v2s.append(v2_v)

            # Any unpaired values we pair up together
            unpaired_v1s = list(pairing_dic.keys())
            unpaired_v1s.sort()
            unpaired_v2s.sort()
            while len(unpaired_v1s) > 0 and len(unpaired_v2s) > 0:
                pairs.append([flag, flag, unpaired_v1s.pop(), unpaired_v2s.pop()])

            # All values for flag are paired
            if len(unpaired_v1s) == 0 and len(unpaired_v2s) == 0:
                paired_up_flags.append(flag)
            # Unpaired values for flag are paired
            else:
                val_dic[V1] = unpaired_v1s
                val_dic[V2] = unpaired_v2s

        # Remove paired-up flags
        for flag in paired_up_flags:
            del self.flag_to_val_dic[flag]

        # Add remaining pairs from flags that have only v1 or v2 values
        remaining_flags = self.flag_to_val_dic.keys()
        rc_v1_flags = [ f for f in remaining_flags if is_flag_reverse_complemented(f) and len(self.flag_to_val_dic[f][V1]) > 0 ]
        rc_v2_flags = [ f for f in remaining_flags if is_flag_reverse_complemented(f) and len(self.flag_to_val_dic[f][V2]) > 0 ]
        fw_v1_flags = [ f for f in remaining_flags if not is_flag_reverse_complemented(f) and len(self.flag_to_val_dic[f][V1]) > 0 ]
        fw_v2_flags = [ f for f in remaining_flags if not is_flag_reverse_complemented(f) and len(self.flag_to_val_dic[f][V2]) > 0 ]

        rc_v1_vals = combine_flag_vals(self.flag_to_val_dic, rc_v1_flags, V1)
        rc_v2_vals = combine_flag_vals(self.flag_to_val_dic, rc_v2_flags, V2)
        fw_v1_vals = combine_flag_vals(self.flag_to_val_dic, fw_v1_flags, V1)
        fw_v2_vals = combine_flag_vals(self.flag_to_val_dic, fw_v2_flags, V2)

        while len(rc_v1_vals) > 0 and len(rc_v2_vals) > 0:
            rc_v1_entry = rc_v1_vals.pop()
            rc_v2_entry = rc_v2_vals.pop()
            pairs.append([rc_v1_entry[0], rc_v2_entry[0], rc_v1_entry[1], rc_v2_entry[1]])
        while len(fw_v1_vals) > 0 and len(fw_v2_vals) > 0:
            fw_v1_entry = fw_v1_vals.pop()
            fw_v2_entry = fw_v2_vals.pop()
            pairs.append([fw_v1_entry[0], fw_v2_entry[0], fw_v1_entry[1], fw_v2_entry[1]])

        if len(rc_v1_vals) > 0 or len(rc_v2_vals) > 0:
            total_missing_reads += (len(rc_v2_vals) + len(rc_v1_vals))
            ERRORS.append("UNPAIRED [RC]: V1=[{}] V2=[{}]".format(
                ",".join([str(f[0]) + ":" + str(f[1]) for f in rc_v1_vals]),
                ",".join([str(f[0]) + ":" + str(f[1]) for f in rc_v2_vals])
            ))
        if len(fw_v1_vals) > 0 or len(fw_v2_vals) > 0:
            total_missing_reads += (len(fw_v1_vals) + len(fw_v2_vals))
            ERRORS.append("UNPAIRED [FW]: V1=[{}] V2=[{}]".format(
                ",".join([str(f[0]) + ":" + str(f[1]) for f in fw_v1_vals]),
                ",".join([str(f[0]) + ":" + str(f[1]) for f in fw_v2_vals])
            ))
        return pairs

def is_flag_reverse_complemented(flag):
    # 5th bit indicates reverse complement
    fifth_bit = get_kth_bit(flag, 5)
    return fifth_bit == 1

def is_flag_duplicate(flag):
    # 11th bit indicates a duplicate
    eleventh_bit = get_kth_bit(flag, 11)
    return eleventh_bit == 1

def get_kth_bit(n, k):
    # https://www.geeksforgeeks.org/find-value-k-th-bit-binary-representation/
    return (n & (1 << (k - 1))) >> (k - 1)

def fail(err_msg = None):
    """ Exits Application and logs error message """
    print("Usage: %s" % USAGE)
    if(err_msg):
        print("ERROR: " + err_msg)
    sys.exit(1)

def write_file(file_name, entry_list):
    """ Writes @contents to @file_name
    :param file_name:
    :param contents:
    :return:
    """
    merge_commands_file = open(file_name, "a")
    merge_commands_file.truncate(0)  # Delete any old data
    merge_commands_file.write("\n".join(entry_list))
    merge_commands_file.close()

def get_flag_entry(flag):
    is_rc = is_flag_reverse_complemented(flag)
    is_dup = is_flag_duplicate(flag)


def extract_bam_entries(bam_file, entry_map, type):
    samfile = pysam.AlignmentFile(bam_file, "rb")
    aln_segments = samfile.fetch()
    for aln_seg in aln_segments:
        #aln_seg.is_reverse
        read_id = aln_seg.query_name
        flag = aln_seg.flag
        score = aln_seg.mapping_quality
        if read_id not in entry_map:
            entry_map[read_id] = Entry(read_id)
        entry = entry_map[read_id]

        if type == CONTROL_BAM:
            entry.add_v1(flag, score)
        elif type == TARGET_BAM:
            entry.add_v2(flag, score)

    return entry_map

def parse_entries(b1, b2):
    entry_map = {}
    entry_map = extract_bam_entries(b1, entry_map, CONTROL_BAM)
    entry_map = extract_bam_entries(b2, entry_map, TARGET_BAM)


    entries = ["template,flag,v1,v2,v1-v2,v2-v1"]
    added = []

    for read_id, entry in entry_map.items():
        fv_pairs = entry.return_flag_val_pairs()
        for fv_pair in fv_pairs:
            f1 = fv_pair[0]
            f2 = fv_pair[1]
            v1 = str(fv_pair[2])
            v2 = str(fv_pair[3])
            d1 = str(fv_pair[2] - fv_pair[3])
            d2 = str(fv_pair[3] - fv_pair[2])
            line = "{},{},{},{},{},{},{}".format(read_id,f1,f2,v1,v2,d1,d2)
            entries.append(line)

    if total_missing_reads > 0:
        print("Missing %d read(s)" % total_missing_reads)

    return entries

def main():
    if len(sys.argv) != 3:
        fail("Please specify two bam files")

    inputs = sys.argv[1:]
    for input in inputs:
        if not os.path.isfile(input):
            fail("%s is not a valid file" % input)

    b1 = inputs[0]
    b2 = inputs[1]

    basename = "{}____{}".format(b1.split("/")[-1].split(".")[0], b2.split("/")[-1].split(".")[0])
    output_file = "{}___bam_differences.csv".format(basename)

    print("%s=%s\n%s=%s\nOUTPUT=%s" % (CONTROL_BAM, b1, TARGET_BAM, b2, output_file))
    entries = parse_entries(b1, b2)

    print("WRITING %s" % output_file)
    write_file(output_file, entries)
    print("ERRORS")
    print("\n".join(ERRORS))

if __name__ == '__main__':
    main()
