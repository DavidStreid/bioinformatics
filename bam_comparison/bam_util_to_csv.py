#!/usr/bin/env python

import sys
import os
import pysam

USAGE="python bam_util_to_csv.py control.bam experimental.bam"

# Save all errors until end and then log
ERRORS = []
TOTAL_MISSING_READS = 0

# Constants to track value differences between the two BAMs
CONTROL_BAM = 'CONTROL_BAM'
TARGET_BAM = 'TARGET_BAM'

def combine_flag_vals(dic, flag_list, val_type):
    """Extends all the values across all flags of an Entry:flag_to_val_dic field

    :param Entry:flag_to_val_dic dic:   flag-to-entries dictionary
    :param int[] flag_list:             list of flags, e.g. [ 63, 1023,... ]
    :param string, val_type:            CONTROL_BAM/TARGET_BAM
    :return: int[ int[] ]               list of [flag, value] pairs
    """
    flag_vals = []
    for f in flag_list:
         for val in dic[f][val_type]:
            entries = [ [f,v] for v in dic[f][val_type] ]
            flag_vals.extend(entries)
    return flag_vals

class Entry:
    """ Tracks all the records for an input read name (SAM QNAME)

    :field read:                        QNAME of SAM record
    :field flag_to_val_dic:             flag -> { CONTROL/TARGET -> [ MAPQ_values ] }
    """
    read = None
    flag_to_val_dic = {}

    def __init__(self, read):
        self.read = read
        self.flag_to_val_dic = {}

    def add_flag_paths(self, flag):
        """ Safely adds path for a flag (populates if the path in the dictionary doesn't exist)
        """
        if flag not in self.flag_to_val_dic:
             self.flag_to_val_dic[flag] = {}
        val_dic = self.flag_to_val_dic[flag]
        if CONTROL_BAM not in val_dic:
            val_dic[CONTROL_BAM] = []
        if TARGET_BAM not in val_dic:
            val_dic[TARGET_BAM] = []

    def add_v1(self, flag, v1):
        """ Adds the CONTROL_BAM flag-MAPQ_value pair
        """
        int_v1 = int(v1)
        self.add_flag_paths(flag)
        self.flag_to_val_dic[flag][CONTROL_BAM].append(int_v1)

    def add_v2(self, flag, v2):
        """ Adds the TARGET flag-MAPQ_value pair
        """
        int_v2 = int(v2)
        self.add_flag_paths(flag)
        self.flag_to_val_dic[flag][TARGET_BAM].append(int_v2)

    def return_flag_val_pairs(self):
        global TOTAL_MISSING_READS

        # Note - this permanently modifies the flag_to_val_dic field
        pairs = []

        paired_up_flags = []
        # Add pairs

        for flag, val_dic in self.flag_to_val_dic.items():
            v1_vals = val_dic[CONTROL_BAM]
            v2_vals = val_dic[TARGET_BAM]

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
                val_dic[CONTROL_BAM] = unpaired_v1s
                val_dic[TARGET_BAM] = unpaired_v2s

        # Remove paired-up flags
        for flag in paired_up_flags:
            del self.flag_to_val_dic[flag]

        # Add remaining pairs from flags that have only v1 or v2 values
        remaining_flags = self.flag_to_val_dic.keys()
        rc_v1_flags = [ f for f in remaining_flags if is_flag_reverse_complemented(f) and len(self.flag_to_val_dic[f][CONTROL_BAM]) > 0 ]
        rc_v2_flags = [ f for f in remaining_flags if is_flag_reverse_complemented(f) and len(self.flag_to_val_dic[f][TARGET_BAM]) > 0 ]
        fw_v1_flags = [ f for f in remaining_flags if not is_flag_reverse_complemented(f) and len(self.flag_to_val_dic[f][CONTROL_BAM]) > 0 ]
        fw_v2_flags = [ f for f in remaining_flags if not is_flag_reverse_complemented(f) and len(self.flag_to_val_dic[f][TARGET_BAM]) > 0 ]

        rc_v1_vals = combine_flag_vals(self.flag_to_val_dic, rc_v1_flags, CONTROL_BAM)
        rc_v2_vals = combine_flag_vals(self.flag_to_val_dic, rc_v2_flags, TARGET_BAM)
        fw_v1_vals = combine_flag_vals(self.flag_to_val_dic, fw_v1_flags, CONTROL_BAM)
        fw_v2_vals = combine_flag_vals(self.flag_to_val_dic, fw_v2_flags, TARGET_BAM)

        while len(rc_v1_vals) > 0 and len(rc_v2_vals) > 0:
            rc_v1_entry = rc_v1_vals.pop()
            rc_v2_entry = rc_v2_vals.pop()
            pairs.append([rc_v1_entry[0], rc_v2_entry[0], rc_v1_entry[1], rc_v2_entry[1]])
        while len(fw_v1_vals) > 0 and len(fw_v2_vals) > 0:
            fw_v1_entry = fw_v1_vals.pop()
            fw_v2_entry = fw_v2_vals.pop()
            pairs.append([fw_v1_entry[0], fw_v2_entry[0], fw_v1_entry[1], fw_v2_entry[1]])

        if len(rc_v1_vals) > 0 or len(rc_v2_vals) > 0:
            TOTAL_MISSING_READS += (len(rc_v2_vals) + len(rc_v1_vals))
            t_or_c = CONTROL_BAM if len(rc_v1_vals) > 0 else TARGET_BAM
            ERRORS.append("UNPAIRED_RC,{},{},{},{}".format(
                self.read,
                t_or_c,
                ",".join([str(f[0]) + ":" + str(f[1]) for f in rc_v1_vals]),
                ",".join([str(f[0]) + ":" + str(f[1]) for f in rc_v2_vals])
            ))
        if len(fw_v1_vals) > 0 or len(fw_v2_vals) > 0:
            TOTAL_MISSING_READS += (len(fw_v1_vals) + len(fw_v2_vals))
            t_or_c = CONTROL_BAM if len(fw_v1_vals) > 0 else TARGET_BAM
            ERRORS.append("UNPAIRED_FW,{},{},{},{}".format(
                self.read,
                t_or_c,
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
        read_id = aln_seg.query_name
        flag = aln_seg.flag
        score = aln_seg.mapping_quality

        # TODO - Add refernece name. It would be interesting to see if reads are aligned to different scaffolds
        # ref = aln_seg.reference_name

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

    if TOTAL_MISSING_READS > 0:
        print("Missing %d read(s)" % TOTAL_MISSING_READS)

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
    missing_file = "{}___missing.csv".format(basename)

    print("%s=%s\n%s=%s\nOUTPUT=%s" % (CONTROL_BAM, b1, TARGET_BAM, b2, output_file))
    entries = parse_entries(b1, b2)

    print("WRITING %s" % output_file)
    write_file(output_file, entries)
    print("ERRORS")
    print("\n".join(ERRORS))
    write_file(missing_file,ERRORS)

if __name__ == '__main__':
    main()
