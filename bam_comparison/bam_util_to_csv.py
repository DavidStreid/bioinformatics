#!/usr/bin/env python

import sys
import os
USAGE="python bam_util_to_csv.py /PATH/TO/bam_util_output"

ERRORS = []

class Entry:
    read = None
    v1 = None
    v2 = None
    flag = None

    def __init__(self, read, flag):
        self.read = read
        self.flag = flag

    def set_v1(self, v1):
        if not v1.isdigit():
            print("Read %s v1 not numeric: %s" % (self.read, v1))
            sys.exit(1)
        int_v1 = int(v1)
        if self.v1 and int_v1 != self.v1:
            err = "Read %s has more than two reads for v2 with different scores: first - %s, update - %s" % (self.read, self.v1, int_v1)
            ERRORS.append(err)
        self.v1 = int_v1

    def set_v2(self, v2):
        if not v2.isdigit():
            print("Read %s v2 not numeric: %s" % (self.read, v2))
            sys.exit(1)
        int_v2 = int(v2)
        if self.v2 and int_v2 != self.v2:
            err = "Read %s has more than two reads for v2 with different scores: first - %s, update - %s" % (self.read, self.v2, int_v2)
            ERRORS.append(err)
        self.v2 = int_v2

    def is_complete(self):
        return self.v2 != None and self.v1 != None
    def to_string(self):
        v1_minus_v2 = self.v1 - self.v2
        v2_minus_v1 = self.v2 - self.v1
        return "%s,%s,%d,%d,%d,%d" % (self.read, self.flag, self.v1, self.v2, v1_minus_v2, v2_minus_v1)

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

def parse_entries(bam_util_output):
    lines = open(bam_util_output, "r")
    curr_read_id=None
    entry_map = {}
    for line in lines:
        lead_char = line[0]

        if lead_char != "<" and lead_char != ">":
            # READ_ID LINE - first line of any entry
            curr_read_id = line.strip()
        else:
            # Populate score entry for curr_read_id & flag
            score = line.split()[-1]
            flag = line.split()[1]

            # Get flag_to_readId dictionary for read_id, or populate w/ empty dic
            if curr_read_id not in entry_map:
                entry_map[curr_read_id] = {}
            read_entry = entry_map[curr_read_id]

            # Get entry for flag of read_id, or populate w/ Entry
            if flag not in read_entry:
                entry = Entry(curr_read_id, flag)
                read_entry[flag] = entry
            else:
                entry = read_entry[flag]

            # Populate score for entry
            if lead_char == "<":
                entry.set_v1(score)
            elif lead_char == ">":
                entry.set_v2(score)
            else:
                print("This shouldn't happen...: %s" % line)
                sys.exit(1)

    entries = ["template,flag,v1,v2,v1-v2,v2-v1"]
    total_missing_reads = 0
    for flag_to_rIds, entry in entry_map.items():
        for flag, entry in flag_to_rIds.items():
            if entry.is_complete():
                entries.append(entry.to_string())
            else:
                total_missing_reads += 1
                err = "Read: %s was not complete - r1: %s, r2: %s" % (entry.read, entry.v1, entry.v2)
                ERRORS.append(err)

    if total_missing_reads > 0:
        print("Missing %d read(s)" % total_missing_reads)

    return entries

def main():
    if len(sys.argv) != 2:
        fail("Please specify the bam_util_output")

    inputs = sys.argv[1:]
    for input in inputs:
        if not os.path.isfile(input):
            fail("%s is not a valid file" % input)

    bam_util_output = inputs[0]

    basename = bam_util_output.split(".")[0]

    output_file = "{}___bam_differences.csv".format(basename)

    print("INPUT=%s OUTPUT=%s" % (bam_util_output, output_file))
    entries = parse_entries(bam_util_output)
    print("WRITING %s" % output_file)
    write_file(output_file, entries)
    print("ERRORS")
    print("\n".join(ERRORS))

if __name__ == '__main__':
    main()


