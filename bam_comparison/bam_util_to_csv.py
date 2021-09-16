#!/usr/bin/env python

import sys
import os
USAGE="python bam_util_to_csv.py /PATH/TO/bam_util_output"
class Entry:
    read = None
    v1 = None
    v2 = None

    def __init__(self, read):
        self.read = read

    def set_v1(self, v1):
        if self.v1:
            print("Read %s has at least two reads for v1: %s" % (read, self.v1))
            sys.exit(1)
        if not v1.isdigit():
            print("Read %s v1 not numeric: %s" % (read, v1))
            sys.exit(1)
        self.v1 = int(v1)

    def set_v2(self, v2):
        if self.v2:
            print("Read %s has at least two reads for v2: %s" % (read, self.v2))
            sys.exit(1)
        if not v2.isdigit():
            print("Read %s v2 not numeric: %s" % (read, v2))
            sys.exit(1)
        self.v2 = int(v2)
    def is_complete(self):
        return self.v2 and self.v1
    def to_string(self):
        v1_minus_v2 = self.v1 - self.v2
        v2_minus_v1 = self.v2 - self.v1
        return "%s,%d,%d,%d,%d" % (self.read, self.v1, self.v2, v1_minus_v2, v2_minus_v1)

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
    entries = ["template,v1,v2,v1-v2,v2-v1"]

    total_missing_reads = 0
    lines = open(bam_util_output, "r")
    first = True
    entry = Entry('HEADER')
    for line in lines:
        lead_char = line[0]
        if lead_char != "<" and lead_char != ">":
            if entry.is_complete():
                entries.append(entry.to_string())
            elif first:
                first = False
            else:
                total_missing_reads += 1
                print("Read: %s was not complete - r1: %s, r2: %s" % (entry.read, entry.v1, entry.v2))
            entry = Entry(line.strip())
        #  If the record is from the first file (--in1), it begins with a '<'. If the record is from the 2nd file (--in2), it begins with a '>'.
        else:
            score = line.split()[-1]
            if lead_char == "<":
                entry.set_v1(score)
            elif lead_char == ">":
                entry.set_v2(score)
            else:
                print("This shouldn't happen...: %s" % line)
                sys.exit(1)

    if entry.is_complete():
        entries.append(entry.to_string())
    else:
        total_missing_reads += 1
        print("Read: %s was not complete - r1: %s, r2: %s" % (entry.read, entry.v1, entry.v2))

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

    output_file = "bam_differences.csv"
    bam_util_output = inputs[0]

    print("INPUT=%s OUTPUT=%s" % (bam_util_output, output_file))
    entries = parse_entries(bam_util_output)
    print("WRITING %s" % output_file)
    write_file(output_file, entries)

if __name__ == '__main__':
    main()

