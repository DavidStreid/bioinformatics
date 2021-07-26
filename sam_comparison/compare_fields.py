#!/usr/bin/env python

import sys
import os
from difflib import SequenceMatcher

USAGE="python compare_fields.py ${CSV1} ${CSV2}"
R1='R1'
R2='R2'

def fail(err_msg = None):
    """ Exits Application and logs error message """
    print("Usage: %s" % USAGE)
    if(err_msg):
        print("ERROR: " + err_msg)
    sys.exit(1)

def write_results_to_file(file_name, results_list):
    """ Writes @contents to @file_name
    :param file_name:
    :param contents:
    :return:
    """

    results_str = "\n".join(results_list)
    merge_commands_file = open(file_name, "a")
    merge_commands_file.truncate(0)  # Delete any old data
    merge_commands_file.write(results_str)
    merge_commands_file.close()

def is_r1(flag):
    # BIT           DESCRIPTION
    # 64    0x40    the first segment in the template
    # 128   0x80    the last segment in the template
    shift_7 = int(flag) >> 6
    if shift_7 % 2 == 1:
        return True
    return False


def create_results_dic(f, k, dic):
    lines = open(f)
    for line in lines:
        vals = line.split()
        qname = vals[0]
        flag = vals[1]
        if not flag.isdigit():
            print("Flag: %s is not a number" % flag)
            continue
        flag = int(flag)
        score = vals[2]

        if qname not in dic:
            dic[qname] = {}
        if k not in dic[qname]:
            dic[qname][k] = {}

        if flag in dic[qname][k] and score != dic[qname][k][flag]:
            print("Overriding %s => %s for %s in %s (FLAG: %d)" % (dic[qname][k][flag], score, qname, k, flag))
        dic[qname][k][flag] = score

def get_val(qname,f1,f2,files_dic,k1,k2):
    if f1 in files_dic[k1] and f2 in files_dic[k2]:
        v1 = files_dic[k1][f1]
        v2 = files_dic[k2][f2]

        if v1.isdigit() and v2.isdigit():
            diff = int(v1) - int(v2)
            return "%s,%d,%d,%s,%s,%d" % (qname,f1,f2,v1,v2,diff)
        else:
            print("Qname %s for %d & %d had non-numberic values: %s, %s" % (qname, f1, f2, v1, v2))

def get_results(k1, k2, dic):
    results = []
    for qname, files in dic.items():
        fnames = files.keys()
        if len(fnames) != 2 or k1 not in fnames or k2 not in fnames:
            print("Invalid number of files for %s: %s" % (qname, ",".join(fnames)))
            continue

        k1_flags = set(files[k1].keys())
        k2_flags = set(files[k2].keys())

        shared = k1_flags.union(k2_flags)
        for flag in list(shared):
            line = get_val(qname,flag,flag,files,k1,k2)
            if line:
                results.append(line)

        k1_uniq = list(k1_flags.difference(k2_flags))
        k2_uniq = list(k2_flags.difference(k1_flags))

        k1_r1_flags = [ f for f in k1_uniq if is_r1(f) ]
        k1_r2_flags = [ f for f in k1_uniq if not is_r1(f) ]

        k2_r1_flags = [ f for f in k2_uniq if is_r1(f) ]
        k2_r2_flags = [ f for f in k2_uniq if not is_r1(f) ]

        while len(k1_r1_flags) > 0 and len(k2_r1_flags) > 0:
            f1 = k1_r1_flags.pop()
            f2 = k2_r1_flags.pop()
            print("Adding missing flags: %d, %d" % (f1, f2))
            line = get_val(qname,f1,f2,files,k1,k2)
            if line:
                results.append(line)
        if len(k1_r1_flags) > 0 or len(k2_r1_flags) > 0:
            print("Did not add all from %s. F1: %s, F2: %s" % (qname, ",".join([ str(f) for f in k1_r1_flags]), ",".join([ str(f) for f in k2_r1_flags])))

        while len(k1_r2_flags) > 0 and len(k2_r2_flags) > 0:
            f1 = k1_r2_flags.pop()
            f2 = k2_r2_flags.pop()
            print("Adding missing flags: %d, %d" % (f1, f2))
            line = get_val(qname,f1,f2,files,k1,k2)
            if line:
                results.append(line)
        if len(k1_r2_flags) > 0 or len(k2_r2_flags) > 0:
            print("Did not add all from %s. F1: %s, F2: %s" % (qname, ",".join([ str(f) for f in k1_r2_flags]), ",".join([ str(f) for f in k2_r2_flags])))
    return results


def main():
    if len(sys.argv) != 3:
        fail("Please specify 2 CSV files")

    inputs = sys.argv[1:]

    for input in inputs:
        if not os.path.isfile(input):
            fail("%s is not a valid file" % input)

    f1 = inputs[0]
    f2 = inputs[1]

    k1 = os.path.basename(f1)
    k2 = os.path.basename(f2)

    results_dic = {}
    create_results_dic(f1, k1, results_dic)
    create_results_dic(f2, k2, results_dic)

    results = get_results(k1, k2, results_dic)
    results_file = "results___%s__%s" % (k1, k2)

    write_results_to_file(results_file, results)

if __name__ == '__main__':
    main()
