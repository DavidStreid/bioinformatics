#!/bin/bash

tsv=$1
bam=$2

if [[ ! -f ${tsv} || ! -f ${bam} ]]; then
  echo "Provide valid .tsv & .bam - tsv=${tsv} bam=${bam}"
  exit 1
fi

work_dir="$(basename ${bam})"
echo "Writing to ${work_dir}"
if [[ ! -d ${work_dir} ]]; then
  mkdir ${work_dir}
else
  printf "\t...already exists\n"
fi

cd ${work_dir}

python3 get_bed_intersection.py ${tsv}
bed=$(ls *.bed)

out_sam="partial.sam"
echo "Extracting ${bed} from ${bam}"
echo "Writing: $(realpath ${out_sam})"
sambamba view ${bam} -h -F "" -L ${bed} -o ${out_sam} -t 8

echo "Processing ${out_sam}..."
printf "\tAS\n"
as_out="as_scores.txt"
echo "AS" > ${as_out}
grep -oP "(?<=\t)AS[:a-zA-z0-9]+" ${out_sam} | cut -d':' -f3 >> ${as_out}

printf "\tMAPQ\n"
mapq_out="mapq_scores.txt"
echo "mapq" > ${mapq_out}
grep -oP "(?<=\t)MQ[:a-zA-z0-9]+" ${out_sam} | cut -d':' -f3 >> ${mapq_out}

printf "\tCIGAR\n"
cigar_out="cigar_out.txt"
echo "cigar" > ${cigar_out}
grep -Ev "^@" ${out_sam} | cut -f6 >> ${cigar_out}

printf "\t\tdouble-sided soft-clipped\n"
# filter the cigar scores out
ds_sc_out="double_stranded_soft_clipped.txt"
grep -oE "^[0-9]+S.*S$" ${cigar_out} > ${ds_sc_out}

printf "\t\tsingle-sided soft-clipped\n"
ss_sc_out="single_stranded_soft_clipped.txt"
grep -E "^[0-9]+S" ${cigar_out} | grep -Ev "S$" > ${ss_sc_out}
grep -E "S$" ${cigar_out} | grep -Ev "^[0-9]+S" >> ${ss_sc_out}

printf "\tClipping\n"
printf "\t\tTotal Reads: $(wc -l ${cigar_out} | grep -oE '[0-9]+')\n"
printf "\t\tsingle_soft_clipped: $(wc -l ${ss_sc_out} | grep -oE '[0-9]+')\n"
printf "\t\tdouble_soft_clipped: $(wc -l ${ds_sc_out} | grep -oE '[0-9]+')\n"

cd -
