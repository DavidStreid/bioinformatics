#/bin/bash

SOURCE_MAP=$1
LIFTOVER_F=$2

if [[ ! -f ${SOURCE_MAP} || ! -f ${LIFTOVER_F} ]]; then
  echo "Please provide a source map & liftover file"
  exit 1
fi

# 
TEMP_DIR=temp
TEMP_MAP=./${TEMP_DIR}/temp.map
TEMP_BED=./${TEMP_DIR}/temp.bed
NEW_BED=./${TEMP_DIR}/liftOver.bed
NEW_MAP=liftOver.map

mkdir -p ${TEMP_DIR}

echo "Preparing Input Map and new Map: ${NEW_MAP}"
HEADER_LINE="#"
sed -e '1,/#/d' ${SOURCE_MAP} >> ${TEMP_MAP}
sed '/#/q' ${SOURCE_MAP} >> ${NEW_MAP}

echo "Converting ${SOURCE_MAP} to temporary BED file: ${TEMP_BED}"
while IFS=$'\n' read -r -a line; do
  chr=$(echo $line | cut -d' ' -f1)
  start=$(echo $line | cut -d' ' -f2)
  end=$(expr ${start} + 1)
  rest=$(echo $line | cut -d' ' -f3-)
  printf "${chr}\t${start}\t${end}\t${rest}\n" >> ${TEMP_BED}
done < ${TEMP_MAP}

echo "Lifting over $(realpath ${TEMP_BED}) using ${LIFTOVER_F}"
liftOver ${TEMP_BED} ${LIFTOVER_F} ${NEW_BED} unMapped

echo "Converting $(realpath ${NEW_BED}) to MAP"
while IFS=$'\n' read -r -a line; do
  chr=$(echo $line | cut -d' ' -f1)
  start=$(echo $line | cut -d' ' -f2)
  rest=$(echo $line | cut -d' ' -f4-)
  printf "${chr}\t${start}\t${rest}\n" >> ${NEW_MAP}
done < ${NEW_BED}


