#!/bin/bash


basedir=$1

if [ ! ${basedir} ]; then
	basedir="."
fi

if [ ! -d ${basedir} ]; then
	echo "Erroooou!! Dir not found"
fi

echo "      Using basedir: ${basedir}/"

current_dir=`pwd`

##################      Lists and variables     ######################


list2=(
"GEM1"
"GEM2"
"GEM3"
"GEM4"
"GEM5"
"GEM6"
"GEM7"
"GEM8"
"GEM9"
"GEM10"
"GEM11"
"GEM12"
"GEM13"
"GEM14"
"GEM15"
"GEM16"
"GEM17"
"GEM18"
"GEM19"
"GEM20"
"GEM21"
"GEM22"
"GEM23"
"GEM24"
"GEM25"
)



#####################  cellranger count  #####################



mkdir -p ${basedir}/cellranger
cd cellranger/

for sc in ${list2[@]}
do
	mkdir -p cellranger/${sc}
	cd cellranger
	cellranger count --localmem 120 --localcores 15 --id=${sc} --fastqs=/path/to/raw_data/${sc} --transcriptome=/path/to/ref_data/refdata-gex-GRCh38-2020-A
	cd ${current_dir}
done

cd ${current_dir}

exit
