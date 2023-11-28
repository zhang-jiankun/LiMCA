#!/bin/bash

usage() {
	cat <<-EOF
		RNAC_for_olfactory.sh
		Usage: `basename $0` 
		Options:
		    -g: genome.fa
		    -f: R1.fq
		    -r: R2.fq
		    -t: N, threads
		    -s: GENOME ASSEMBLY
		    -v: SNP FILE 
		    -h: help information
		
		Action:
		    -A: Reads mapping
		    -P: Contacts parsing, haplotype imputation
		    -S: Reconstruction of 3D genome structures, structure filtering
		    -C: Contact analysis of Greek islands and OR genes for olfactory neurons
		    -D: Similar to -C, but 3D structure analysis for olfactory neurons

		See also:
		    GEO README (ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE121nnn/GSE121791/suppl/GSE121791%5F00README%2Emd%2Etxt)
	EOF
	exit 1;
}

GENOME=""
R1=""
R2=""
nThreads=8
PREFIX=""
OUTPUT=""
ASSEMBLY=""
PAR="/share/home/zhangjk/scRNAC/MOE/snps/mm10.par.bed"
or_con="/dshare/xielab/analysisdata/ThreeDGene/zhangjk/scRNAC/Olfactory/proc_data/RNA/ORs/olfactory_receptors.con"
or_data="/dshare/xielab/analysisdata/ThreeDGene/zhangjk/scRNAC/Olfactory/proc_data/RNA/ORs"

# extract options and their arguments into variables
while getopts "g:f:r:o:t:p:s:v:hAPCSD" options; do
	case ${options} in 
		g) GENOME=${OPTARG};;
		f) R1=${OPTARG};;
		r) R2=${OPTARG};;
		o) OUTPUT=${OPTARG};;
		t) nThreads=${OPTARG};;
		p) PREFIX=${OPTARG};;
		s) ASSEMBLY=${OPTARG};;
		v) SNP=${OPTARG};;
		A) ALIGN=1;;
		P) PARSE=1;;
		S) MODEL=1;;
		C) CONTACT=1;;
		D) SPATIAL=1;;
		h) usage;;
		*) exit 1;;
	esac
done

##### configure software path
dipc="/share/home/zhangjk/hic_tools/dip-c/scripts"
folder="/share/home/zhangjk/hic_tools/dip-c/color"
or="/share/home/zhangjk/scRNAC/MOE/or_data"
juicer="java -Xmx2g -jar /share/home/zhangjk/software/juicer_tools_1.22.01.jar"

[[ -z ${OUTPUT} ]] && { echo "OUTPUT: NONE!"; exit 1; }
[ -d ${OUTPUT} ] || { mkdir -p ${OUTPUT}; }
[ -d ${OUTPUT}/logs ] || mkdir -p ${OUTPUT}/logs

################################# 1. alignment #################################

if [ ${ALIGN-0} == 1 ];then

	echo  -e "$(date)\n"
	echo  -e "GENOME: ${GENOME}"
	echo  -e "R1: ${R1}"
	echo  -e "R2: ${R2}"
	echo  -e "OUTPUT: ${OUTPUT}"
	echo  -e "PREFIX: ${PREFIX}"
	echo  -e "---------------------------------------------------\n"

	echo -e "FASTQ stats:\n"
	
	echo -e "R1(${R1})\t${PREFIX}\t$(kseq_fastq_base ${R1})"
	echo -e "R2(${R2})\t${PREFIX}\t$(kseq_fastq_base ${R2})"
	echo -e "---------------------------------------------------\n"

	[[ -z ${GENOME} || ! -s ${GENOME} ]] && { echo "GENOME: not found!"; exit 1; }
	[[ -z ${R1} || ! -s ${R1} ]] && { echo "R1: not found!"; exit 1; }
	[[ -z ${R2} || ! -s ${R2} ]] && { echo "R2: not found!"; exit 1; }

	seqtk mergepe ${R1} ${R2} | pre-meta - | bwa mem -5SP -p ${GENOME} -t ${nThreads} - 2>${OUTPUT}/logs/${PREFIX}_mapping.log | gzip > ${OUTPUT}/${PREFIX}.sam.gz

fi

################################# 2. parse contacts #################################

if [ ${PARSE-0} == 1 ];then

	[[ -z ${SNP} || ! -s ${SNP} ]] && { echo "SNP file: not found!"; exit 1; }
	echo "SNP file: ${SNP}"
	
	## male
	hickit.js sam2seg -v ${SNP} ${OUTPUT}/${PREFIX}.sam.gz 2>${OUTPUT}/logs/${PREFIX}_sam2seg.log | hickit.js chronly - | hickit.js bedflt ${PAR} - | gzip > ${OUTPUT}/${PREFIX}.allValidSeg.gz
	hickit -i ${OUTPUT}/${PREFIX}.allValidSeg.gz -o - 2>${OUTPUT}/logs/${PREFIX}_seg2pairs.log | sed "s/chromosome/chromsize/g" | bgzip > ${OUTPUT}/${PREFIX}.pairs.gz
	zcat ${OUTPUT}/${PREFIX}.pairs.gz | awk -v FS='\t' -v OFS='\t' '{if(substr($0, 1, 1)=="#"){print}else{$6="+";$7="+";print}}' | bgzip > ${OUTPUT}/${PREFIX}.temp.pairs.gz
	hickit --dup-dist=1000 -i ${OUTPUT}/${PREFIX}.temp.pairs.gz -o - 2>${OUTPUT}/logs/${PREFIX}_seg2pairs.log2 | bgzip > ${OUTPUT}/${PREFIX}.allValidPairs.gz

	# convert to .hic format
	$juicer pre -n ${OUTPUT}/${PREFIX}.allValidPairs.gz ${OUTPUT}/${PREFIX}.hic ${ASSEMBLY} > ${OUTPUT}/logs/${PREFIX}_juicer.log
	cut -f1-2 ${GENOME}.fai | grep "^chr" > ${OUTPUT}/${ASSEMBLY}.chrom.sizes
	pairix -p pairs ${OUTPUT}/${PREFIX}.allValidPairs.gz
	pairtools stats -o ${OUTPUT}/logs/${PREFIX}_allValidPairs_stats.txt ${OUTPUT}/${PREFIX}.allValidPairs.gz
	pairsqc.py -p ${OUTPUT}/${PREFIX}.allValidPairs.gz -c ${OUTPUT}/${ASSEMBLY}.chrom.sizes -t P -O ${OUTPUT}/${PREFIX}

	[ -d ${OUTPUT}/structure ] || mkdir -p ${OUTPUT}/structure
	
	ln -fs ${OUTPUT}/${PREFIX}.allValidPairs.gz ${OUTPUT}/structure/contacts.pairs.gz
	# resolve haplotypes via imputation
	hickit -i ${OUTPUT}/structure/contacts.pairs.gz -u -o - 2>${OUTPUT}/logs/${PREFIX}_impute.log | bgzip > ${OUTPUT}/structure/${PREFIX}.impute.pairs.gz
	ln -fs ${OUTPUT}/structure/${PREFIX}.impute.pairs.gz ${OUTPUT}/structure/impute.pairs.gz
	${dipc}/hickit_pairs_to_con.sh ${OUTPUT}/structure/contacts.pairs.gz
	${dipc}/hickit_impute_pairs_to_con.sh ${OUTPUT}/structure/impute.pairs.gz
	
	# haplotype-resolved contacts
	${dipc}/con_imputed_to_juicer_pre_short.sh ${OUTPUT}/structure/impute.con.gz
fi

################################# 3. modeling #################################

if [ ${MODEL-0} == 1 ];then
	# submit jobs to slurm scheduler
	for rep in `seq 1 3`;do
		sbatch -J ${PREFIX} -o ${OUTPUT}/logs/${PREFIX}_modeling_rep${rep}.log <<-EOF
			#!/bin/bash
			#SBATCH -c 1
			#SBATCH --mem 8G
			#SBATCH --time 2-00:00:00

			hickit -s${rep} -M -i ${OUTPUT}/structure/${PREFIX}.impute.pairs.gz -Sr1m -c1 -r10m -c2 -b4m \
				-b1m -O ${OUTPUT}/structure/${PREFIX}_1Mb_${rep}.3dg \
				-b200k -O ${OUTPUT}/structure/${PREFIX}_200kb_${rep}.3dg \
				-D5 -b50k -O ${OUTPUT}/structure/${PREFIX}_50kb_${rep}.3dg \
				-D5 -b20k -O ${OUTPUT}/structure/${PREFIX}_20kb_${rep}.3dg 

			RESARR=("1Mb" "200kb" "50kb" "20kb")
			for RES in \${RESARR[@]};do
				${dipc}/hickit_3dg_to_3dg_rescale_unit.sh ${OUTPUT}/structure/${PREFIX}_\${RES}_${rep}.3dg
				dip-c clean3 -c ${OUTPUT}/structure/impute.con.gz ${OUTPUT}/structure/${PREFIX}_\${RES}_${rep}.dip-c.3dg > ${OUTPUT}/structure/${PREFIX}_\${RES}_${rep}.clean.3dg
			done
		EOF
	done

fi

################################# 4. contact map analysis #################################

if [ ${CONTACT-0} == 1 ]; then

	[ -s ${OUTPUT}/features ] || mkdir -p ${OUTPUT}/features
	[ -s ${OUTPUT}/or ] || mkdir -p ${OUTPUT}/or # stands for olfactory receptor 
	
	$juicer pre -n ${OUTPUT}/structure/impute.juicer.txt.gz ${OUTPUT}/${PREFIX}.impute.hic ${folder}/mm10.chr.hom.len
	dip-c color2 -b1000000 -H -c ${folder}/mm10.cpg.1m.txt -s ${OUTPUT}/structure/contacts.con.gz > ${OUTPUT}/features/${PREFIX}_cpg_b1m.color2

	## Enhancers and olfactory receptors analyses 

	echo -e "\n$(date) analyze OR .......... \n"
	dip-c ard -d10000000 -c ${or_con} -h100000 ${OUTPUT}/structure/contacts.con.gz > ${OUTPUT}/or/${PREFIX}_ors_d10m_h100k.ard
	dip-c ard -d10000000 -c ${or_con} -n ${OUTPUT}/structure/contacts.con.gz > ${OUTPUT}/or/${PREFIX}_ors_d10m_n.ard
	dip-c ard -d200000 -c ${or_con}  -n ${OUTPUT}/structure/contacts.con.gz > ${OUTPUT}/or/${PREFIX}_ors_d200k_n.ard

	echo -e "\n$(date) analyze all enhancers ..........\n"
	dip-c ard -d10000000 -c ${or_data}/GI_peaks.con -h100000 ${OUTPUT}/structure/contacts.con.gz > ${OUTPUT}/or/${PREFIX}_enhancers_d10m_h100k.ard
	dip-c ard -d10000000 -c ${or_data}/GI_peaks.con -n ${OUTPUT}/structure/contacts.con.gz > ${OUTPUT}/or/${PREFIX}_enhancers_d10m_n.ard
	dip-c ard -d100000 -c ${or_data}/GI_peaks.con -n ${OUTPUT}/structure/contacts.con.gz > ${OUTPUT}/or/${PREFIX}_enhancers_d100k_n.ard
	
	echo -e "\n$(date) analyze candidate enhancers ..........\n"
	dip-c ard -d10000000 -c ${or_data}/candidate_peaks.con -h100000 ${OUTPUT}/structure/contacts.con.gz > ${OUTPUT}/or/${PREFIX}_candidate_enhancers_only_d10m_h100k.ard
	dip-c ard -d10000000 -c ${or_data}/candidate_peaks.con -n ${OUTPUT}/structure/contacts.con.gz > ${OUTPUT}/or/${PREFIX}_candidate_enhancers_only_d10m_n.ard
	dip-c ard -d100000 -c ${or_data}/candidate_peaks.con -n ${OUTPUT}/structure/contacts.con.gz > ${OUTPUT}/or/${PREFIX}_candidate_enhancers_only_d100k_n.ard

	echo -e "\n$(date) analyze canonical enhancers ..........\n"
	dip-c ard -d10000000 -c /share/home/zhangjk/scRNAC/MOE/or_data/GSE121791_genomic_data.enhancers.con -h100000 ${OUTPUT}/structure/contacts.con.gz > ${OUTPUT}/or/${PREFIX}_canonical_enhancers_only_d10m_h100k.ard
	dip-c ard -d10000000 -c /share/home/zhangjk/scRNAC/MOE/or_data/GSE121791_genomic_data.enhancers.con -n ${OUTPUT}/structure/contacts.con.gz > ${OUTPUT}/or/${PREFIX}_canonical_enhancers_only_d10m_n.ard
	dip-c ard -d100000 -c /share/home/zhangjk/scRNAC/MOE/or_data/GSE121791_genomic_data.enhancers.con -n ${OUTPUT}/structure/contacts.con.gz > ${OUTPUT}/or/${PREFIX}_canonical_enhancers_only_d100k_n.ard
			
	echo -e "\n$(date) analyze subset 0 enhancers ..........\n"
	dip-c ard -d10000000 -c ${or_data}/subset_27GIs_0.con -h100000 ${OUTPUT}/structure/contacts.con.gz > ${OUTPUT}/or/${PREFIX}_subset_enhancers_0_d10m_h100k.ard
	dip-c ard -d10000000 -c ${or_data}/subset_27GIs_0.con -n ${OUTPUT}/structure/contacts.con.gz > ${OUTPUT}/or/${PREFIX}_subset_enhancers_0_d10m_n.ard
	dip-c ard -d100000 -c ${or_data}/subset_27GIs_0.con -n ${OUTPUT}/structure/contacts.con.gz > ${OUTPUT}/or/${PREFIX}_subset_enhancers_0_d100k_n.ard

	echo -e "\n$(date) analyze subset 1 enhancers ..........\n"
	dip-c ard -d10000000 -c ${or_data}/subset_27GIs_1.con -h100000 ${OUTPUT}/structure/contacts.con.gz > ${OUTPUT}/or/${PREFIX}_subset_enhancers_1_d10m_h100k.ard
	dip-c ard -d10000000 -c ${or_data}/subset_27GIs_1.con -n ${OUTPUT}/structure/contacts.con.gz > ${OUTPUT}/or/${PREFIX}_subset_enhancers_1_d10m_n.ard
	dip-c ard -d100000 -c ${or_data}/subset_27GIs_1.con -n ${OUTPUT}/structure/contacts.con.gz > ${OUTPUT}/or/${PREFIX}_subset_enhancers_1_d100k_n.ard

	echo -e "\n$(date) analyze subset 2 enhancers ..........\n"
	dip-c ard -d10000000 -c ${or_data}/subset_27GIs_2.con -h100000 ${OUTPUT}/structure/contacts.con.gz > ${OUTPUT}/or/${PREFIX}_subset_enhancers_2_d10m_h100k.ard
	dip-c ard -d10000000 -c ${or_data}/subset_27GIs_2.con -n ${OUTPUT}/structure/contacts.con.gz > ${OUTPUT}/or/${PREFIX}_subset_enhancers_2_d10m_n.ard
	dip-c ard -d100000 -c ${or_data}/subset_27GIs_2.con -n ${OUTPUT}/structure/contacts.con.gz > ${OUTPUT}/or/${PREFIX}_subset_enhancers_2_d100k_n.ard  
fi

if [ ${SPATIAL-0} == 1 ];then

	[ -s ${OUTPUT}/mmCIF ] || mkdir -p ${OUTPUT}/mmCIF
	
	dip-c color -n ${folder}/mm10.chr.txt ${OUTPUT}/conformation/${PREFIX}_20kb_1.clean.3dg | dip-c vis -c /dev/stdin ${OUTPUT}/conformation/${PREFIX}_20kb_1.clean.3dg > ${OUTPUT}/mmCIF/${PREFIX}_20kb_1.clean.n.cif 
	dip-c color -c ${folder}/mm10.cpg.20k.txt ${OUTPUT}/conformation/${PREFIX}_20kb_1.clean.3dg | dip-c vis -M -c /dev/stdin ${OUTPUT}/conformation/${PREFIX}_20kb_1.clean.3dg > ${OUTPUT}/mmCIF/${PREFIX}_cpg.cif 

	# calculate raidal positioning
	dip-c color -C ${OUTPUT}/conformation/${PREFIX}_20kb_1.clean.3dg > ${OUTPUT}/features/${PREFIX}_C.color
	dip-c color -R -c ${folder}/mm10.cpg.20k.txt ${OUTPUT}/conformation/${PREFIX}_20kb_1.clean.3dg > ${OUTPUT}/features/R_cpg.color
	dip-c color -R --min-num=1 -c ${or_data}/olfactory_receptors.binary.20k.txt ${OUTPUT}/conformation/${PREFIX}_20kb_1.clean.3dg > ${OUTPUT}/or/R_olfactory_receptors.color  # need rerun

	dip-c color -C ${output}/conformation/${sample}_1Mb_1.clean.3dg > ${output}/features/${sample}_C.1Mb.color 
	dip-c color -R -c ${folder}/mm10.cpg.1m.txt  ${output}/conformation/${sample}_1Mb_1.clean.3dg > ${output}/features/R_cpg.1Mb.color
	dip-c color -R --min-num=1 -c ${or_data}/olfactory_receptors.binary.1m.txt ${output}/conformation/${sample}_1Mb_1.clean.3dg > ${output}/or/R_olfactory_receptors.1Mb.color

	### Enhancers and olfactory receptors are analyses from 3D genomes with dip-c:
	# enhancer near each enhancer
	dip-c pd -1 ${or_data}/enhancers.hom.leg.txt ${OUTPUT}/conformation/${PREFIX}_20kb_1.clean.3dg > ${OUTPUT}/or/${PREFIX}.enhancers.hom.pd.matrix
	
	python ${dipc}/network_around.py 2.5 ${OUTPUT}/or/${PREFIX}.enhancers.hom.pd.matrix \
		${or_data}/enhancers.hom.leg.txt ${or_data}/enhancers.hom.leg.txt ${or_data}/enhancers.hom.name.txt ${or_data}/enhancers.hom.name.txt \
		> ${OUTPUT}/or/${PREFIX}.d2.5.enhancers.network.txt

	python ${dipc}/network_around.py 5 ${OUTPUT}/or/${PREFIX}.enhancers.hom.pd.matrix \
		${or_data}/enhancers.hom.leg.txt ${or_data}/enhancers.hom.leg.txt ${or_data}/enhancers.hom.name.txt ${or_data}/enhancers.hom.name.txt \
		> ${OUTPUT}/or/${PREFIX}.d5.enhancers.network.txt

	# distance_threshold = 2.5
	python ${dipc}/network.py ${OUTPUT}/or/${PREFIX}.enhancers.hom.pd.matrix \
		${or_data}/enhancers.hom.leg.txt ${or_data}/enhancers.hom.name.txt > ${OUTPUT}/or/${PREFIX}.d2.5.enhancers.network.comp.txt

	python ${dipc}/network_5.py ${OUTPUT}/or/${PREFIX}.enhancers.hom.pd.matrix \
		${or_data}/enhancers.hom.leg.txt ${or_data}/enhancers.hom.name.txt > ${OUTPUT}/or/${PREFIX}.d5.enhancers.network.comp.txt

	##### enhancers near each olfactory receptors
	dip-c pd -1 ${or_data}/olfactory_receptors.hom.leg.txt -2 ${or_data}/enhancers.hom.leg.txt \
		${OUTPUT}/conformation/${PREFIX}_20kb_1.clean.3dg > ${OUTPUT}/or/${PREFIX}.olfactory_receptors_vs_enhancers.hom.pd.matrix

	python ${dipc}/network_around.py 2.5 \
		${OUTPUT}/or/${PREFIX}.olfactory_receptors_vs_enhancers.hom.pd.matrix \
		${or_data}/olfactory_receptors.hom.leg.txt ${or_data}/enhancers.hom.leg.txt \
		${or_data}/olfactory_receptors.hom.name.txt ${or_data}/enhancers.hom.name.txt \
		> ${OUTPUT}/or/${PREFIX}.d2.5.olfactory_receptors_vs_enhancers.network.txt

	python ${dipc}/network_around.py 5 \
		${OUTPUT}/or/${PREFIX}.olfactory_receptors_vs_enhancers.hom.pd.matrix \
		${or_data}/olfactory_receptors.hom.leg.txt ${or_data}/enhancers.hom.leg.txt \
		${or_data}/olfactory_receptors.hom.name.txt ${or_data}/enhancers.hom.name.txt \
		> ${OUTPUT}/or/${PREFIX}.d5.olfactory_receptors_vs_enhancers.network.txt

	##### Olfactory receptors near each olfactory receptor

	dip-c pd -1 ${or_data}/olfactory_receptors.hom.leg.txt \
		${OUTPUT}/conformation/${PREFIX}_20kb_1.clean.3dg > ${OUTPUT}/or/${PREFIX}.olfactory_receptors.hom.pd.matrix

	python ${dipc}/network.py ${OUTPUT}/or/${PREFIX}.olfactory_receptors.hom.pd.matrix \
		${or_data}/olfactory_receptors.hom.leg.txt ${or_data}/olfactory_receptors.hom.name.txt \
		> ${OUTPUT}/or/${PREFIX}.d2.5.olfactory_receptors.network.comp.txt

	python ${dipc}/network_5.py ${OUTPUT}/or/${PREFIX}.olfactory_receptors.hom.pd.matrix \
		${or_data}/olfactory_receptors.hom.leg.txt ${or_data}/olfactory_receptors.hom.name.txt \
		> ${OUTPUT}/or/${PREFIX}.d5.olfactory_receptors.network.comp.txt

	python ${dipc}/network_around.py 2.5 \
		${OUTPUT}/or/${PREFIX}.olfactory_receptors.hom.pd.matrix \
		${or_data}/olfactory_receptors.hom.leg.txt ${or_data}/olfactory_receptors.hom.leg.txt \
		${or_data}/olfactory_receptors.hom.name.txt ${or_data}/olfactory_receptors.hom.name.txt \
		> ${OUTPUT}/or/${PREFIX}.d2.5.olfactory_receptors.network.txt

	python ${dipc}/network_around.py 5 \
		${OUTPUT}/or/${PREFIX}.olfactory_receptors.hom.pd.matrix \
		${or_data}/olfactory_receptors.hom.leg.txt ${or_data}/olfactory_receptors.hom.leg.txt \
		${or_data}/olfactory_receptors.hom.name.txt ${or_data}/olfactory_receptors.hom.name.txt \
		> ${OUTPUT}/or/${PREFIX}.d5.olfactory_receptors.network.txt

	## See https://github.com/tanlongzhi/dip-c/issues/20
	## awk: awk -F"," '{OFS="\t"}{$1=$3?$1"(mat)":$1"(pat)";print $1, $2}' 
	## awk -F"," '{print $1}' ${or_data}/enhancers.hom.leg.txt | sed "s/chr//g; s/X/20/g; s/Y/21/g;"  > enhancers.hom.color
	## enhancers cif

	dip-c pos -l ${or_data}/enhancers.hom.leg.txt ${OUTPUT}/conformation/${PREFIX}_20kb_1.clean.3dg > ${OUTPUT}/or/${PREFIX}.enhancers.hom.pos
	paste ${or_data}/enhancers.hom.name.txt ${or_data}/enhancers.hom.color ${OUTPUT}/or/${PREFIX}.enhancers.hom.pos \
		| python ${dipc}/name_color_x_y_z_to_cif.py /dev/stdin > ${OUTPUT}/mmCIF/${PREFIX}.enhancers.hom.cif
	
	# olfactory receptors cif
	dip-c pos -l ${or_data}/olfactory_receptors.hom.leg.txt ${OUTPUT}/conformation/${PREFIX}_20kb_1.clean.3dg > ${OUTPUT}/or/${PREFIX}.olfactory_receptors.hom.pos
	paste ${or_data}/olfactory_receptors.hom.name.txt  ${or_data}/olfactory_receptors.hom.color ${OUTPUT}/or/${PREFIX}.olfactory_receptors.hom.pos \
		| python ${dipc}/name_color_x_y_z_to_cif.py /dev/stdin > ${OUTPUT}/mmCIF/${PREFIX}.olfactory_receptors.hom.cif

	# calculate single-cell chromatin compartment values along the genome
	dip-c color -c ${folder}/mm10.cpg.20k.txt -s3 ${OUTPUT}/conformation/${PREFIX}_20kb_1.clean.3dg > ${OUTPUT}/features/${PREFIX}_cpg_s3.color
fi

echo "$(date) Done!"
