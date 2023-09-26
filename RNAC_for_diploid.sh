#!/bin/bash

usage() {
		cat <<-EOF
		Usage: `basename $0` 
		Options:
		    -g: genome.fa
		    -f: R1.fq
		    -r: R2.fq
		    -t: threads
		    -h: help information
	EOF
	exit 1;
}

GENOME=""
R1=""
R2=""
nThreads=8
PREFIX=""
OUTPUT=""
SNP=""
ASSEMBLY=""

# extract options and their arguments into variables
while getopts "g:1:2:o:t:p:s:v:iAPMDR" options; do
	case ${options} in 
		g) GENOME=${OPTARG};;
		1) R1=${OPTARG};;
		2) R2=${OPTARG};;
		o) OUTPUT=${OPTARG};;
		p) PREFIX=${OPTARG};;
		t) nThreads=${OPTARG};;
		v) SNP=${OPTARG};;
		s) ASSEMBLY=${OPTARG};;
		A) ALIGN=1;;
		P) PARSE=1;;
		M) MODEL=1;;
		D) DIPC=1;;
		R) RNAC=1;;
		i) usage;;
		*) exit 1;;
	esac
done

dipc="/share/home/zhangjk/hic_tools/dip-c/scripts"
JUICER="java -Xmx4g -jar /share/home/zhangjk/software/juicer_tools_1.14.08.jar"
colorfolder="/share/home/zhangjk/scRNAC/RNAC_Data/hg38"

# 1. the first step
if [ ${PARSE-0} == 1 ];then

	[[ -z ${OUTPUT} ]] && { echo "OUTPUT: NONE!"; exit 1; }
	[ -d ${OUTPUT} ] || { mkdir -p ${OUTPUT}; }
	[ -d ${OUTPUT}/logs ] || mkdir -p ${OUTPUT}/logs

	echo "$(date)"
	echo "RNAC pipeline for GM12878"
	echo "\n---------------------------------------------------"
	echo "GENOME: ${GENOME}"
	echo "R1: ${R1}"
	echo "R2: ${R2}"
	echo "OUTPUT: ${OUTPUT}"
	echo "PREFIX: ${PREFIX}"
	echo "SNP: ${SNP}"
	echo "\n---------------------------------------------------"

	[[ -z ${GENOME} || ! -s ${GENOME} ]] && { echo "GENOME: not found!"; exit 1; }
	[[ -z ${R1} || ! -s ${R1} ]] && { echo "R1: not found!"; exit 1; }
	[[ -z ${R2} || ! -s ${R2} ]] && { echo "R2: not found!"; exit 1; }

	echo "FASTQ stats"
	echo -e "R1(${R1})\t${PREFIX}\t$(kseq_fastq_base ${R1})"
	echo -e "R2(${R2})\t${PREFIX}\t$(kseq_fastq_base ${R2})"

	# skip aligning reads to reference genome
	if [ ${ALIGN-0} == 1 ];then
		seqtk mergepe ${R1} ${R2} | \
			pre-meta - | \
			bwa mem -5SP -p ${GENOME} -t ${nThreads} - 2>${OUTPUT}/logs/${PREFIX}_mapping.log | \
			gzip > ${OUTPUT}/${PREFIX}.sam.gz
	fi

	# deduplicate regardless of strand information, important!
	hickit.js sam2seg -v ${SNP} ${OUTPUT}/${PREFIX}.sam.gz 2>${OUTPUT}/logs/${PREFIX}_sam2seg.log | hickit.js chronly -y - |  gzip > ${OUTPUT}/${PREFIX}.allValidSeg.gz
	hickit -i ${OUTPUT}/${PREFIX}.allValidSeg.gz -o - 2>${OUTPUT}/logs/${PREFIX}_seg2pairs.log | sed "s/chromosome/chromsize/g" | bgzip > ${OUTPUT}/${PREFIX}.pairs.gz
	zcat ${OUTPUT}/${PREFIX}.pairs.gz | awk -v FS='\t' -v OFS='\t' '{if(substr($0, 1, 1)=="#"){print}else{$6="+";$7="+";print}}' | bgzip > ${OUTPUT}/${PREFIX}.temp.pairs.gz
	hickit --dup-dist=1000 -i ${OUTPUT}/${PREFIX}.temp.pairs.gz -o - 2>${OUTPUT}/logs/${PREFIX}_seg2pairs.log2 | sed "s/chromosome/chromsize/g" | bgzip > ${OUTPUT}/${PREFIX}.allValidPairs.gz

	cut -f1-2 ${GENOME}.fai | grep "^chr" > ${OUTPUT}/hg38.chrom.sizes
	awk '{OFS="\t"}{print $1"(pat)",$2"\n"$1"(mat)",$2}' ${OUTPUT}/hg38.chrom.sizes > ${OUTPUT}/hg38.chr.hom.len
	pairix -p pairs ${OUTPUT}/${PREFIX}.allValidPairs.gz
	pairtools stats -o ${OUTPUT}/logs/${PREFIX}_allValidPairs_stats.txt ${OUTPUT}/${PREFIX}.allValidPairs.gz
	pairsqc.py -p ${OUTPUT}/${PREFIX}.allValidPairs.gz -c ${OUTPUT}/hg38.chrom.sizes -t P -O ${OUTPUT}/${PREFIX}

	# calculate 3dg
	[ -d ${OUTPUT}/structure ] || mkdir -p ${OUTPUT}/structure
	# readlink: get absolute path of output folder
	ln -fs $(readlink -f $OUTPUT)/${PREFIX}.allValidPairs.gz $(readlink -f $OUTPUT)/structure/contacts.pairs.gz
	hickit -i $OUTPUT/structure/contacts.pairs.gz -u -o - 2>$OUTPUT/logs/${PREFIX}_impute.log | bgzip > $OUTPUT/structure/${PREFIX}.impute.pairs.gz 
	#!
	ln -fs $(readlink -f $OUTPUT)/structure/${PREFIX}.impute.pairs.gz $(readlink -f $OUTPUT)/structure/impute.pairs.gz 

	# make .hic and .cool file
	$JUICER pre -n  ${OUTPUT}/${PREFIX}.allValidPairs.gz ${OUTPUT}/${PREFIX}.hic hg38 > ${OUTPUT}/logs/${PREFIX}_juicer.log
	hic2cool convert -r 0 ${OUTPUT}/${PREFIX}.hic ${OUTPUT}/${PREFIX}.mcool > $OUTPUT/logs/${PREFIX}_hic2cool.log
	
	#! 
	${dipc}/hickit_impute_pairs_to_con.sh $OUTPUT/structure/impute.pairs.gz
	${dipc}/hickit_pairs_to_con.sh $OUTPUT/structure/contacts.pairs.gz
	#!
	${dipc}/con_imputed_to_juicer_pre_short.sh ${OUTPUT}/structure/impute.con.gz
	$JUICER pre -n  ${OUTPUT}/structure/impute.juicer.txt.gz $OUTPUT/${PREFIX}.impute.hic ${OUTPUT}/hg38.chr.hom.len > $OUTPUT/logs/${PREFIX}_impute_juicer.log  
	
	source activate dipc && export PYTHONPATH=/share/home/zhangjk/anaconda3/envs/dipc/lib/python2.7/site-packages/
	mkdir -p $OUTPUT/compartment
	dip-c info $OUTPUT/structure/impute.con.gz >$OUTPUT/logs/${PREFIX}.info.log 2>/dev/null
	dip-c info $OUTPUT/structure/contacts.con.gz >>$OUTPUT/logs/${PREFIX}.info.log 2>/dev/null
	dip-c color2 -b 1000000 -H -c $colorfolder/hg38.cpg.1Mb.txt  -s $OUTPUT/structure/contacts.con.gz > $OUTPUT/compartment/${PREFIX}_cpg_b1Mb.color2
	dip-c color2 -b 200000 -H -c $colorfolder/hg38.cpg.200Kb.txt -s $OUTPUT/structure/contacts.con.gz > $OUTPUT/compartment/${PREFIX}_cpg_b200Kb.color2
	dip-c color2 -b 50000 -H -c $colorfolder/hg38.cpg.50Kb.txt -s $OUTPUT/structure/contacts.con.gz > $OUTPUT/compartment/${PREFIX}_cpg_b50Kb.color2
	dip-c color2 -b 20000 -H -c $colorfolder/hg38.cpg.20Kb.txt -s $OUTPUT/structure/contacts.con.gz > $OUTPUT/compartment/${PREFIX}_cpg_b20Kb.color2
	# clean some intermediate files
	rm  ${OUTPUT}/${PREFIX}.temp.pairs.gz
	
fi


if [ ${MODEL-0} == 1 ];then

echo $(date)
echo $OUTPUT
echo $PREFIX

for rep in `seq 4 5`;do
sbatch -J $PREFIX -o $OUTPUT/logs/${PREFIX}_rep${rep}.log <<EOF
#!/bin/bash
#SBATCH -c 1
#SBATCH -J $PREFIX
#SBATCH --mem 8G
#SBATCH -J ${PREFIX}_${rep}
#SBATCH --time 3-00:00:00

hickit -s${rep} -M -i ${OUTPUT}/structure/${PREFIX}.impute.pairs.gz -Sr1m -c1 -r10m -c2 -b4m \
	-b1m -O ${OUTPUT}/structure/${PREFIX}_1Mb_${rep}.3dg \
	-b200k -O ${OUTPUT}/structure/${PREFIX}_200Kb_${rep}.3dg \
	-D5 -b50k -O ${OUTPUT}/structure/${PREFIX}_50Kb_${rep}.3dg \
	-D5 -b20k -O ${OUTPUT}/structure/${PREFIX}_20Kb_${rep}.3dg

# RESARR=("1Mb" "200Kb" "50Kb" "20Kb")
# for RES in \${RESARR[@]};do
#     # ${dipc}/hickit_3dg_to_3dg_rescale_unit.sh ${OUTPUT}/structure/${PREFIX}_\${RES}_${rep}.3dg
#     dip-c clean3 -c ${OUTPUT}/structure/impute.con.gz ${OUTPUT}/structure/${PREFIX}_\${RES}_${rep}.dip-c.3dg > ${OUTPUT}/structure/${PREFIX}_\${RES}_${rep}.clean.3dg
# done
EOF
done

fi

if [ ${DIPC-0} == 1 ];then

	RESARR=("1Mb" "200Kb" "50Kb" "20Kb")
	CON=`ls $OUTPUT/structure/contacts.con.gz`

	# 3.1 dip-c align (RMS RMSD)
	mkdir -p $OUTPUT/alignment 
	for RES in ${RESARR[@]};do
		dip-c align -o $OUTPUT/alignment/bs_${RES}_ $OUTPUT/structure/*_${RES}_[1-3].clean.3dg 2>$OUTPUT/logs/${PREFIX}_rmsd_${RES}.logs > $OUTPUT/alignment/${PREFIX}_${RES}_rmsd.color 
	done

	# 3.2 (relate to Fig. S11. Tan et al., Science 2018)
	mkdir -p $OUTPUT/features

	for RES in 200Kb 50Kb 20Kb;do
		
		S3=`ls $OUTPUT/structure/${PREFIX}_${RES}_1.clean.3dg`
		binsize=`echo $RES | sed "s/Mb/000000/g; s/Kb/000/g;"`

		### Calculate radial positioning
		dip-c color -C $S3 > $OUTPUT/features/${PREFIX}_${RES}_C.color

		### Calculate single-cell chromatin compartment
		dip-c color -c $colorfolder/hg38.cpg.${RES}.txt -s3 $S3 > $OUTPUT/compartment/${PREFIX}_${RES}_cpg_s3.color
		dip-c color2 -b $binsize -H -c $colorfolder/hg38.cpg.${RES}.txt -s $CON > $OUTPUT/compartment/${PREFIX}_cpg_b${RES}.color2

		### Calculate intrachromosomal and interchromosomal neigbors as proxy of chromosome surface
		dip-c color -i3 -S 1000000 $S3 > $OUTPUT/features/${PREFIX}_${RES}_i3_S1mb.color
		dip-c color -I3 -S 1000000  $S3 > $OUTPUT/features/${PREFIX}_${RES}_I3_S1mb.color
		dip-c color -d3  $S3 > $OUTPUT/features/${PREFIX}_${RES}_d3.color 
		dip-c color -r3 $S3 > $OUTPUT/features/${PREFIX}_${RES}_r3.color 

		# imprinting ?
		dip-c color -h $S3 > $OUTPUT/features/${PREFIX}_${RES}_h.color

	done

	# 3.3 dip-c vis
	VISFOLDER="$OUTPUT/visualization"
	VISPREFIX="${PREFIX}_20Kb_1.clean"
	S3=`ls $OUTPUT/structure/${PREFIX}_20Kb_1.clean.3dg`
	
	mkdir -p $VISFOLDER
	DISTANCE=3

	dip-c color -n $colorfolder/hg38.chr.txt $S3 |  dip-c vis -c /dev/stdin $S3 > $VISFOLDER/$VISPREFIX.n.cif
	dip-c color -l $colorfolder/hg38.chr.len $S3 | dip-c vis -c /dev/stdin $S3 > $VISFOLDER/$VISPREFIX.l.cif
	dip-c color -c $colorfolder/hg38.cpg.20k.txt $S3 | dip-c vis -M -c /dev/stdin $S3 > $VISFOLDER/$VISPREFIX.cpg.cif

	# color by arm locus divided by arm length
	dip-c color -L $colorfolder/hg38.chr.cen $S3 | dip-c vis -c /dev/stdin $S3 > $VISFOLDER/$VISPREFIX.cen.cif 
	dip-c vis -M -c $OUTPUT/features/${PREFIX}_20Kb_h.color $S3 > $VISFOLDER/$VISPREFIX.h.cif
	dip-c vis -M -c $OUTPUT/features/${PREFIX}_20Kb_i${DISTANCE}_S1mb.color $S3 > $VISFOLDER/$VISPREFIX.i${DISTANCE}.cif
	dip-c vis -M -c $OUTPUT/features/${PREFIX}_20Kb_I${DISTANCE}_S1mb.color $S3 > $VISFOLDER/$VISPREFIX.I${DISTANCE}.cif
	dip-c vis -M -c $OUTPUT/features/${PREFIX}_20Kb_d${DISTANCE}.color $S3 > $VISFOLDER/$VISPREFIX.d${DISTANCE}.cif
	dip-c vis -M -c $OUTPUT/features/${PREFIX}_20Kb_r${DISTANCE}.color $S3 > $VISFOLDER/$VISPREFIX.r${DISTANCE}.cif

	# expand a nucleus into seprate chromosomes
	dip-c exp $S3 > $VISFOLDER/$VISPREFIX.exp.3dg 2>$VISFOLDER/$VISPREFIX.exp.py
	dip-c color -n  $colorfolder/hg38.chr.txt $VISFOLDER/$VISPREFIX.exp.3dg | dip-c vis -c /dev/stdin $VISFOLDER/$VISPREFIX.exp.3dg > $VISFOLDER/$VISPREFIX.exp.n.cif

fi

if [ ${RNAC-0} == 1 ];then
	
	set -x 

	###
	### Expressed genes taken from corresponding RNA-seq dataset.
	### Expression value: log10(FPKM + 1)
	### 

	# gene's position
	pos="/share/home/zhangjk/scRNAC/RNAC_Data/hg38/hg38.genes.coords.csv"
	random_sampler="python /share/home/zhangjk/scRNAC/scripts/random_sampler.py"
	inner_join="python /share/home/zhangjk/scRNAC/scripts/inner_join.py"

	[ -d $OUTPUT/genes ] || mkdir -p $OUTPUT/genes
	mkdir -p $OUTPUT/visualization
	S3=`ls $OUTPUT/structure/${PREFIX}_20Kb_1.clean.3dg`

	### 1. as Tan et al., NSMB (2019)

	## $OUTPUT/features/${PREFIX}.expr.genes.info.tsv
	## | 1st col | 2nd col | 3rd col | 4th col | 5th col |
	## | gene_id | log10(FPKM) | chrom | tss | gene_name | 
	
	cat $OUTPUT/RNA/${PREFIX}.genes.results | cut -f1,7 | awk '(NR==1 || $2>1)' \
		| tr '\t' ',' | ${inner_join} $pos /dev/stdin | tr ',' '\t' \
		| tail -n+2 | sort -k3,3V -k4,4n | grep -v chrY | awk '{OFS="\t"}{print $1, log($2+1)/log(10), $3, $4, $5}' > $OUTPUT/genes/${PREFIX}.expr.genes.info.tsv
	
	### create leg (expressed)
	awk 'BEGIN{FS="\t"}{print $5"(pat)"; print $5"(mat)"}' $OUTPUT/genes/${PREFIX}.expr.genes.info.tsv > $OUTPUT/genes/${PREFIX}.expr.genes.name
	awk 'BEGIN{FS="\t"}{print $3","$4",0"; print $3","$4",1"}' $OUTPUT/genes/${PREFIX}.expr.genes.info.tsv > $OUTPUT/genes/${PREFIX}.expr.genes.leg
	### create expression color file
	awk 'BEGIN{FS="\t"}{print $2; print $2}' $OUTPUT/genes/${PREFIX}.expr.genes.info.tsv > $OUTPUT/genes/${PREFIX}.expr.genes.fpkm.color
	### Find  position
	dip-c pos -l $OUTPUT/genes/${PREFIX}.expr.genes.leg $S3 > $OUTPUT/genes/${PREFIX}.expr.genes.pos
	paste $OUTPUT/genes/${PREFIX}.expr.genes.name $OUTPUT/genes/${PREFIX}.expr.genes.fpkm.color $OUTPUT/genes/${PREFIX}.expr.genes.pos \
		| python $dipc/name_color_x_y_z_to_cif.py /dev/stdin > $OUTPUT/visualization/${PREFIX}.expr.genes.fpkm.cif

	### 2. structure analysis for genes

	### Expressed 
	awk -F"," '{OFS="\t"}{ hom=$3==0?"(pat)":"(mat)"; print $1hom, $2 }' $OUTPUT/genes/${PREFIX}.expr.genes.leg \
		| paste /dev/stdin $OUTPUT/genes/${PREFIX}.expr.genes.name \
		| paste /dev/stdin $OUTPUT/genes/${PREFIX}.expr.genes.pos \
		| sort -k1,1V -k2,2n \
		| tee >(cut -f3 > $OUTPUT/genes/${PREFIX}.expr.genes.sorted.name) \
		| cut -f1,2,4-6 > $OUTPUT/genes/${PREFIX}.expr.genes.3dg

	paste $OUTPUT/genes/${PREFIX}.expr.genes.sorted.name $OUTPUT/genes/${PREFIX}.expr.genes.3dg \
		| awk '{print $1, $4, $5, $6, $2$3}' \
		| sed '1i gene_name x y z locus' | tr ' ' ',' > $OUTPUT/genes/${PREFIX}.expr.genes.sorted.csv

	dip-c color -i 3 $OUTPUT/genes/${PREFIX}.expr.genes.3dg > $OUTPUT/genes/${PREFIX}.expr.genes.i3.color 
	dip-c color -d 3 $OUTPUT/genes/${PREFIX}.expr.genes.3dg > $OUTPUT/genes/${PREFIX}.expr.genes.d3.color
	dip-c color -r 3 $OUTPUT/genes/${PREFIX}.expr.genes.3dg > $OUTPUT/genes/${PREFIX}.expr.genes.r3.color

	### create mmcif file
	for f in $OUTPUT/genes/${PREFIX}.expr.genes.[idr]3.color;do
		awk 'BEGIN{OFS=","; print "locus,color"}{print $1$2, $3}' $f \
			| ${inner_join} $OUTPUT/genes/${PREFIX}.expr.genes.sorted.csv /dev/stdin \
			| tail -n+2 | tr ',' '\t' | sed 's/"//g' \
			| tee >(awk '{OFS="\t"}{print $3,$2}'>${f/.color/.tsv}) \
			| awk '{OFS="\t"}{print $3, $2, $4, $5, $6}' \
			| python $dipc/name_color_x_y_z_to_cif.py /dev/stdin > ${f/.color/.cif}
		mv ${f/.color/.tsv} $OUTPUT/features/
		mv ${f/.color/.cif} $OUTPUT/visualization
	done

	### 3. Repeat procedures
	N=$(cat $OUTPUT/genes/${PREFIX}.expr.genes.info.tsv | wc -l)

	### Random subset
	( head -n1 $OUTPUT/RNA/${PREFIX}.genes.results && ${random_sampler} $OUTPUT/RNA/${PREFIX}.genes.results $N ) | cut -f1,7 | tr '\t' ',' \
		| ${inner_join} $pos /dev/stdin | sed 's/"//g'  | tr ',' '\t' | tail -n+2 \
		| sort -k3,3V -k4,4n | grep -v chrY | awk '{OFS="\t"}{print $1, log($2+1)/log(10), $3, $4, $5}' > $OUTPUT/genes/${PREFIX}.random.run01.genes.tsv
	### format to leg
	awk 'BEGIN{FS="\t"}{print $3","$4",0"; print $3","$4",1"}' $OUTPUT/genes/${PREFIX}.random.run01.genes.tsv > $OUTPUT/genes/${PREFIX}.random.run01.genes.leg
	awk 'BEGIN{FS="\t"}{print $5"(pat)"; print $5"(mat)"}' $OUTPUT/genes/${PREFIX}.random.run01.genes.tsv > $OUTPUT/genes/${PREFIX}.random.run01.genes.name
	awk 'BEGIN{FS="\t"}{print $2; print $2}' $OUTPUT/genes/${PREFIX}.random.run01.genes.tsv > $OUTPUT/genes/${PREFIX}.random.run01.genes.fpkm.color
	dip-c pos -l $OUTPUT/genes/${PREFIX}.random.run01.genes.leg $S3 > $OUTPUT/genes/${PREFIX}.random.run01.genes.pos
	paste $OUTPUT/genes/${PREFIX}.random.run01.genes.name $OUTPUT/genes/${PREFIX}.random.run01.genes.fpkm.color $OUTPUT/genes/${PREFIX}.random.run01.genes.pos \
		| python $dipc/name_color_x_y_z_to_cif.py /dev/stdin > $OUTPUT/visualization/${PREFIX}.random.run01.genes.fpkm.cif

	awk -F"," '{OFS="\t"}{ hom=$3==0?"(pat)":"(mat)"; print $1hom, $2 }' $OUTPUT/genes/${PREFIX}.random.run01.genes.leg \
		| paste /dev/stdin $OUTPUT/genes/${PREFIX}.random.run01.genes.name \
		| paste /dev/stdin $OUTPUT/genes/${PREFIX}.random.run01.genes.pos \
		| sort -k1,1V -k2,2n \
		| tee >(cut -f3 > $OUTPUT/genes/${PREFIX}.random.run01.genes.sorted.name) \
		| cut -f1,2,4-6 > $OUTPUT/genes/${PREFIX}.random.run01.genes.3dg

	paste $OUTPUT/genes/${PREFIX}.random.run01.genes.sorted.name $OUTPUT/genes/${PREFIX}.random.run01.genes.3dg \
		| awk '{print $1, $4, $5, $6, $2$3}' \
		| sed '1i gene_name x y z locus' | tr ' ' ',' > $OUTPUT/genes/${PREFIX}.random.run01.genes.sorted.csv

	dip-c color -i 3 $OUTPUT/genes/${PREFIX}.random.run01.genes.3dg > $OUTPUT/genes/${PREFIX}.random.run01.genes.i3.color 
	dip-c color -d 3 $OUTPUT/genes/${PREFIX}.random.run01.genes.3dg > $OUTPUT/genes/${PREFIX}.random.run01.genes.d3.color
	dip-c color -r 3 $OUTPUT/genes/${PREFIX}.random.run01.genes.3dg > $OUTPUT/genes/${PREFIX}.random.run01.genes.r3.color   

	for f in $OUTPUT/genes/${PREFIX}.random.run01.genes.[idr]3.color;do
		awk 'BEGIN{OFS=","; print "locus,color"}{print $1$2, $3}' $f \
			| ${inner_join}  $OUTPUT/genes/${PREFIX}.random.run01.genes.sorted.csv /dev/stdin \
			| tail -n+2 | tr ',' '\t' | sed 's/"//g' \
			| tee >(awk '{OFS="\t"}{print $3,$2}'>${f/.color/.tsv}) \
			| awk '{OFS="\t"}{print $3, $2, $4, $5, $6}' \
			| python $dipc/name_color_x_y_z_to_cif.py /dev/stdin > ${f/.color/.cif}
		mv ${f/.color/.tsv} $OUTPUT/features/
		mv ${f/.color/.cif} $OUTPUT/visualization/
	done
	set +x
fi

echo "$(date) Done!"
