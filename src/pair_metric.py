#!/usr/bin/env python 
#SBATCH -c 2
#SBATCH --mem 32G
#SBATCH --time 5-00:00:00
#SBATCH -J JJJ
#SBATCH -o temp/slurm-%j.out 


import os, sys
import numpy  as np 
import pandas as pd 

def read_pairs(pairs_path):
	pairs = pd.read_csv(pairs_path, compression='gzip', sep='\t', comment='#', header=None)
	pairs.columns = [ 'readID', 'chrom1', 'pos1', 'chrom2', 'pos2', 'strand1', 'strand2', 'phase1', 'phase2']
	return pairs 

def read_repli_data(repli_path='/dshare/xielab/analysisdata/ThreeDGene/zhangjk/3DAnalysis/2023-12-02_Conc_dipc2/proc_data/repli_data/GM12878_repli_chip.bg'):
	repli_data = pd.read_csv(repli_path, sep='\t',
		names=['chrom', 'start', 'end', 'score'])
	repli_data['score'] = repli_data['score'] - np.quantile(repli_data['score'], 0.7)
	return repli_data 

def bin_pair_dist(pairs):
	def dist_windows(base=2, nbins=150):
		return [ 1000*base**(0.125*i) for i in range(nbins) ]
	dist_bins = dist_windows()
	cis_pairs = pairs[pairs.chrom1 == pairs.chrom2]
	cis_pairs['distance'] = cis_pairs['pos2'] - cis_pairs['pos1']
	hist_count, dist_hist = np.histogram(cis_pairs['distance'], dist_bins)
	return np.asarray(list(zip(dist_hist[1:], hist_count)))

def compute_repli_score(pairs, repli_data, extended_range=1e4, chromSize='/share/home/zhangjk/data/hsp_genome/hg38.chrom.size'):

	import pybedtools as pbt 

	repli_pbt = pbt.BedTool.from_dataframe(repli_data)

	anchorOne = pairs[['chrom1', 'pos1']]
	anchorOne.columns = ['chrom', 'start']
	anchorTwo = pairs[['chrom2', 'pos2']]
	anchorTwo.columns = ['chrom', 'start']
	anchors = pd.concat([anchorOne, anchorTwo], ignore_index=True)

	anchors = anchors.sort_values(['chrom', 'start'], ascending=True).drop_duplicates(subset=['chrom', 'start'], keep='first')
	anchors['end'] = anchors['start'] + 1
	anchors_pbt = pbt.BedTool.from_dataframe(anchors).slop(b=extended_range, g=chromSize)
	overlaps = anchors_pbt.intersect(repli_pbt, wao=True).to_dataframe(names=['chrom', 'start', 'end', 'chrom1', 'start1', 'end1', 'score', 'overlapping_base'])
	overlaps = overlaps[['chrom', 'start', 'end', 'score']]
	overlaps['score'] = overlaps['score'].replace('.', np.nan).astype(float)
	mean_repli_score = overlaps.groupby(['chrom', 'start', 'end'])['score'].mean().reset_index()
	repli_scores = mean_repli_score['score'].dropna()
	repli_status = round((repli_scores>0).sum() / float(repli_scores.shape[0]) * 100, 3)

	return repli_status 

def summary_pair_dist(pairs_path):

	pairs = read_pairs(pairs_path)
	samplename = os.path.basename(pairs_path).replace('.allValidPairs.gz', '')
	cis_pairs = pairs[pairs.chrom1 == pairs.chrom2].copy(deep=True)
	cis_pairs['distance'] = np.abs(cis_pairs['pos2'] - cis_pairs['pos1'])

	n_total = pairs.shape[0]
	n_cis   = cis_pairs.shape[0]
	n_near    = cis_pairs[(cis_pairs.distance >= 25000) & (cis_pairs.distance < 2e6)].shape[0]
	n_mitotic = cis_pairs[(cis_pairs.distance >= 2e6) & (cis_pairs.distance < 1.2e7)].shape[0]
	perc_near = n_near / float(n_near + n_mitotic)
	perc_mitotic = n_mitotic / float(n_near + n_mitotic)
	farAvgDist = np.mean(cis_pairs[cis_pairs.distance > 4.5e6]['distance'])
	repli_data = read_repli_data()
	repli_score = compute_repli_score(pairs, repli_data)
	
	return np.asarray(list([samplename, n_total, n_cis, perc_near, perc_mitotic, farAvgDist, repli_score]))

# /dshare/xielab/analysisdata/ThreeDGene/zhangjk/3DAnalysis/2023-12-02_Conc_dipc2/dock/cell_012/cell_012.pairs.gz
sample = sys.argv[1]
pairs_path = f'/dshare/xielab/analysisdata/ThreeDGene/zhangjk/3DAnalysis/2023-12-02_Conc_dipc2/dock/{sample}/{sample}.allValidPairs.gz'
# dist_count = bin_pair_dist(pairs_path)
# np.savetxt(f'/dshare/xielab/analysisdata/ThreeDGene/zhangjk/3DAnalysis/2023-12-02_Conc_dipc2/dock/{sample}/{sample}_pairs_dist_binning.tsv',
# 	dist_count, delimiter='\t')
dist_summary = summary_pair_dist(pairs_path)
np.savetxt(f'/dshare/xielab/analysisdata/ThreeDGene/zhangjk/3DAnalysis/2023-12-02_Conc_dipc2/dock/{sample}/{sample}_pairs_dist_summary.tsv',
	dist_summary, delimiter='\t', fmt='%s')
