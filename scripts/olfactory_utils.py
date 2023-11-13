#!/usr/bin/env python

import os
import numpy as np
import pandas as pd
from scipy import stats
import numpy.ma as ma
from itertools import combinations, product

def analyse_interacting_network(node='GI', distance=2.5):
	assert node == 'GI' or node == 'OR'
	##### change file path of OR expression table 
	expr_df = pd.read_csv('olfr_expression_qualified_meta.with_rank.tsv', sep='\t')
	n = 0
	for sample, group in expr_df.groupby('Cell'):
		group = group[group.homolog!='unknown']
		if group.shape[0] == 0: continue # cells without any OR that can be determined homolog
		if not os.path.exists(f'../port/{sample}/conformation'): continue # cells without structures available
		group['#name'] = group['gene_name'] + '(mat)'
		group.loc[group.homolog!='B6','#name'] = group['gene_name'] + '(pat)'
		if node == 'OR':
			##### change network file of OR or GIs 
			network = pd.read_csv(f'../port/{sample}/or/{sample}.d{distance}.olfactory_receptors.network.txt', sep='\t')
		if node == 'GI':
			network = pd.read_csv(f'../port/{sample}/or/{sample}.d{distance}.olfactory_receptors_vs_enhancers.network.txt', sep='\t')
		network_active = network.merge(group, on='#name', how='inner')
		if n == 0: network_df = network_active
		else: network_df = pd.concat([network_df, network_active], ignore_index=False)
		n += 1
	return network_df

def expressed_or_comp(node='OR', distance=2.5):
	expr_df = pd.read_csv('../proc_data/OR_determining/OR_expression_homologs_determination.with_rank.tsv', sep='\t')
	res = pd.DataFrame(columns=['Cell', 'gene_name', 'FPKM', 'aCount', 'bCount', 'n_variants','homolog', 'Cluster', 'or_type', 'rank', 'or_gene'] + 
		['#n_node', 'n_edge_inter','n_edge_short_intra',  'n_edge_long_intra'])
	for sample in np.unique(expr_df.Cell):
		sample_df = expr_df[expr_df.Cell==sample]
		if not os.path.exists(f'../port/{sample}/conformation'): continue # cells without structures available  
		sample_df = sample_df[sample_df.homolog!='unknown']
		or_genes = [  f'{record.gene_name}(mat)' if record.homolog == 'B6' else  f'{record.gene_name}(pat)' for row, record in sample_df.iterrows() ]
		sample_df['or_gene'] = or_genes 
		network = pd.read_csv(f'../port/{sample}/or/{sample}.d2.5.olfactory_receptors.network.comp.txt', sep='\t')
		for row, component in network.iterrows():
			for or_gene in or_genes:
				if or_gene in component.nodes:
					tmp = sample_df[sample_df['or_gene']==or_gene]
					tmp['#n_node'] = component['#n_node']
					tmp['n_edge_inter'] = component['n_edge_inter']
					tmp['n_edge_short_intra'] = component['n_edge_short_intra']
					tmp['n_edge_long_intra'] = component['n_edge_long_intra']
					res = pd.concat([res, tmp], ignore_index=True)
	return res

def find_enhancer_aggregates(distance=2.5):

	port_path = '/dshare/xielab/analysisdata/ThreeDGene/zhangjk/scRNAC/Olfactory/port'
	n = 0
	for sample in os.listdir(port_path):
		if not os.path.exists(f'{port_path}/{sample}/conformation'): continue # cells without structures available
		network = os.path.join(f'{port_path}/{sample}/or', f'{sample}.d{distance}.olfactory_receptors_vs_enhancers.network.txt')
		network_df = pd.read_csv(network, sep='\t')
		if not np.any(network_df.others): 
			# print(sample)
			continue # no OR-GI interactions
		network_df['n_other'] = network_df.n_other_inter + network_df.n_other_intra
		network_df = network_df.sort_values(['n_other', 'n_other_inter', 'n_other_hom'], ascending=False) # sorting
		network_df.index = range(network_df.shape[0])
		first_cluster = network_df.loc[0] # first enhancer aggregate
		first_cluster_enhancers = set(first_cluster.others.split(','))
		for idx in range(1, network_df.shape[0]):
			candidate_cluster = network_df.loc[idx]
			candidate_cluster_enhancers = set(candidate_cluster .others.split(','))
			enhancer_intersections = first_cluster_enhancers.intersection(candidate_cluster_enhancers)
			if len(enhancer_intersections) == 0:
				second_cluster = candidate_cluster
				break
		first_cluster_df = pd.DataFrame(first_cluster).T
		first_cluster_df['aggregate'] = 'first'
		second_cluster_df = pd.DataFrame(second_cluster).T
		second_cluster_df['aggregate'] = 'second'
		aggregate_df = pd.concat([first_cluster_df, second_cluster_df], ignore_index=True)
		aggregate_df = aggregate_df.rename(columns={'#name': 'center_OR'})
		aggregate_df['Cell'] = sample
		if n == 0: aggregate_tab = aggregate_df
		else: aggregate_tab = pd.concat([aggregate_tab, aggregate_df], ignore_index=True)
		n += 1
	return aggregate_tab

def first_olfr_position(network, gi_cluster, olfr_type='1_dominant'):
	subset_df = network[ (network.olfr_type==olfr_type) & (network['rank']==1) ]
	outputs = list()
	for _, record in subset_df.iterrows():
		first_aggregates  = set(gi_cluster[ (gi_cluster.Cell==record.Cell) & (gi_cluster['aggregate']=='first')].others.values[0].split(','))
		second_aggregates = set(gi_cluster[ (gi_cluster.Cell==record.Cell) & (gi_cluster['aggregate']=='second')].others.values[0].split(','))

		if record.n_other == 0:
			record_output = np.asarray([record.Cell, record['#name'], len(first_aggregates), 0, len(second_aggregates), 0, 0])
		else:
			enhancers = set(record.others.split(','))
			intersection_first  = enhancers.intersection(first_aggregates)
			intersection_second = enhancers.intersection(second_aggregates)
			record_output = np.asarray([record.Cell, record['#name'], len(first_aggregates), len(intersection_first), len(second_aggregates), len(intersection_second), record.n_other])
		outputs.append(record_output)
	outputs_df = pd.DataFrame(np.asarray(outputs), columns=['Cell', '#name', 'first', 'first_intersection', 'second', 'second_intersection', 'n_other'])
	return outputs_df

def compute_gi_strength(data):
	"""
		See 20220301_Candidate_GI/11_enhancers_interchromosomal_contact_strength.html;
		exactly the same as method defined by Tan
	"""
	_, N = data.shape
	x, y = N//2-1, N//2+1
	central_pixel = np.mean(data[x:y, x:y])
	surrounding_mask = np.zeros((N, N))
	surrounding_mask[x:y, x:y] = 1
	mx = ma.masked_array(data, surrounding_mask)
	surrounding_pixel = mx.mean()
	return central_pixel, surrounding_pixel, float(central_pixel)/surrounding_pixel

def compute_or_strength(data):
	""" 
		See 20220301_Candidate_GI/11_enhancers_interchromosomal_contact_strength.html; 
		exactly the same as method defined by Tan 
	"""
	_, N = data.shape
	x, y = N//2-2, N//2+2
	central_pixel = np.mean(data[x:y, x:y])
	surrounding_mask = np.zeros((N, N))
	surrounding_mask[x:y, x:y] = 1
	mx = ma.masked_array(data, surrounding_mask)
	surrounding_pixel = mx.mean()
	return central_pixel, surrounding_pixel, float(central_pixel)/surrounding_pixel

