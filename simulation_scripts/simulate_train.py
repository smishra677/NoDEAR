import msprime as msp
import numpy as np
import os
import pandas as pd
import matplotlib.pyplot as plt
from collections import Counter
#msp.simulate(seed,samplesize,ne,mutation,recob)
#ts=msp.simulate(random_seed=42,sample_size=2,Ne=7000,length=10,mutation_rate=10e-2,recombination_rate=10e-2)


def np_har(n):
    su_=[]
    for i in range(1,n,1):
        su_.append(1/i)
    return sum(su_)

sim={'train_10k':10000,'train_50k':10000,'train_50k_100':10000,'train_100k':10000,'train_150k':10000,'train_200k':10000}


for i in sim:
	os.chdir('./'+i)
	file_rate =str(i)+'_rate.npy'
	file_rate_recom =str(i)+'_rec_rate.npy'
	file_rate_tree =str(i)+'_tree_rate.npy'
	
	lis=[]
	lia=[]
	lis1=[]
	lis_tree=[]


	
	rate_array_3 = np.random.uniform(low=5e-12, high=1e-11, size=1000)
	rate_array_4 = np.random.uniform(low=1e-11, high=5e-11, size=1000)
	rate_array_5 = np.random.uniform(low=5e-11, high=1e-10, size=1000)
	rate_array_6 = np.random.uniform(low=1e-10, high=5e-10, size=1000)
	rate_array_7 = np.random.uniform(low=5e-10, high=1e-9, size=1000)
	rate_array_8 = np.random.uniform(low=1e-9, high=5e-9, size=1000)
	rate_array_9 = np.random.uniform(low=5e-9, high=1e-8, size=1000)
	rate_array_10 = np.random.uniform(low=1e-8, high=5e-8, size=1000)
	rate_array_11 = np.random.uniform(low=5e-8, high=1e-7, size=1000)
	rate_array_12 = np.random.uniform(low=1e-7, high=5e-7, size=1000)
	
	
	rate_array = np.concatenate([rate_array_3, rate_array_4, rate_array_5, rate_array_6, rate_array_7, rate_array_8, rate_array_9, rate_array_10, rate_array_11, rate_array_12], axis=0)

	
	n_pop=20
	for simNum in range(sim[i]):
		print(i, simNum)


		file_h = str(simNum)+'_haps.npy'
		file_P = str(simNum)+'_pos.npy'
	


		if i.split('_')[1]=='10k':
			print('10k')
			tsa=msp.sim_ancestry(samples=20,ploidy=1,population_size=70000,sequence_length=10000,recombination_rate=rate_array[simNum],model=msp.StandardCoalescent(),additional_nodes=msp.NodeType.RECOMBINANT | msp.NodeType.GENE_CONVERSION,coalescing_segments_only=False)
			ts= msp.sim_mutations(tsa,rate=1e-8,model=msp.InfiniteSites())
		elif i=='train_50k_100':
			print('train_50k_100')
			tsa=msp.sim_ancestry(samples=100,ploidy=1, population_size=70000,sequence_length=50000,recombination_rate=rate_array[simNum],model=msp.StandardCoalescent(),additional_nodes=msp.NodeType.RECOMBINANT | msp.NodeType.GENE_CONVERSION,coalescing_segments_only=False)
			ts= msp.sim_mutations(tsa,rate=1e-8,model=msp.InfiniteSites())
			n_pop=100

		elif i=='train_50k':
			print('50k')
			tsa=msp.sim_ancestry(samples=20,ploidy=1, population_size=70000,sequence_length=50000,recombination_rate=rate_array[simNum],model=msp.StandardCoalescent(),additional_nodes=msp.NodeType.RECOMBINANT | msp.NodeType.GENE_CONVERSION,coalescing_segments_only=False)
			ts= msp.sim_mutations(tsa,rate=1e-8,model=msp.InfiniteSites())
		elif i.split('_')[1]=='100k':
			print('100k')
			tsa=msp.sim_ancestry(samples=20,ploidy=1,population_size=70000,sequence_length=100000,recombination_rate=rate_array[simNum],model=msp.StandardCoalescent(),additional_nodes=msp.NodeType.RECOMBINANT | msp.NodeType.GENE_CONVERSION,coalescing_segments_only=False)
			ts= msp.sim_mutations(tsa,rate=1e-8,model=msp.InfiniteSites())
		elif i.split('_')[1]=='150k':
			print('150k')
			tsa=msp.sim_ancestry(samples=20,ploidy=1,population_size=70000,sequence_length=150000,recombination_rate=rate_array[simNum],model=msp.StandardCoalescent(),additional_nodes=msp.NodeType.RECOMBINANT | msp.NodeType.GENE_CONVERSION,coalescing_segments_only=False)
			ts= msp.sim_mutations(tsa,rate=1e-8,model=msp.InfiniteSites())
		elif i.split('_')[1]=='200k':
			print('200k---')
			tsa=msp.sim_ancestry(samples=20,ploidy=1,population_size=70000,sequence_length=200000,recombination_rate=rate_array[simNum],model=msp.StandardCoalescent(),additional_nodes=msp.NodeType.RECOMBINANT | msp.NodeType.GENE_CONVERSION,coalescing_segments_only=False)
			ts= msp.sim_mutations(tsa,rate=1e-8,model=msp.InfiniteSites())



		H= ts.genotype_matrix()


		



		np.save(file_h,H)
		lia.append(ts.num_sites)
		
		P = np.array([s.position for s in ts.sites()],dtype='float32')
		print(ts.num_sites)

		node_table = ts.tables.nodes

		flags = node_table.flags

		recombinant_nodes = np.where(flags & msp.NodeType.RECOMBINANT.value)[0]


		gc_nodes = np.where(flags & msp.NodeType.GENE_CONVERSION.value)[0]

		rho = (len(recombinant_nodes)/2)/np_har(n_pop)
		#rho= ((ts.num_trees)/np_har(n_pop))

		with open(str(simNum)+".vcf", "w") as vcf_file:

			vcf_file.write("##fileformat=VCFv4.1\n")
			vcf_file.write("##source=msprime\n")

			vcf_file.write("##contig=<ID=1,length=50000>\n")  
			vcf_file.write("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n")
			vcf_file.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT")
			for ii in range(ts.num_samples // 2):
				vcf_file.write(f"\tdiploid{ii}")
			vcf_file.write("\n")
			
			
			haplotypes = list(ts.haplotypes())
			for variant in ts.variants():
				
				vcf_file.write(f"1\t{int(variant.site.position)}\t.\t{variant.alleles[0]}\t{variant.alleles[1]}\t.\t.\t.\tGT")
				for ij in range(0, ts.num_samples, 2):
					
					diplotype = f"{variant.genotypes[ij]}|{variant.genotypes[ij + 1]}"
					vcf_file.write(f"\t{diplotype}")
				vcf_file.write("\n")

		#with open(str(simNum)+'.vcf', 'w') as file:
			#for var in ts.variants():
				#file.write(f"{var.site.position}\t{var.alleles}\t{var.genotypes}\n")


		lis.append(rho)
		lis_tree.append(ts.num_trees)
		lis1.append(rate_array[simNum])
		#lis.append(rate_array[simNum])
		np.save(file_h,H)
		#np.save(file_P,P)
	
	#import matplotlib.pyplot as plt
	plt.clf()
	plt.scatter(lis,lia)
	plt.savefig(str(i)+'.png')
	#plt.show()
	arr=np.array(lis)
	np.save(file_rate,arr)

	arr2=np.array(lis_tree)
	np.save(file_rate_tree,arr2)


	arr1=np.array(lis1)
	np.save(file_rate_recom,arr1)
	os.chdir('..')
