import msprime as msp
import numpy as np
import os
import pandas as pd
import matplotlib.pyplot as plt
from IPython.display import SVG, display

N0 = 70000
sim = {'popspliteq_10k': 120,'popspliteq_50k': 120,'popspliteq_100k': 120,'popspliteq_150k': 120,'popspliteq_200k': 120}

def np_har(n):
    su_ = []
    for i in range(1, n, 1):
        su_.append(1/i)
    return sum(su_)

for i in sim:
    print(i)
    os.chdir('./' + i)
    file_rate = str(i) + '_rate.npy'
    file_rate_recom = str(i) + '_rec_rate.npy'
    lis = []
    lis1 = []
    lia = []
    #rate_array_1 = np.random.uniform(low=10e-13, high=10e-12, size=20)
    #rate_array_2 = np.random.uniform(low=10e-12, high=10e-11, size=20)
    rate_array_3 = np.random.uniform(low=10e-11, high=10e-10, size=40)
    rate_array_4 = np.random.uniform(low=10e-10, high=10e-09, size=40)
    rate_array_5 = np.random.uniform(low=10e-09, high=10e-08, size=40)
    #rate_array_6 = np.random.uniform(low=10e-08, high=10e-07, size=20)
    rate_array = np.concatenate([rate_array_3, rate_array_4, rate_array_5], axis=0)

    for simNum in range(sim[i]):
        changes = []
        changes_rate = str(simNum) + '.csv'
        N0 = 70000
        demography = msp.Demography()
        demography.add_population(name="A", initial_size=70000.0)
        demography.add_population(name="B", initial_size=70000)
        demography.add_population(name="C", initial_size=70000)
        demography.add_population_split(time=10000, derived=["B", "C"], ancestral="A")
        #ts = msprime.sim_ancestry(samples={"A": 1, "B": 1}, demography=demography, random_seed=12)


        print(i, simNum)

        file_h = str(simNum) + '_haps.npy'
        file_P = str(simNum) + '_pos.npy'


        #tsa = msp.sim_ancestry(samples={"B": 10, "C": 10}, ploidy=1, demography=demography, sequence_length=50000, recombination_rate=rate_array[simNum], model=msp.StandardCoalescent())
        #ts = msp.sim_mutations(tsa, rate=1e-8, model=msp.InfiniteSites())

        if i.split('_')[1]=='10k':
            print('10k')
            tsa=msp.sim_ancestry(samples={"B": 10, "C": 10},ploidy=1,demography=demography,sequence_length=10000,recombination_rate=rate_array[simNum],model=msp.StandardCoalescent())
            ts= msp.sim_mutations(tsa,rate=1e-8,model=msp.InfiniteSites())
        elif i.split('_')[1]=='50k':
            print('50k')
            tsa=msp.sim_ancestry(samples={"B": 10, "C": 10},ploidy=1, demography=demography,sequence_length=50000,recombination_rate=rate_array[simNum],model=msp.StandardCoalescent())
            ts= msp.sim_mutations(tsa,rate=1e-8,model=msp.InfiniteSites())
        
        elif i.split('_')[1]=='100k':
            print('100k')
            tsa=msp.sim_ancestry(samples={"B": 10, "C": 10},ploidy=1,demography=demography,sequence_length=100000,recombination_rate=rate_array[simNum],model=msp.StandardCoalescent())
            ts= msp.sim_mutations(tsa,rate=1e-8,model=msp.InfiniteSites())
        
        elif i.split('_')[1]=='150k':
            print('150k')
            tsa=msp.sim_ancestry(samples={"B": 10, "C": 10},ploidy=1,demography=demography,sequence_length=150000,recombination_rate=rate_array[simNum],model=msp.StandardCoalescent())
            ts= msp.sim_mutations(tsa,rate=1e-8,model=msp.InfiniteSites())
            
        elif i.split('_')[1]=='200k':
            print('200k---')
            tsa=msp.sim_ancestry(samples={"B": 10, "C": 10},ploidy=1,demography=demography,sequence_length=200000,recombination_rate=rate_array[simNum],model=msp.StandardCoalescent())
            ts= msp.sim_mutations(tsa,rate=1e-8,model=msp.InfiniteSites())
            
        H = ts.genotype_matrix()
        lia.append(ts.num_sites)
        P = np.array([s.position for s in ts.sites()], dtype='float32')

        rho = ((ts.num_trees) / np_har(20))

        with open(str(simNum) + ".vcf", "w") as vcf_file:
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

        print(rho)
        lis1.append(rate_array[simNum])
        lis.append(rho)
        np.save(file_h, H)

    plt.clf()
    plt.scatter(lis, lia)
    plt.savefig(str(i) + '.png')
    arr = np.array(lis)
    np.save(file_rate, arr)

    arr1 = np.array(lis1)
    np.save(file_rate_recom, arr1)
    os.chdir('..')



