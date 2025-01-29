import os



'''
Please use the lookup table produced by pyrho_test_eq.py
'''
sim = {'geneconversioneq_50k': 120}




for ii in sim:
    os.chdir('./'+ii)
    for simNum in range(sim[ii]):
        sim_file='geneconversioneq_demography'
        
        print(simNum)

        os.system('pyrho optimize --tablefile geneconversioneq_demography_lookuptable.hdf  --vcffile '+str(simNum)+'.vcf --outfile '+str(simNum)+'.rmap --blockpenalty 50 --windowsize 50 --logfile .')

