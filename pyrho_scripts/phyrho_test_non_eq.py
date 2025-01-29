
import os



sim = {'noneqtes_50k': 120}




for ii in sim:
    os.chdir('./'+ii)
    for simNum in range(sim[ii]):

        sim_file=str(simNum)
        
        print(simNum)
        os.system('pyrho make_table -n 20 -N 21 --numthreads 100 --mu 1e-8 --logfile . --outfile '+str(simNum)+'_lookuptable.hdf --approx --smcpp_file '+sim_file+'.csv --decimate_rel_tol 0.1')
        os.system('pyrho hyperparam -n 20 --numthreads 100 --mu 1e-8 --blockpenalty 50,100 --windowsize 25,50 --logfile . --tablefile '+str(simNum)+'_lookuptable.hdf --num_sims 6 --smcpp_file '+sim_file+'.csv --outfile '+str(simNum)+'_hyperparam_results.txt')
        os.system('pyrho optimize  --numthreads 50 --tablefile '+str(simNum)+'_lookuptable.hdf  --vcffile '+str(simNum)+'.vcf --outfile '+str(simNum)+'.rmap --blockpenalty 50 --windowsize 50 --logfile .')
        #os.system('pyrho compute_r2 --quantiles .25,.5,.75 --compute_mean --samplesize 20 --tablefile '+str(simNum)+'_lookuptable.hdf --outfile '+str(simNum)+'_r2.txt')

