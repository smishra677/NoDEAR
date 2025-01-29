
import os
import time


sim = {'eqtes_50k': 120}


time_=[]


for ii in sim:
    
    os.chdir('./'+ii)
    time_start= time.time()
    demography='eq_pop_size'
    os.system('pyrho make_table -n 20 -N 21 --numthreads 100 --mu 1e-8 --logfile . --outfile eq_lookuptable.hdf --approx --smcpp_file '+demography+'.csv --decimate_rel_tol 0.1')
    os.system('pyrho hyperparam -n 20 --numthreads 100 --mu 1e-8 --blockpenalty 50,100 --windowsize 25,50 --logfile . --tablefile eq_lookuptable.hdf --num_sims 6 --smcpp_file '+demography+'.csv --outfile eq_hyperparam_results.txt')
    time_.append(time.time()-time_start)
    time_start= time.time()

    for simNum in range(sim[ii]):

        
        print(simNum)

        os.system('pyrho optimize --tablefile eq_lookuptable.hdf  --vcffile '+str(simNum)+'.vcf --outfile '+str(simNum)+'.rmap --blockpenalty 50 --windowsize 50 --logfile .')
       
    time_.append(time.time()-time_start)

print(time_)