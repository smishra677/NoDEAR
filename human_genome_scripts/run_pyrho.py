import os
import time
import pandas as pd


time_=[]

'''
Demographic history: ACB_pop_sizes.csv 


'''

start_time= time.time()
os.system('pyrho make_table -n 20 -N 21 --numthreads 100 --mu 1e-8 --logfile . --outfile Lookup_human_fin_lookuptable.hdf --approx --smcpp_file ACB_pop_sizes.csv --decimate_rel_tol 0.1')
os.system('pyrho hyperparam -n 20 --numthreads 100 --mu 1e-8 --blockpenalty 50,100 --windowsize 25,50 --logfile . --tablefile Lookup_human_fin_lookuptable.hdf --num_sims 6 --smcpp_file ACB_pop_sizes.csv --outfile Lookup_human_fin_hyperparam_results.txt')
time_.append(time.time()-start_time)
print(time_)
print('$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$')

start_time= time.time()
for k in range(3408):
    print(k,str(k*(50000))+'-'+str((k+1)*50000))

    
    simNum='output_'+str(k*(50000))+'_'+str((k+1)*50000)
    os.system('pyrho optimize --tablefile Lookup_human_fin_lookuptable.hdf  --vcffile '+str(simNum)+'.vcf --outfile ./results/'+str(simNum)+'.rmap --blockpenalty 50 --windowsize 50 --logfile .')
    
time_.append(time.time()-start_time)

print(time_)
pd.DataFrame({'time':time_}).to_csv('time.csv')

