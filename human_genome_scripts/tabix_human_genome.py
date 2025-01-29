import os

oep =open('/N/slate/samishr/Recombination_landscape_9_10/XGBOOST/NoDear/url_list_fin.txt')
oep_r =oep.read()
oep.close()
url=[]
for ie in oep_r.strip().split('\n'):
    url.append(ie.split(' 	')[0])

urls = []
chromosome_names = []

for ii in url:
    chromosome_names.append(ii.split('chr')[1].split('.')[0])

os.environ["PATH"] = "/local/bin/bcftools:" + os.environ["PATH"]
rejected_interval=[]
time_=[]
for ur in url[:1]:
    pred_hu=[]
    chromosome_names= ur.split('chr')[1].split('.')[0]
    print(chromosome_names,ur)
    if ur[-4:]=='.tbi':
        ur=ur[:-4]
    
    count_snp=[]
    interval=[]   
    window_size=50000
    for k in range(0,3408):
        start= str(k*(50000))
        end=str((k+1)*50000)
                  
        print(str(k*(50000))+'-'+str((k+1)*50000))
        #print('tabix -h  	'+ur+' '+chromosome_names+':'+start+'-'+end+' > /N/slate/samishr/Recombination_landscape_9_10/XGBOOST/NoDear/human_genome_1/output_'+str(k*(50000))+'_'+str((k+1)*50000)+'.vcf')
        os.system('tabix -h  	'+ur+' '+chromosome_names+':'+start+'-'+end+' > /N/slate/samishr/Recombination_landscape_9_10/XGBOOST/NoDear/human_genome_1/output_'+str(k*(50000))+'_'+str((k+1)*50000)+'.vcf')
 
