import os
'''
sample.txt contains the IDs for 10 finish individuals

HG00271
HG00276
HG00310
HG00315
HG00327
HG00334
HG00339
HG00341
HG00346
HG00353

'''


for k in range(3408):
    print(k,str(k*(50000))+'-'+str((k+1)*50000))
    os.system('bgzip -c output_'+str(k*(50000))+'_'+str((k+1)*50000)+ '.vcf > output_'+str(k*(50000))+'_'+str((k+1)*50000)+ '.vcf.gz')
    os.system('bcftools index output_'+str(k*(50000))+'_'+str((k+1)*50000)+ '.vcf.gz')
    os.system('bcftools view -S sample.txt output_'+str(k*(50000))+'_'+str((k+1)*50000)+ '.vcf.gz -o output_'+str(k*(50000))+'_'+str((k+1)*50000)+'.vcf')
    