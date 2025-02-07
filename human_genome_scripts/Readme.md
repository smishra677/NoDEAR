# Human Genome

## Requirements
- Install `bcftools`.

## Steps

1. **Run `tabix_human_genome.py`:** This will retrieve all the 50kb windows for chromosome 6 from the finish dataset.  
   You can find the URL in `url_list_fin.txt` (the first entry).

2. **Run `filter_vcf.py`:** This will extract the VCF for 10 individuals (listed in `sample.txt`) and reindex the VCF.

3. **Run `run_pyrho.py`:** This will run Pyrho on the extracted dataset.  
   `ACB_pop_sizes.csv` contains the demographic history of the population and should be used to generate the lookup table.
