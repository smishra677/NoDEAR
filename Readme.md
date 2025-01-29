# NoDear: No Disequilibrium Estimation of Accurate Recombination

## Requirements
1. XGBoost  
2. msprime  

Please find the model and results in `train-NoDear.ipynb`.  
Scripts to generate the simulated dataset are located in the `./simulation_scripts` folder.  
Scripts to run Pyrho are in the `./pyrho_scripts` folder.  
Scripts to obtain human genome data are in the `./human_genome_scripts` folder.

## Sequence of Actions

1. Simulate data (inside the `./simulation_scripts` folder)  
   *(Please find the README file inside each subfolder.)*

2. Get human genome data (inside the `./human_genome_scripts` folder)  
   *(Please find the README file inside each subfolder.)*

3. Run Pyrho (inside the `./pyrho_scripts` folder)  
   *(Please find the README file inside each subfolder.)*

4. Run `train-NoDear.ipynb` (runs XGBoost, produces results)
