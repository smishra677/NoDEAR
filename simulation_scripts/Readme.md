# Simulations

There are 5 files to simulate.

## Requirement
Please create all the folders inside the main folder. The folder names should be structured as follows:

`Type_of_data_size_of_window`

Examples:
- `train_10k` → Training data with a 10kb window.
- `train_10k`, `train_50k`, `train_100k`, `train_150k`, `train_200k`
- `eqtes_10k`, `eqtes_50k`, `eqtes_100k`, `eqtes_150k`, `eqtes_200k`
- `noneqtes_10k`, `noneqtes_50k`, `noneqtes_100k`, `noneqtes_150k`, `noneqtes_200k`
- `geneconversioneq_10k`, `geneconversioneq_50k`, `geneconversioneq_100k`, `geneconversioneq_150k`, `geneconversioneq_200k`
- `popspliteq_10k`, `popspliteq_50k`, `popspliteq_100k`, `popspliteq_150k`, `popspliteq_200k`

## Steps

1. **Run `simulate_train.py`:** Simulates 10,000 data points with various window sizes.
    - You will also get the following files for each of the 10,000 data points:
        - `_haps.npy`: Haplotype matrix (MxN where M is the number of SNPs in the window and N is the number of samples, N=20).
        - `_rec_rate.npy`: Population recombination rate 'c' for the window.
        - `_tree_rate.npy`: Number of tree changes in the window.
        - `_rate.npy`: Value of population recombination rate Rho for the window.
    - It also produces a VCF file for Pyrho.

2. **Run `eq_test.py`, `gene_conversion_test.py`, and `population_split_test.py`** to generate test data.
