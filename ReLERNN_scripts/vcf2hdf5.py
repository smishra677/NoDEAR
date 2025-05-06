import h5py

import allel

import argparse

 

def main():

    parser = argparse.ArgumentParser()

    parser.add_argument("--vcf", type=str, required=True)

    parser.add_argument("--output", type=str, required=True)

    args = parser.parse_args()

   

    allel.vcf_to_hdf5(args.vcf,args.output,fields="*",overwrite=True)

 

if __name__ == "__main__":

    main()