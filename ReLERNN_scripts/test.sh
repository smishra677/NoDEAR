SIMULATE="ReLERNN_SIMULATE"
TRAIN="ReLERNN_TRAIN"
PREDICT="ReLERNN_PREDICT"
BSCORRECT="ReLERNN_BSCORRECT"
SEED="42"
MU="1e-8"
URTR="9"
DIR="test_3"
OUT_DIR="$DIR/splitVCFs"
GENOME="genome.bed"

 

 

# Prediction
mkdir -p "$DIR/results"

for i in $(seq 0 120); do
    VCF="$i.vcf"

    python ./vcf2hdf5.py --vcf "$VCF" --output "$OUT_DIR/${i}_1:0-50000.hdf5"

    ${PREDICT} \
                                --vcf $$OUT_DIR/${i}_1:0-50000.hdf5 \
                                --projectDir "$DIR" \
                                --seed "$SEED"

                mv "$DIR/${i}_1:0-50000.hdf5" "$DIR/${i}.PREDICT.txt"
            
            ${BSCORRECT} \
			--projectDir ${DIR} \
			--nSlice 100\
            --gpuID "1"\
			--nReps 1000 \
			--seed ${SEED}


        #mv "$NEW_FILE" "$OLD_FILE"

        mv "$DIR/${i}.PREDICT.BSCORRECTED.txt" "$DIR/results/${i}.PREDICT.BSCORRECTED.txt"
        mv "$DIR/${i}.PREDICT.txt" "$DIR/results/${i}.PREDICT.txt"

done
