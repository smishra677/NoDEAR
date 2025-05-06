SIMULATE="ReLERNN_SIMULATE"
TRAIN="ReLERNN_TRAIN"
PREDICT="ReLERNN_PREDICT"
BSCORRECT="ReLERNN_BSCORRECT"
SEED="42"
MU="1e-8"
URTR="9"
DIR="./test_3"
GENOME="./genome.bed"
Train_VCF="./0.vcf"

# Simulate data
${SIMULATE} \
    --vcf ${Train_VCF} \
    --genome ${GENOME} \
    --projectDir ${DIR} \
    --assumedMu ${MU} \
    --upperRhoThetaRatio ${URTR} \
    --nTrain 100000 \
    --nVali 2000\
    --nTest 100 \
    --seed ${SEED}

# Train network
${TRAIN} \
    --projectDir ${DIR} \
	--gpuID "1"\
    --nEpochs "205"\
    --seed ${SEED}


