#!/bin/bash
case $1 in
"mkdir")
    mkdir raw ft sip fasta
    ;;
"convert")
    source ~/miniconda3/etc/profile.d/conda.sh
    conda activate mono
    mono YOUR_PATH/Raxport.exe -i raw -o ft -j 40
    ;;
"clean")
    rm -r sip/*
    ;;
"run")
    source ~/miniconda3/etc/profile.d/conda.sh
    conda activate py2
    printf "\n=====Make decoy database=====\n\n"
    python $YOUR_PATH/EnsembleScripts/sipros_prepare_protein_database.py \
        -i fasta/db.faa \
        -o fasta/Decoy.fasta \
        -c SiprosEnsembleConfig.cfg
    printf "\n=====Search database=====\n\n"
    export OMP_NUM_THREADS=10
    files=(ft/*.FT2)
    echo "${files[@]}" | xargs -n 1 -P 9 \
        bash -c 'YOUR_PATH/SiprosEnsembleOMP -f $0 -c SiprosEnsembleConfig.cfg -o sip'
    wait
    printf "\n=====Filter PSM=====\n\n"
    python $YOUR_PATH/EnsembleScripts/sipros_psm_tabulating.py \
        -i sip \
        -c SiprosEnsembleConfig.cfg -o sip
    python $YOUR_PATH/EnsembleScripts/sipros_ensemble_filtering.py \
        -i sip/*.tab \
        -c SiprosEnsembleConfig.cfg \
        -o sip
    printf "\n====Assemble protein=====\n\n"
    python $YOUR_PATH/EnsembleScripts/sipros_peptides_assembling.py \
        -c SiprosEnsembleConfig.cfg \
        -w sip
    printf "\n====Refine FDR=====\n\n"
    conda activate r
    Rscript $YOUR_PATH/refineProteinFDR.R \
        -pro *.pro.txt \
        -psm *.psm.txt \
        -fdr 0.005 \
        -o YOUR_FILE_NAME
    Rscript $YOUR_PATH/getSpectraCountInEachFT.R \
        -pro YOUR_FILE_NAME.proRefineFDR.txt \
        -psm YOUR_FILE_NAME.psm.txt \
        -o YOUR_FILE_NAME
    ;;
*)
    ./cmd.sh "run"
    ;;
esac
