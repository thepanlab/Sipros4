#!/bin/bash
source ~/miniconda3/etc/profile.d/conda.sh
binPath=YOUR_BINARY_PATH
case $1 in
"mkdir")
    mkdir configs raw ft sip fasta
    ${binPath}/configGenerator -i C13.cfg -o configs -e C
    ;;
"mkdb")
    conda activate r
    Rscript ${binPath}/makeDBforLabelSearch.R \
        -pro YOUR_PATH/*.pro.txt \
        -faa YOUR_PATH/db.fasta \
        -o fasta/db.faa
    ;;
"convert")
    source ~/miniconda3/etc/profile.d/conda.sh
    conda activate mono
    mono ${binPath}/Raxport.exe -i raw -o ft
    ;;
"clean")
    rm -r sip/* configs/*
    ;;
"run")
    echo "Search database"
    conda activate mpi
    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${CONDA_PREFIX}/lib
    export OMP_NUM_THREADS=10
    time mpirun -np 16 $binPath/SiprosV3mpi \
        -g configs \
        -w ft \
        -o sip
    echo "Filter PSM"
    conda activate py2
    time python ${binPath}/scripts/sipros_peptides_filtering.py \
        -c C13.cfg \
        -w sip
    wait
    python ${binPath}/scripts/sipros_peptides_assembling.py \
        -c C13.cfg \
        -w sip
    wait
    python ${binPath}/scripts/ClusterSip.py \
        -c C13.cfg \
        -w sip
    wait
    Rscript ${binPath}/refineProteinFDR.R \
        -pro sip/*.pro.txt \
        -psm sip/*.psm.txt \
        -fdr 0.01 \
        -o sip/YOUR_FILE_NAME
    wait
    Rscript ${binPath}/getLabelPCTinEachFT.R \
        -pro sip/YOUR_FILE_NAME.proRefineFDR.txt \
        -psm sip/*.psm.txt \
        -thr 5 \
        -o sip/YOUR_FILE_NAME
    filter() {
        mkdir $1
        cp sip/*$1*.sip $1
        python ${binPath}/scripts/sipros_peptides_filtering.py \
            -c C13.cfg \
            -w $1
        wait
        python ${binPath}/scripts/sipros_peptides_assembling.py \
            -c C13.cfg \
            -w $1
        wait
        python ${binPath}/scripts/ClusterSip.py \
            -c C13.cfg \
            -w $1
        wait
    }
    files=(ft/*.FT2)
    files=($(echo "${files[@]}" | tr ' ' '\n' | cut -d '/' -f2 | tr '\n' ' '))
    files=($(echo "${files[@]}" | tr ' ' '\n' | cut -d '.' -f1 | tr '\n' ' '))
    files=($(echo "${files[@]}" | tr ' ' '\n' | sort -u | tr '\n' ' '))
    for folder in ${files[@]}; do
        filter "$folder" &
    done
    wait
    ;;
"filter")
    conda activate python2
    filter() {
        mkdir $1
        cp sip/*$1*.sip $1
        python ${binPath}/scripts/sipros_peptides_filtering.py \
            -c C13.cfg \
            -w $1
        python ${binPath}/scripts/sipros_peptides_assembling.py \
            -c C13.cfg \
            -w $1
        python ${binPath}/scripts/ClusterSip.py \
            -c C13.cfg \
            -w $1
    }
    files=(raw/*.raw)
    files=($(echo "${files[@]}" | tr ' ' '\n' | cut -d '/' -f2 | tr '\n' ' '))
    files=($(echo "${files[@]}" | tr ' ' '\n' | cut -d '.' -f1 | tr '\n' ' '))
    files=($(echo "${files[@]}" | tr ' ' '\n' | sort -u | tr '\n' ' '))
    for folder in ${files[@]}; do
        filter "$folder" &
    done
    wait
    ;;
*)
    ./SIPcmd.sh "run"
    ;;
esac
