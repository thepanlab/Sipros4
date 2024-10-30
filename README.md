### SiprosToolKits

This repository contains tools for stable isotopic mass spectrometry-based metaproteomics research developed by Sipros team. These include Raxport, SiprosV4, SiprosEnsemble and some scripts

[You can find the simple tutorial for 13C-labeled E. coli on our wiki page](https://github.com/thepanlab/SiprosToolKits/wiki/13C-labeled-E.-coli-SIP-proteomic-search-tutorial)

### Citation

Xiong, Yi, Ryan S. Mueller, Shichao Feng, Xuan Guo, and Chongle Pan. "Proteomic stable isotope probing with an upgraded Sipros algorithm for improved identification and quantification of isotopically labeled proteins." Microbiome 12 (2024).

### install sipros conda environment (recommand)

[sipros on bioconda](https://anaconda.org/bioconda/sipros)

```bash
conda create -n sipros bioconda::sipros
conda activate sipros
# get help
copyConfigTemplate -h
Raxport -h
SiprosEnsembleOMP -h
EnsembleScripts_sipros_prepare_protein_database -h
EnsembleScripts_sipros_psm_tabulating -h
EnsembleScripts_sipros_ensemble_filtering -h
EnsembleScripts_sipros_peptides_assembling -h
V4Scripts_refineProteinFDR -help
V4Scripts_getSpectraCountInEachFT -help
V4Scripts_makeDBforLabelSearch -help
SiprosV4OMP -h
V4Scripts_sipros_peptides_filtering -h
V4Scripts_sipros_peptides_assembling -h
V4Scripts_ClusterSip -h
V4Scripts_getLabelPCTinEachFT -help
```

### Install environment by yourself

Raxport relies on .net. Some scripts rely on python2 and R.

```bash
conda create -n py2 scikit-learn python=2.7
conda create -n mono -c conda-forge mono
conda create -n r -c conda-forge -c bioconda r-base r-stringr r-tidyr bioconductor-biostrings
```

### Convert Raw files

```bash
conda activate mono
# -j is the threads that you want to limit
mono bin/Raxport.exe -i raw -o ft -j 8
```

### Unlabeled search demo

One-key Script is in [EnsembleScripts](EnsembleScripts/cmd.sh)

Slurm Script is in [EnsembleScripts](EnsembleScripts/UnlabelForSlurm.sb)

```bash
# OMP_NUM_THREADS is the threads that you want to limit
export OMP_NUM_THREADS=10
# search the scans against the fasta databse, this command will take a long time
bin/SiprosEnsembleOMP -f demo.FT2 -c configTemplates/SiprosEnsembleConfig.cfg -o sip

conda activate py2
# convert .Spe2Pep.txt file to .tab file
EnsembleScripts/sipros_psm_tabulating.py -i sip -c configTemplates/SiprosEnsembleConfig.cfg -o sip
# filter PSMs, output qualified PSMs to .psm.txt file
EnsembleScripts/sipros_ensemble_filtering.py -i demo.tab configTemplates/SiprosEnsembleConfig.cfg -o sip
# assembly protein groups from peptide, output proteins to .pro.txt
EnsembleScripts/sipros_peptides_assembling.py -c configTemplates/SiprosEnsembleConfig.cfg -w sip

conda activate r
# control FDR, output qualified protein groups to .proRefineFDR.txt
Rscript V4Scripts/refineProteinFDR.R -pro demo.pro.txt -psm demo.psm.txt -fdr 0.005 -o demo
# get spectra count of each protein groups, output spectra count to .SPcount.txt
Rscript V4Scripts/getSpectraCountInEachFT.R -pro dmo.proRefineFDR.txt -psm demo.psm.txt -o demo
```

### Labeled search demo

One-key Script is in [V4Scripts](V4Scripts/SIPcmd.sh)

Slurm Script is in [V4Scripts](V4Scripts/LabelForSlurm.sb)

```bash
# generate configs
mkdir configs raw ft sip fasta
configGenerator -i SiprosV4Config.cfg -o configs -e C

conda activate r

# make db of identified proteins by SiprosEnsemble
Rscript V4Scripts/makeDBforLabelSearch.R \
    -pro YOUR_PATH/*.pro.txt \
    -faa YOUR_PATH/db.fasta \
    -o fasta/db.faa

# search the scans against the fasta database, this command will take a long time
SiprosV4OMP -f demo.FT2 -c C13_10000Pct.cfg -o sip

conda activate py2

# filter PSMs
python V4Scripts/sipros_peptides_filtering.py \
    -c SiprosV4Config.cfg \
    -w sip

# filter proteins
python V4Scripts/sipros_peptides_assembling.py \
    -c SiprosV4Config.cfg \
    -w sip

# cluster SIP abundance of protein
python V4Scripts/ClusterSip.py \
    -c SiprosV4Config.cfg \
    -w sip

conda activate r

# refine protein FDR
Rscript ${binPath}/refineProteinFDR.R \
    -pro sip/*.pro.txt \
    -psm sip/*.psm.txt \
    -fdr 0.01 \
    -o sip/YOUR_FILE_NAME

# get SIP abundance of each protein in each FT2 file
Rscript ${binPath}/getLabelPCTinEachFT.R \
    -pro sip/YOUR_FILE_NAME.proRefineFDR.txt \
    -psm sip/*.psm.txt \
    -thr 5 \
    -o sip/YOUR_FILE_NAME
```

### Compile this project

If you want to use mpi version on your server or make the binary suitable for CPU in other architecture, you can compile it by yourself. The source code is in [configGenerator](./configGenerator/), [siprosEnsembleCmakeAll](./siprosEnsembleCmakeAll/), and [SiprosV4CmakeAll](./SiprosV4CmakeAll/). They are all Cmake project and easy to compile in any IDE.



