### SiprosToolKits

Tools collection of Sipros for stable isotopic mass spectrum meta proteomic research. Raxport, SiprosV3, SiprosEnsemble and some scripts are in this repository.

### Install environment 

```bash
conda create -n py2 scikit-learn python=2.7
conda create -n mono -c conda-forge mono
conda create -n r -c conda-forge -c bioconda r-base r-stringr r-tidyr bioconductor-biostrings
```

### Make temp folder

```bash
cd demo
mkdir fasta raw ft sip configs 
cd ..
```

copy all .raw files to raw

### Convert Raw files

```bash
conda activate mono
# -j is the threads that you want to limit
mono bin/Raxport.exe -i raw -o ft -j 8
```

### Unlabeled search

```bash
# OMP_NUM_THREADS is the threads that you want to limit
export OMP_NUM_THREADS=10
# search the scans against the fasta databse, this command will take a long time
bin/SiprosEnsembleOMP -f demo/ft/demo.FT2 -c configTemplates/SiprosEnsembleConfig.cfg -o demo/sip

conda activate py2
# convert .Spe2Pep.txt file to .tab file
EnsembleScripts/sipros_psm_tabulating.py -i demo/sip -c configTemplates/SiprosEnsembleConfig.cfg -o demo/sip
# filter PSMs, output qualified PSMs to .psm.txt file
EnsembleScripts/sipros_ensemble_filtering.py -i demo.tab configTemplates/SiprosEnsembleConfig.cfg -o demo/sip
# assembly protein groups from peptide, output proteins to .pro.txt
EnsembleScripts/sipros_peptides_assembling.py -c configTemplates/SiprosEnsembleConfig.cfg -w demo/sip

conda activate r
# control FDR, output qualified protein groups to .proRefineFDR.txt
Rscript V3Scripts/refineProteinFDR.R -pro demo.pro.txt -psm demo.psm.txt -fdr 0.005 -o demo
# get spectra count of each protein groups, output spectra count to .SPcount.txt
Rscript V3Scripts/getSpectraCountInEachFT.R -pro dmo.proRefineFDR.txt -psm demo.psm.txt -o demo
```
One-key Script is in [EnsembleScripts](EnsembleScripts/cmd.sh)

Slurm Script is in [EnsembleScripts](EnsembleScripts/UnlabelForSlurm.sb)

### Labeled search

```bash
# generate configs
mkdir configs raw ft sip fasta
configGenerator -i SiprosV3Config.cfg -o configs -e C

conda activate r

# make db of identified proteins by SiprosEnsemble
Rscript V3Scripts/makeDBforLabelSearch.R \
    -pro YOUR_PATH/*.pro.txt \
    -faa YOUR_PATH/db.fasta \
    -o fasta/db.faa

conda activate py2

# filter PSMs
python V3Scripts/sipros_peptides_filtering.py \
    -c SiprosV3Config.cfg \
    -w sip

# filter proteins
python V3Scripts/sipros_peptides_assembling.py \
    -c SiprosV3Config.cfg \
    -w sip

# cluster SIP abundance of protein
python V3Scripts/ClusterSip.py \
    -c SiprosV3Config.cfg \
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

One-key Script is in [V3Scripts](V3Scripts/SIPcmd.sh)

Slurm Script is in [V3Scripts](V3Scripts/LabelForSlurm.sb)

### Compile this project

If you want to use mpi version on your server or make the binary suitable for CPU in other architecture, you can compile it by yourself. The source code is in [configGenerator](./configGenerator/), [siprosEnsembleCmakeAll](./siprosEnsembleCmakeAll/), and [siprosV3CmakeAll](./siprosV3CmakeAll/). They are all Cmake project and easy to compile in any IDE.



