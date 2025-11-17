### SiprosToolKits

SiprosToolKits hosts the utilities developed by the Sipros team for stable-isotope mass spectrometry-based metaproteomics. The toolbox includes Raxport (RAW file conversion), SiprosV4, SiprosEnsemble, and the automation scripts that connect them into complete unlabeled or SIP-labeled workflows.

#### Sipros 5 (recommended)

Sipros 5 delivers automatic SIP-labeling FDR control, improved identification sensitivity (≈3× more labeled identifications than Sipros 4), and DIA support. Install it directly from Bioconda:

```bash
conda install bioconda::sipros
siproswf -h
sipros -h
```

- [Sipros 5 repository and tutorial](https://github.com/thepanlab/sipros5)  
- [Sipros package on Bioconda](https://anaconda.org/bioconda/sipros)

#### Sipros 4 (legacy release)

Sipros 4 remains available for backward compatibility or when you need the original processing steps. Create the legacy environment with the pinned Bioconda build:

```bash
conda create -n sipros4 bioconda::sipros=4.02
conda activate sipros4

# Useful entry points (all ship with -h / -help documentation)
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

A concise tutorial for a 13C-labeled *E. coli* dataset is available on the [project wiki](https://github.com/thepanlab/SiprosToolKits/wiki/13C-labeled-E.-coli-SIP-proteomic-search-tutorial).

#### Supporting environments

Raxport requires Mono/.NET, several scripts rely on Python 2, and the downstream statistics make use of R. Set up the helper environments once and activate them as needed:

```bash
conda create -n py2 scikit-learn python=2.7
conda create -n mono -c conda-forge mono
conda create -n r -c conda-forge -c bioconda r-base r-stringr r-tidyr bioconductor-biostrings
```

#### Workflow overview

##### Convert vendor RAW files

```bash
conda activate mono
# -j limits Raxport to the desired number of threads
mono bin/Raxport.exe -i raw -o ft -j 8
```

##### Unlabeled search demo

Use the one-key script in `EnsembleScripts/cmd.sh` or the Slurm-friendly submission script in `EnsembleScripts/UnlabelForSlurm.sb`. To run the individual steps manually:

```bash
export OMP_NUM_THREADS=10

# 1. Search against your FASTA database (long-running)
bin/SiprosEnsembleOMP \
    -f demo.FT2 \
    -c configTemplates/SiprosEnsembleConfig.cfg \
    -o sip

conda activate py2

# 2. Convert *.Spe2Pep.txt files to tabular output
EnsembleScripts/sipros_psm_tabulating.py \
    -i sip \
    -c configTemplates/SiprosEnsembleConfig.cfg \
    -o sip

# 3. Filter PSMs and emit demo.psm.txt
EnsembleScripts/sipros_ensemble_filtering.py \
    -i demo.tab \
    -c configTemplates/SiprosEnsembleConfig.cfg \
    -o sip

# 4. Assemble peptide-level results into proteins
EnsembleScripts/sipros_peptides_assembling.py \
    -c configTemplates/SiprosEnsembleConfig.cfg \
    -w sip

conda activate r

# 5. Control protein-group FDR
Rscript V4Scripts/refineProteinFDR.R \
    -pro demo.pro.txt \
    -psm demo.psm.txt \
    -fdr 0.005 \
    -o demo

# 6. Summarize spectra counts per protein group
Rscript V4Scripts/getSpectraCountInEachFT.R \
    -pro demo.proRefineFDR.txt \
    -psm demo.psm.txt \
    -o demo
```

##### Labeled search demo

Use `V4Scripts/SIPcmd.sh` for an automated local run or `V4Scripts/LabelForSlurm.sb` for cluster execution. The manual sequence is outlined below:

```bash
# 1. Prepare directories and templates
mkdir -p configs raw ft sip fasta
configGenerator -i SiprosV4Config.cfg -o configs -e C

conda activate r

# 2. Build a reduced FASTA based on SiprosEnsemble identifications
Rscript V4Scripts/makeDBforLabelSearch.R \
    -pro YOUR_PATH/*.pro.txt \
    -faa YOUR_PATH/db.fasta \
    -o fasta/db.faa

# 3. Search the scans against the labeled database (long-running)
SiprosV4OMP \
    -f demo.FT2 \
    -c C13_10000Pct.cfg \
    -o sip

conda activate py2

# 4. Filter peptide-spectrum matches
python V4Scripts/sipros_peptides_filtering.py \
    -c SiprosV4Config.cfg \
    -w sip

# 5. Assemble proteins from filtered peptides
python V4Scripts/sipros_peptides_assembling.py \
    -c SiprosV4Config.cfg \
    -w sip

# 6. Cluster SIP abundance per protein group
python V4Scripts/ClusterSip.py \
    -c SiprosV4Config.cfg \
    -w sip

conda activate r

# 7. Refine protein FDR
Rscript ${binPath}/refineProteinFDR.R \
    -pro sip/*.pro.txt \
    -psm sip/*.psm.txt \
    -fdr 0.01 \
    -o sip/YOUR_FILE_NAME

# 8. Quantify SIP abundance for every FT2 file
Rscript ${binPath}/getLabelPCTinEachFT.R \
    -pro sip/YOUR_FILE_NAME.proRefineFDR.txt \
    -psm sip/*.psm.txt \
    -thr 5 \
    -o sip/YOUR_FILE_NAME
```

#### Build from source

If you need MPI-enabled binaries or want to target specific CPU architectures, compile the projects in `configGenerator/`, `siprosEnsembleCmakeAll/`, and `siprosV4CmakeAll/`. Each directory is a CMake project that can be configured in your IDE or via standard CMake workflows.

#### Citation

Xiong, Yi, Ryan S. Mueller, Shichao Feng, Xuan Guo, and Chongle Pan. “Proteomic stable isotope probing with an upgraded Sipros algorithm for improved identification and quantification of isotopically labeled proteins.” *Microbiome* 12 (2024).
