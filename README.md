# DiffInvex: Computational framework for quantifying changes in somatic selection during cancer evolution and treatment. 

**DiffInvex (Differential Introns Versus Exons)** is a statistical method to identify conditional selection on point mutations between two or more time points/conditions. **DiffInvex** first estimates the local baseline mutation rate (BMR) based on the intronic mutations in addition to UTRs and flanking intergenic mutations, controlling for the heterogeneous mutational landscape across the genome, without the need for inferring BMR from covariate information such as gene expression and replication time. Additionally, **DiffInvex** employs a locus sampling approach  to control for within-gene variation in mutation risk by matching the trinucleotide and pentanucleotide composition as well as the DNA methylation status between target (exonic) and background (intronic and intergenic) regions. Next, **DiffInvex** utilized a Poisson regression model, further regularized by a weakly-informative prior, to determine the differential excess of point mutations in target regions over the baseline regions. This regression model can control for confounding factors between conditions such as tumor types and for technical variation between data sourced from different cohorts.

**For a full description of the method and applications, please visit [DiffInvex Manuscript](https://www.biorxiv.org/content/10.1101/2024.06.17.599362v1).**
  
## Contents
- [Download](#Download)
- [Directory Setup](#directory_setup)
- [Annotations](#annotations)
- [Input Preparation](#input_preparation)
- [Parameters](#parameters)
- [Usage](#usage)
  
     
### <a name="Download"></a>Download
```bash
cd ~
git clone https://github.com/AISKhalil/HiCNAtra.git
```
   
     
### <a name="directory_setup"></a>Directory Setup
After downloading the **HiCNAtra** directory, all the annotation files (reference genome sequence, mappability track, and GC track) should be downloaded and allocated to their corresponding sub-directories inside the **HiCNAtraTool** directory:
- The annotations directory structure will look like this:

```
    HiCNAtraTool/
    +- @HiCNAtra/
    +- Annotations/
    |  +- hg19/
    |  |  +- UCSC_chromFa/
    |  |  |  +- chr1.fa
    |  |  |  +- chr2.fa
    |  |  |  +- . . .
    |  |  |
    |  |  +- Anshul_UniqueMappability/
    |  |  |  +- globalmap_k20tok81/
    |  |  |  |  +- chr1.uint8.unique
    |  |  |  |  +- chr2.uint8.unique    
    |  |  |  |  +- . . .
    |  |  |  |
    |  |  |  +- globalmap_k101tok101/    
    |  |  |  +- . . .    
    |  |  |
    |  |  +- ChrisaMiller_GCContents/
    |  |  |  +- gcWinds.readLength100/
    |  |  |  |  +- chr1.gc
    |  |  |  |  +- chr2.gc
    |  |  |  |  +- . . .    
    |  |  |  |
    |  |  |  +- gcWinds.readLength50/
    |  |  |  +- . . .    
    |  |  +- UCSC_Centromeres.txt
    |  |  +- UCSC_Telomeres.txt
    |  |  +- UCSC_gapRegions.txt
    |  |  +- Anshul_wgEncodeHg19ConsensusSignalArtifactRegions.bed
    |  |  
    |  +- hg18/
    |  +- hg38/
    |  +- . . .
```
  
  
### <a name="annotations"></a>Annotations
- **HiCNAtra** uses the reference genome sequence (e.g. hg19) for computing the effective (restriction fragment) lengths that are used for correcting the Hi-C/3C-seq contact map based on the experiment (choice of restriction enzyme).
Please download and extract the reference genome sequence [hg19 reference genome sequence](http://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/chromFa.tar.gz) in `...HiCNAtra/HiCNAtraTool/Annotations/hg19/UCSC_chromFa/` sub-directory. 
  
- **HiCNAtra** also uses the unique mappability tracks for computing the mappability scores that are used for correcting the Hi-C/3C-seq contact map  [Unique mappability tracks for several species](https://sites.google.com/site/anshulkundaje/projects/mappability). This includes per-base unique mappability tracks for a large range of read lengths for several key species [Umap and Bismap: quantifying genome and methylome mappability](https://academic.oup.com/nar/article/46/20/e120/5086676). Please download and extract the mappability tracks [globalmap_k101tok101](https://personal.broadinstitute.org/anshul/projects/umap/encodeHg19Male/globalmap_k101tok101.tgz) and [globalmap_k20tok81](https://personal.broadinstitute.org/anshul/projects/umap/encodeHg19Female/globalmap_k20tok81.tgz) in `...HiCNAtra/HiCNAtraTool/Annotations/hg19/Anshul_UniqueMappability/` sub-directory.

- Additionally, **HiCNAtra** computes the GC score, for correcting the Hi-C/3C-seq contact maps, from the reference genome sequence (e.g. hg19). Alteratively, **HiCNAtra** has the option to use the Chris Miller's pre-calculated tracks [GC tracks](https://xfer.genome.wustl.edu/gxfer1/project/cancer-genomics/readDepth/index.html). This includes the GC tracks for a large range of read lengths for human genome [ReadDepth](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0016327). Please download and extract the GC tracks based on the read length [gcWinds.readLength100.hg19](https://xfer.genome.wustl.edu/gxfer1/project/cancer-genomics/readDepth/gcWinds.readLength100.hg19.tar), [gcWinds.readLength200.hg19](https://xfer.genome.wustl.edu/gxfer1/project/cancer-genomics/readDepth/gcWinds.readLength200.hg19.tar), [gcWinds.readLength76.hg19](https://xfer.genome.wustl.edu/gxfer1/project/cancer-genomics/readDepth/gcWinds.readLength76.hg19.tar), [gcWinds.readLength50.hg19](https://xfer.genome.wustl.edu/gxfer1/project/cancer-genomics/readDepth/gcWinds.readLength50.hg19.tar), [gcWinds.readLength36.hg19](https://xfer.genome.wustl.edu/gxfer1/project/cancer-genomics/readDepth/gcWinds.readLength36.hg19.tar), and [gcWinds.readLength27.hg19](https://xfer.genome.wustl.edu/gxfer1/project/cancer-genomics/readDepth/gcWinds.readLength27.hg19.tar) in `...HiCNAtra/HiCNAtraTool/Annotations/hg19/ChrisaMiller_GCContents/` sub-directory. 
User can select the GC calculation method by setting the `gcCalculationMethod` parameter.
  
  
### <a name="input_preparation"></a>Input Preparation
**HiCNAtra** input is HDF5 files that are generated by [hiclib](https://mirnylab.bitbucket.io/hiclib/index.html?) after applying the [iterative mapping](https://mirnylab.bitbucket.io/hiclib/tutorial/01_iterative_mapping.html) module only. They include the information needed for **HiCNAtra** in a dict-like structure with the main keys `'chrms1', 'chrms2', 'cuts1', 'cuts2', 'rfragIdxs1', 'rfragIdxs2', 'strands1', 'strands2','rsites1','rsites2'`.

  - **(1)** install the [hiclib](https://mirnylab.bitbucket.io/hiclib/index.html?).  
  - **(2)** edit the iterative mapping module based on the read length and restriction enzyme information [Mapping.py](./Scripts/Mapping.py)  
    ` python ./Scripts/Mapping.py -i inputFastQFilePath -o outputFileName`.
  - **(3)** extract the main keys of Hi-C data and convert from `h5dict` to `HDF5` format that can be read by HiCNAtra using [h5dictToHDF5.py](./Scripts/h5dictToHDF5.py)  
    ` python ./Scripts/h5dictToHDF5.py -i inputDictFile -o outputDictFile`.
  
  
### <a name="parameters"></a>Parameters
The main analysis parameters of **HiCNAtra**:

    HDF5Files                - the input HDF5 files {'/input/file1.hdf5','/input/file2.hdf5', ..}. 
                               This HDF5 files are the output of hiclib tool after running the iterative aligning module.
   
    readLength               - the short sequencing read length.
 
    restrictionEnzyme        - the name of the restriction enzyme that is used for the Hi-C/3C-seq experiment.

    maximumMoleculeLength    - the maximum molecule length (in bps). 

    referenceGenome          - the reference genome (e.g. hg19).
      
    binSize                  - the bin size of the RD signal (default = 5Kb).

    contactMapBinSize        - the bin size of the contact map (default = 100Kb).

    outputDirectory          - the directory that is used for saving all CNV information, raw contact map, and corrected contact map.

    RDmethod                 - the method to be used for computing the RD signal (default = 1): (1) "entire restriction fragment"
                               counting (for Hi-C data), (2) Paired-end method (for 3C-seq), (3) Exact-cut (position) approach,
                               (4) Midpoint (of restriction fragment) approach.
    
    gcCalculationMethod      - the method to be used for computing the GC scores (default = 2): (1) from Christopher A. Miller's 
                               pre-calculated GC tracks, (2) from the reference genome sequence.
                               
    ploidyLevel              - the whole-genome ploidy level of the cell line {'diploid', 'triploid', 'tetraploid', 'free'(default)}.
    
    cisOnly                  - a flag to compute and normalize the cis interaction frequencies only (cisOnly = 1), or both cis and trans 
                               interaction frequencies (cisOnly = 0).
 
   
     
### <a name="usage"></a>Usage 
Here, we use **GM06990** small sample as an example [GSM455133](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM455133). The input HDF5 file [GM06990_SRR027956_Input.hdf5](Example/GM06990_SRR027956_Input.hdf5) is provided in `Example/` sub-directory. Output files are saved in `Example/GM06990_HiC` sub-directory.   
   
**N.B.** this sample has only ~5 million reads, we use it for testing the HiCNAtra installation. For accurate analysis, you can incorporate more biological replicates [GSE18199](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE18199).
   
Start Matlab, then edit and run the following set of commands based on your data [runHiCNAtraScript.m](./Scripts/runHiCNAtraScript.m).
```
% Define the required information for creating HiCNAtra object based on the experiment and reference genome 
restrictionEnzyme = 'hindIII';
maximumMoleculeLength = 500;
readLength = 76;
referenceGenome = 'hg19';

% Add 'HiCNAtraTool' directory to Matlab search path
HiCNAtraDirectory = 'HiCNAtraTool';
addpath(HiCNAtraDirectory);

% Define the input HDF5 file(s)
HDF5Files = {'Example/GM06990_SRR027956_Input.hdf5'};

% Create HiCNAtra object 'GM06990_HiC' with the defined parameters
GM06990_HiC = HiCNAtra(HDF5Files, HiCNAtraDirectory, readLength, restrictionEnzyme, maximumMoleculeLength, referenceGenome);

% Set more parameters (optional)
GM06990_HiC.contactMapBinSize = 500000;
GM06990_HiC.ploidyLevel = 'diploid';
GM06990_HiC.outputDirectory = 'Example/GM06690_HiC';

% run 'RD calculator' module (Pipeline stage 1)
GM06990_HiC.RDcalculator;

% run 'CNV caller' module (Pipeline stage 2)
GM06990_HiC.ploidyLevel = 'diploid';
GM06990_HiC.CNVcaller;

% run 'contact map corrector' (Pipeline stage 3) module that compute and correct the contact map
GM06990_HiC.contactMapCorrector;

% save the HiCNAtraObject, so you can load it directly for further analysis.
save('Example/GM06990_HiC.mat');

% plot the CNV tracks (e.g chr11)
chrNumber = 11;
GM06990_HiC.CNVsTrackPlot('plot',chrNumber);

% plot the raw contact map (e.g. chr1 )
chrNumber = 1;
GM06990_HiC.rawContactMapPlot('plot',chrNumber);

% plot the HiCNAtra-corrected contact map (e.g. chr1 )
chrNumber = 1;
GM06990_HiC.normContactMapPlot('plot',chrNumber); 

```

## HiCNAtra user manual . . . Coming soon!  
