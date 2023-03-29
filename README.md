# MetaCC allows scalable and integrative analyses of both long-read and short-read metagenomic Hi-C data

- [Overview](#overview)
- [System Requirements](#system-requirements)
- [Installation Guide](#installation-guide)
- [A test dataset to demo MetaCC](#a-test-dataset-to-demo-metacc)
- [Instruction to process raw data](#instruction-to-process-raw-data)
- [MetaCC analysis](#metacc-analysis)
- [Instruction of reproducing results in MetaCC paper](https://github.com/dyxstat/Reproduce_MetaCC)
- [Contacts and bug reports](#contacts-and-bug-reports)
- [Copyright and License Information](#copyright-and-license-information)
- [Issues](https://github.com/dyxstat/MetaCC/issues)

# Overview
`MetaCC` is an efficient and integrative framework for analyzing both long-read and short-read metaHi-C datasets.
In the `MetaCC` framework, raw metagenomic Hi-C contacts are first efficiently and effectively normalized 
by a new normalization method, `NormCC`. Leveraging NormCC-normalized Hi-C contacts, 
the binning module in `MetaCC` enables the retrieval of high-quality MAGs and downstream analyses.


* **If you want to reproduce results in our MetaCC paper, please read our instructions [here](https://github.com/dyxstat/Reproduce_MetaCC).**

* **Scripts to process the intermediate data and plot figures of our MetaCC paper are available [here](https://github.com/dyxstat/Reproduce_MetaCC/tree/main/Scripts).**


# System Requirements
## Hardware requirements
`MetaCC` requires only a standard computer with enough RAM to support the in-memory operations.

## Software requirements
### OS Requirements
`MetaCC` v1.0.0 is supported and tested in *MacOS* and *Linux* systems.

### Python and R Dependencies
`MetaCC` mainly depends on the Python scientific stack:

```
numpy
scipy
pysam
scikit-learn
pandas
Biopython
igraph
leidenalg
```
and R scientific package:

```
MASS
```



# Installation Guide
We recommend using [**conda**](https://docs.conda.io/projects/conda/en/latest/user-guide/install/download.html) to install `MetaCC`. 
Typical installation time is 1-5 minutes depending on your system.

### Clone the repository with git
```
git clone https://github.com/dyxstat/MetaCC.git
```

Once complete, enter the repository folder and then create a `MetaCC` environment using conda.


### Enter the MetaCC folder
```
cd MetaCC
```

### Add dependencies of external softwares
Since MetaCC needs to execute the external softwares in the folder [Auxiliary](https://github.com/dyxstat/MetaCC/tree/main/Auxiliary), you can run the following commands to make sure that all external softwares are executable:
```
chmod +x Auxiliary/test_getmarker.pl
chmod +x Auxiliary/FragGeneScan/FragGeneScan
chmod +x Auxiliary/FragGeneScan/run_FragGeneScan.pl
chmod +x Auxiliary/hmmer-3.3.2/bin/hmmsearch
```

### Construct the conda environment in the linux or MacOS system
```
conda env create -f MetaCC_linux_env.yaml
```
or
```
conda env create -f MetaCC_osx_env.yaml
```

### Enter the conda environment
```
conda activate MetaCC_env
```

### Install the R package

```
# Enter the R
R

# Download the R package and you may need to select a CRAN mirror for the installation
install.packages('MASS')
```



# A test dataset to demo MetaCC
We provide a small dataset, located under the [Test](https://github.com/dyxstat/MetaCC/tree/main/Test) directory, to test the software:
```
python ./MetaCC.py test
```


# Instruction to process raw data
Follow the instructions in this section to process the raw shotgun and Hi-C data and generate the input for the `MetaCC` framework:

### Clean raw shotgun and Hi-C reads

Adaptor sequences are removed by `bbduk` from the `BBTools` suite with parameter `ktrim=r k=23 mink=11 hdist=1 minlen=50 tpe tbo` and reads are quality-trimmed using `bbduk` with parameters `trimq=10 qtrim=r ftm=5 minlen=50`. Additionally, the first 10 nucleotides of Hi-C reads are trimmed by `bbduk` with parameter `ftl=10`. Identical PCR optical and tile-edge duplicates for Hi-C reads were removed by the script `clumpify.sh` from `BBTools` suite.

### Assemble shotgun reads

For the shotgun library, de novo metagenome assembly is produced by an assembly software, such as MEGAHIT.
```
megahit -1 SG1.fastq.gz -2 SG2.fastq.gz -o ASSEMBLY --min-contig-len 1000 --k-min 21 --k-max 141 --k-step 12 --merge-level 20,0.95
```

### Align Hi-C paired-end reads to assembled contigs

Hi-C paired-end reads are aligned to assembled contigs using a DNA mapping software, such as BWA MEM. Then, samtools with parameters ‘view -F 0x904’ is applied to remove unmapped reads, supplementary alignments, and secondary alignments. BAM file needs to be sorted **by name** using 'samtools sort'.
```
bwa index final.contigs.fa
bwa mem -5SP final.contigs.fa hic_read1.fastq.gz hic_read2.fastq.gz > MAP.sam
samtools view -F 0x904 -bS MAP.sam > MAP_UNSORTED.bam
samtools sort -n MAP_UNSORTED.bam -o MAP_SORTED.bam
```


# MetaCC analysis
## Implement the NormCC normalization module
```
python ./MetaCC.py norm [Parameters] FASTA_file BAM_file OUTPUT_directory
```
### Parameters
```
-e (required): Case-sensitive enzyme name. Use multiple times for multiple enzymes 
--min-len: Minimum acceptable contig length (default 1000)
--min-mapq: Minimum acceptable alignment quality (default 30)
--min-match: Accepted alignments must be at least N matches (default 30)
--min-signal: Minimum acceptable signal (default 2)
--thres: the fraction of discarded normalized Hi-C contacts 
         (default 0.05, which means discarding the lowest 5% of normalized Hi-C contacts as spurious)
--cover (optional): Cover existing files. Otherwise, an error will be returned if the output file is detected to exist.
-v (optional): Verbose output about more specific details of the procedure.
```
### Input File

* **FASTA_file**: a fasta file of the assembled contig (e.g. Test/final.contigs.fa)
* **BAM_file**: a bam file of the Hi-C alignment (e.g. Test/MAP_SORTED.bam)


### Output File

* **contig_info.csv**: information of assembled contigs with three columns (contig name, the number of restriction sites on contigs, and contig length)
* **Normalized_contact_matrix.npz**: a sparse matrix of normalized Hi-C contact maps in csr format and can be reloaded using  **scipy.sparse.load_npz('Normalized_contact_matrix.npz')**
* **NormCC_normalized_contact.gz**: Compressed format of the normalized contacts and contig information by pickle. 
This file can further serve as the input of MetaCC binning module.
* **MetaCC.log**: the specific implementation information of NormCC normalization module


### Example
```
python ./MetaCC.py norm -v final.contigs.fa MAP_SORTED.bam out_directory
```


## Implement the MetaCC binning module
**MetaCC binning module is based on the NormCC-normalized Hi-C contacts and thus must be implemented after the NormCC normalization module.**
```
python ./MetaCC.py cluster [Parameters] FASTA_file OUTPUT_directory
```
### Parameters
```
--min-binsize: Minimum bin size used in output (default 150000)
--num-gene (optional): Number of marker genes detected. If there is no input, 
                       the number of marker genes will be automatically detected.
--random-seed (optional): seed for the Leiden clustering. If there is no input, a random seed will be employed.
--cover (optional): Cover existing files. Otherwise, an error will be returned if the output file is detected to exist.
-v (optional): Verbose output about more specific details of the procedure.
```
### Input File

* **FASTA_file**: a fasta file of the assembled contig (e.g. Test/final.contigs.fa)
* **OUTPUT_directory**: please pay attention that the output directory of the MetaCC binning module should be the same as that of the NormCC normalization module.

### Output File

* **BIN**: folder containing the fasta files of initial draft genomic bins
* **MetaCC.log**: the specific implementation information of NormCC normalization module




# Contacts and bug reports
If you have any questions or suggestions, welcome to contact Yuxuan Du (yuxuandu@usc.edu).


# Copyright and License Information
This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program. If not, see http://www.gnu.org/licenses/.







