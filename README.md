# pXg: proteomics X genomics
<img src="https://github.com/progistar/pXg/blob/main/img/intro.png"/>

---
- [About pXg](#about-pxg)
- [Usage](#usage)
  - [Input format](#input-format)
  - [Output format](#output-format)
  - [Command-line interface](#command-line-interface)
    - [List of parameters](#list-of-parameters)
    - [Basic command](#basic-command)
- [Tutorial](#tutorial)
  - [RNA-Seq alignment](#rna-seq-alignment)
  - [SAM preparation](#sam-preparation)
  - [Toy example](#toy-example)
  - [Run pXg](#run-pxg)
  - [Run Percolator using the pXg results](#run-percolator-using-the-pxg-results)
  - [IGV viewer](#igv-viewer)
- [TODO](#todo)
  - [Integrative pipeline](#integrative-pipeline)
---

## About pXg

pXg (proteomics X genomics), a software tool that enables the reliable identification of both cMAPs and ncMAPs from de novo peptide sequencing by utilizing RNA-Seq data.
<br>

## Usage
pXg can be integrated with any search engines such as PEAKS and pNovo3.
It was developed for the reliable identification of ncMAPs from de novo peptide sequencing; however, it can also be used to capture the number of reads mapped to a peptide sequence.
### Input format
|Input    | Description    | Format    | Mandatory   |
| :---:   | :---:       | :---:     | :---:       |
| Searh result       | A list of PSMs identified from a search engine (e.g. <a href="https://www.bioinfor.com/peaks-studio/" target="_blank">PEAKS</a>, <a href="http://pfind.org/software/pNovo/index.html" target="_blank">pNovo3</a>, <a href="https://github.com/Noble-Lab/casanovo" target="_blank">Casanovo</a>)     | TSV or CSV | Yes   |
| Gene annotation    | It must be the same file used in the read alignment (e.g. <a href="https://www.gencodegenes.org/" target="_blank">Gencode</a>, <a href="https://ensemblgenomes.org/" target="_blank">Ensembl</a>)       | GTF        | Yes   |
| RNA-Seq reads      | Mapped and unmapped RNA-Seq reads. The file must be sorted by coordinates | SAM        | Yes   |
| Protein sequences  | Canonical and contaminant protein sequences (e.g. UniProt)                        | Fasta      | No    |

*pXg is not applicable to the flat formatted output in pNovo3. A user must convert the flat format to CSV or TSV.<br>

### Output format
|Output    | Description    | Format   | Mandatory   |
| :---:   | :---:       | :---:     | :---:       |
| pXg result                | This is a main output file and contains a list of identification as TSV format         | TSV         | Yes   |
| pXg result for Percolator | This is a main output file and contains a list of identification as PIN format         | PIN         | Yes   |
| Unknown sequences          | A list of unmapped reads matching to peptides                            | Flat        | Yes   |
| Matched reads*             | Matched reads to peptides passing all filters                            | SAM         | No    |
| Matched peptides*          | A list of peptides passing all filters                                   | GTF         | No    |

*Although the pXg result contains PSM information with corresponding RNA-Seq counts, it is not suitable for visualization. <br>
 Two output files (matched reads and peptides) are available for direct use in <a href="https://software.broadinstitute.org/software/igv/" target="_blank">IGV</a>, making visualization easier.

### Command-line interface
#### List of Parameters
|Option    | Description    | Mandatory   |
| :---:   | :---:       | :---:     |
| gtf_file       | GTF file path. We recommand to use the same gtf corresponding to alignment |Yes   |
| sam_file       | SAM file path. The sam file must be sorted by coordinate |Yes   |
| psm_file       | PSM file path. It is expected that the psm file is derived from proteomics search by de novo or database search engine |Yes   |
| pept_col       | Peptide column index in the psm file |Yes   |
| scan_col       | Scan number index in the psm file    |Yes   |
| output         | Base output name of pXg |Yes   |
| sep            | Specify the column separator. Possible values are csv or tsv. Default is csv |No   |
| mode           | Specify the method of translation nucleotides. 3 for three-frame and 6 for six-frame. Default is 6 |No   |
| add_feat_cols  | Specify the indices for additional features to generate PIN file. Several features can be added by comma separator. ex> 5,6,7|No  |
| ileq           | Controls whether pXg treats isoleucine (I) and leucine (L) as the same/equivalent with respect to a peptide identification. Default is true |No   |
| lengths        | Range of peptide length to consider. Default is 8-15. You can write in this way (min-max, both inclusive) : 8-13 |No   |
| fasta_file     | Canonical sequence database to report conservative assignment of noncanonical PSMs |No   |
| rank           | How many candidates will be considered per a scan. Default is 100 (in other words, use all ranked candidates) |No   |
| out_sam        | Report matched reads as SAM format (true or false). Default is true |No   |
| out_gtf        | Report matched peptides as GTF format (true or false). Default is true |No   |
| out_canonical  | Report caonical peptides for SAM and/or GTF formats (true or false). Default is true |No   |
| out_noncanonical| Report noncaonical peptides for SAM and/or GTF formats (true or false). Default is true |No   |
| penalty_mutation   | Penalty per a mutation. Default is 1 |No   |
| penalty_AS         | Penalty for alternative splicing. Default is 10 |No   |
| penalty_5UTR       | Penalty for 5`-UTR. Default is 20 |No   |
| penalty_3UTR       | Penalty for 3`-UTR. Default is 20 |No   |
| penalty_ncRNA      | Penalty for noncoding RNA. Default is 20 |No   |
| penalty_FS         | Penalty for frame shift. Default is 20 |No   |
| penalty_IR         | Penalty for intron region. Default is 30 |No   |
| penalty_IGR        | Penalty for intergenic region. Default is 30 |No   |
| penalty_asRNA      | Penalty for antisense RNA. Default is 30 |No   |
| penalty_unknown    | Penalty for unmapped reads. Default is 100 |No   |
| gtf_partition_size*       | The size of treating genomic region at once. Default is 5000000 |No   |
| sam_partition_size*       | The size of treating number of reads at once. Default is 1000000 |No   |
| threads*                  | The number of threads. Default is 4|No   |

*size parameters can effect memory usage and time. If your machine does not have enough memory, then decrease those values.

#### Basic command
```bash
java -Xmx30G -jar pXg.jar \
--gtf_file [gene annotation file path] \
--sam_file [sorted SAM file path] \
--psm_file [de novo result file path] \
--fasta_file [protein sequence fasta file paht] \
--pept_col [index of peptide column] \
--score_col [index of search score column] \
--scan_col [index of scan number column] \
--out_canonical false \
--output [base output file name]
```

## Tutorial
This tutorial aims to understand how to run pXg and estimate FDR from the result. It contains 1) running <a href="https://github.com/alexdobin/STAR" target="_blank">STAR2</a> aligner with 2-pass parameter, 2) preparing SAM file from the alignment, 3) running pXg and 4) several post-processing including <a href="https://github.com/percolator/percolator" target="_blank">Percolator</a>, merging pXg result with the result of Percolator and estimating separated FDR.
Note that it neither contains how to run de novo peptide sequencing engines such as <a href="https://www.bioinfor.com/peaks-studio/" target="_blank">PEAKS</a>, <a href="http://pfind.org/software/pNovo/index.html" target="_blank">pNovo3</a> and <a href="https://github.com/Noble-Lab/casanovo" target="_blank">Casanovo</a> AND how to create deep learning based features.

### RNA-Seq alignment
We recommand to align fastq files using STAR2 with The Cancer Genome Atlas (TCGA) <a href="https://docs.gdc.cancer.gov/Data/Bioinformatics_Pipelines/Expression_mRNA_Pipeline" target="_blank">two-pass alignment option</a>. 

### SAM preparation
Once you get the aligned BAM or SAM file, you MUST sort the file by chromosomal coordinates and convert to SAM.
Currently, pXg only supports SAM file.

We provide a code for preprocessing SAM file using <a href="http://www.htslib.org/" target="_blank">SAMtools</a> below:
```bash
samtools sort -o in.sorted.sam in.sam -@ 8
```

The "in.sorted.sam" is used for pXg input.

### Toy example
In this tutorial, toy datasets including 1) de novo results, 2) in.sorted.sam, 3) gene annotation (GTF) and 4) protein sequence fasta file are provided in the <a href="https://github.com/progistar/pXg/tree/main/tutorial" target="_blank">tutorial</a> folder so that a user can try to run the pXg pipeline. 

### Run pXg
Using the toy datasets, you can run the pXg pipline using following command: <br>
```bash
java -Xmx2G -jar pXg.v2.0.1.jar \
--gtf_file toy.gtf \
--sam_file toy.sorted.sam \
--psm_file toy.psm.csv \
--fasta_file toy.fasta \
--output toy \
--scan_col 5 \
--file_col 2 \
--pept_col 4 \
--score_col 8 \
--charge_col 11 \
--add_feat_cols 15 \
--sep csv \
--mode 6 \
--threads 2
```
This may take about 2 mins.<br>

Note that the memory option "-Xmx50G" depends on the size of SAM file. In our experience, "-Xmx30G" is enough to deal with ~20G file. 

### Run Percolator using the pXg results
Once you get the pXg result, you can add more features such as spectral similarity and delta retention time described in our manuscript. Without the additional features, still it is possible to run <a href="https://github.com/percolator/percolator" target="_blank">Percolator</a> and estimate FDR from the pXg results.<br>
We recommand to use Percolator version >= v3.06.1 because former versions have an issue to print proteinIds.<br>
Post processing codes are also provided in the tutorial folder (post_process.ipynb).

### IGV viewer
<img src="https://github.com/progistar/pXg/blob/main/img/toy.viewer.png"/>
When pXg finishes identifying peptides, the resulting GTF and SAM files are immediately available in the <a href="https://software.broadinstitute.org/software/igv/" target="_blank">IGV viewer</a>.

## TODO
### Integrative pipeline
* The pXg pipline will be intragted with SOTA de novo sequencing tools and RNA-Seq aligners to deal with 1) MS/MS spectra, 2) fastq (RNA-Seq), 3) Gene annotation, 4) Genome sequence, and 5) Protein sequence database at once!<br>
* For example, Casanovo is the SOTA now; therefore, we are implementing to deal with that of input file.