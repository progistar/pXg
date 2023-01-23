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
    - [Match-only command](#match-only-command)
- [3rd-party application](#3rd-party-application)
  - [IGV viewer](#igv-viewer)
  - [Percolator](#percolator)
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
| Searh result       | A list of PSMs identified from a search engine (e.g. PEAKS, pNovo3*)     | TSV or CSV | Yes   |
| Gene information   | It must be the same file used in the read alignment (e.g. STAR2)         | GTF        | Yes   |
| RNA-Seq reads      | Mapped and unmapped RNA-Seq reads                                        | SAM        | Yes   |
| Protein sequences  | Canonical and contaminant protein sequences (e.g. UniProt)               | Fasta      | No    |

*pXg is not applicable to the flat formatted output in pNovo3. A user must convert the flat format to CSV or TSV.<br>
 We provide a convertor for pNovo3 v3.13 (see "Command-line interface").
### Output format
|Output    | Description    | Format   | Mandatory   |
| :---:   | :---:       | :---:     | :---:       |
| pXg result              | This is a main output file and contains a list of identification         | TSV         | Yes   |
| Read distributions      | Distributions for peptides per target/decoy read count                   | Flat        | Yes   |
| PSM score distributions | Target and decoy PSM score distributions                                 | Flat        | Yes   |
| Unmapped-reads match    | A list of unmapped reads matching to peptides                            | Flat        | Yes   |
| Matched reads*          | Matched reads to peptides passing all filters                            | SAM         | No    |
| Matched peptides*       | A list of peptides passing all filters                                   | GTF         | No    |

*Although the main output, pXg result contains PSM information with corresponding RNA-Seq counts, it is not suitable for visualization. <br>
 Two output files (matched reads and peptides) are available for direct use in IGV, making visualization easier.
### Command-line interface
#### List of Parameters
|Option    | Description    | Mandatory   |
| :---:   | :---:       | :---:     |
| gtf       | GTF file path. We recommand to use the same gtf corresponding to alignment |Yes   |
| sam       | SAM file path. The sam file must be sorted by coordinate |Yes   |
| psm       | PSM file path. It is expected that the psm file is derived from proteomics search by de novo or database search engine |Yes   |
| pept_col       | Peptide column index in the psm file |Yes   |
| scan_cols       | Scan identifier indices in the psm file. Multiple columns are also possible because sometimes it is not enough to distinguish scan by only scan id. You can write multiple indices such like that: 1,2,5 |Yes   |
| out       | Output path of pXg |Yes   |
| sep       | Specify the column separator. Possible values are csv or tsv. Default is csv |No   |
| pval       | p-value cutoff of randomly matched peptide-read pairs. Default is 0.01 |No   |
| fdr       | FDR cutoff to discard low-quality peptide-spectrum matches. Default is 0.1 |No   |
| ileq       | Controls whether pXg treats isoleucine (I) and leucine (L) as the same/equivalent with respect to a peptide identification. Default is true |No   |
| length       | Range of peptide length to consider. Default is 8-15. You can write in this way (min-max, both inclusive) : 8-13 |No   |
| fasta       | Canonical sequence database to report conservative assignment of noncanonical PSMs |No   |
| rank       | How many candidates will be considered per a scan. Default is 100 (in other words, use all ranked candidates) |No   |
| out_sam       | Report matched reads as SAM format (true or false). Default is true |No   |
| out_gtf       | Report matched peptides as GTF format (true or false). Default is true |No   |
| out_canonical       | Report caonical peptides for SAM and/or GTF formats (true or false). Default is true |No   |
| out_noncanonical       | Report noncaonical peptides for SAM and/or GTF formats (true or false). Default is true |No   |
| pMut       | Penalty per a mutation. Default is 1 |No   |
| pAS       | Penalty for alternative splicing. Default is 10 |No   |
| p5UTR       | Penalty for 5`-UTR. Default is 20 |No   |
| p3UTR       | Penalty for 3`-UTR. Default is 20 |No   |
| pncRNA       | Penalty for noncoding RNA. Default is 20 |No   |
| pFS       | Penalty for frame shift. Default is 20 |No   |
| pIR       | Penalty for intron region. Default is 30 |No   |
| pIGR       | Penalty for intergenic region. Default is 30 |No   |
| pasRNA       | Penalty for antisense RNA. Default is 30 |No   |
| punmap       | Penalty for unmapped reads. Default is 100 |No   |
| gtf_partition_size*       | The size of treating genomic region at once. Default is 5000000 |No   |
| sam_partition_size*       | The size of treating number of reads at once. Default is 1000000 |No   |

*size parameters can effect memory usage and time. If your machine does not have enough memory, then decrease those values.

#### Basic command
```bash
java -Xmx30G -jar pXg.jar -gtf gencode.gtf -sam aligned.sorted.sam -psm peaks.result -fasta uniprot_contam.fasta -pept_col 4 -score_col 8 -scan_cols 1,2,5  -pval 0.01 -fdr 0.1 -out_canonical false -out peaks.pXg
```

#### Match-only command
If one wants to map the identified peptides to RNA-Seq reads, we recommand to use below command line.
```bash
java -Xmx30G -jar pXg.jar -gtf gencode.gtf -sam aligned.sorted.sam -psm peaks.result -fasta uniprot_contam.fasta -pept_col 4 -score_col 8 -scan_cols 1,2,5 -ileq false  -pval 0.01 -fdr 1 -out peaks.pXg
```

## 3rd-party application
### IGV viewer

### Percolator

