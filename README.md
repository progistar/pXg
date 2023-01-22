# pXg: proteomics X genomics
---
- [About pXg](#about-pxg)
- [Usage](#usage)
  - [Input format](#input-format)
  - [Output format](#output-format)
  - [Command-line interface](#command-line-interface)
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
| :---:   | :------:       | :---:     | :---:       |
| Searh result       | A list of PSMs identified from a search engine (e.g. PEAKS, pNovo3*)     | TSV or CSV | Yes   |
| Gene information   | It must be the same file used in the read alignment (e.g. STAR2)         | GTF        | Yes   |
| RNA-Seq reads      | Mapped and unmapped RNA-Seq reads                                        | SAM        | Yes   |
| Protein sequences  | Canonical and contaminant protein sequences (e.g. UniProt)               | Fasta      | No    |

*pXg is not applicable to the flat formatted output in pNovo3. A user must convert the flat format to CSV or TSV.<br>
 We provide a convertor for pNovo3 v3.13 (see "Command-line interface").
### Output format
|Output    | Description    | Format   | Mandatory   |
| :---:   | :------:       | :---:     | :---:       |
| pXg result              | This is a main output file and contains a list of identification         | TSV         | Yes   |
| Read distributions      | Distributions for peptides per target/decoy read count                   | Flat        | Yes   |
| PSM score distributions | Target and decoy PSM score distributions                                 | Flat        | Yes   |
| Unmapped-reads match    | A list of unmapped reads matching to peptides                            | Flat        | Yes   |
| Matched reads*          | Matched reads to peptides passing all filters                            | SAM         | No    |
| Matched peptides*       | A list of peptides passing all filters                                   | GTF         | No    |

*Although the main output, pXg result contains PSM information with corresponding RNA-Seq counts, it is not suitable for visualization. <br>
 Two output files (matched reads and peptides) are available for direct use in IGV, making visualization easier.
### Command-line interface

## 3rd-party application
### IGV viewer
### Percolator

