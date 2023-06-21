# bio-data
To process original .fasta file from PDB, generate the antibody's corresponding antigen .fasta and .pdb file

## 1. Download original fast and pdb files

You can use your own methods to download these files, usually using the download links from PDB and SAbDab. Then, store these files in under `\ori_fasta` and `\ori_pdb` folder.

## 2. Change the parameters in `config.py`

In the configuration file `config.py `, there are several parameters.

| Parameters  | Meaning                                               |
| ----------- | ----------------------------------------------------- |
| `file_path` | the path where you store data                         |
| `tsv_files` | the summary files you download by searching on SAbDab |
| `num_schm`  | numbering scheme                                      |
| `cdr_schm`  | the CDR definition scheme                             |

Note:

- If you use other `num_schm` and `cdr_schm`, you need to check and add/change the `process.py`. The following image is for your reference on `cdr_schm`.

## 3. Error check
From the two summray files generated, you can check for errors that may have occurred during processing.
