# bio-data
To process original .fasta file from PDB, generate the antibody's corresponding antigen .fasta and .pdb file

## 1. Download original fast and pdb files

You can use your own methods to download these files, usually using the download links from PDB and SAbDab. Then, store these files in under `\ori_fasta` and `\ori_pdb` folder.

## 2. Change the parameters in `config.py`

In the configuration file `config.py` , there are several parameters.

| Parameters  | Meaning                                               |
| ----------- | ----------------------------------------------------- |
| `file_path` | the path where you store data                         |
| `tsv_files` | the summary files you download by searching on SAbDab |
| `num_schm`  | numbering scheme                                      |
| `cdr_schm`  | the CDR definition scheme                             |

Note:

- If you use other `num_schm` and `cdr_schm`, you need to check and add/change the `process.py`. The following image is for your reference on `cdr_schm`.

<img width="1185" alt="image" src="https://github.com/songjie-guo/bio-data/assets/69680257/3328ae24-6022-4042-9edc-79de6b40ed1e">


## 3. Error check
From the two summray files generated, you can check for errors that may have occurred during processing.

## 4. Additional Note

If you simply want to get antibody-antigen pairs, run the `write_all_pairs.py` will be good enough. Nothing else needed.

Example data you will get are in the data foler for your reference.
