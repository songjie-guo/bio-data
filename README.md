# Biology Data Processing (from SAbDAb & Protein Data Bank)

Generate the antibody's corresponding antigen .fasta and .pdb file (can be used for antibody)

## 1. Framework
1. Select antibody from SAbDAb (e.g. structure search), get "xxx.tsv" file.
2. From "xxx.tsv" file, get all pdb codes (e.g.'1a2y')
3. Download and store all original fasta files and pdb files (**NOT** included in this repo)
4. Use the code in this repo
5. Check the error log. **Manually fix errors.**

## 2. Implementation 
### 2.1 Preparation
You first need to download and store all original fasta files and pdb files. You may use your own methods X)

### 2.2 File Structre
Your folders look like this:

```
code/
  antigen.py
  process.py
  main.py

bio-data/ (you store original fasta & pdb files)
  ori_fasta/
    1a2y.fasta
    ...
  ori_pdb/
    1a2y.pdb 
    ...

bio-data-test/ (change the path if needed and change corresponding codes)
  antigen_fasta/
    1a2y_1.fasta 
    ...
  antigen_pdb/
    1a2y_1.pdb 
    ...
  log/
    error.log
```

Note that all data are stored in a seperate folder, for better management

## 3. Ckeck Errors
The `error.log` shows all possible errors that may occur during the processing, so that you may take a further look and manually adjust it.

Explanation:
|error|meaning|
| ------------- | ------------- |
|'No valid ori_fasta.'|The fasta file you download is not valid/empty.|
|'No valid ori_pdb.'|The pdb file you download is not valid/empty.|
|'No fasta generated.'|package `anarci` cannot find an antigen chain.|
|'No pdb generated.'|The antigen chain name found from fasta file cannot correspond to any chain_id in a pdb file. (Their naming may be totally different)|

## 4. Notes
`antigen.py` defines a class called "**Antigen**". Actually, by adding several lines, you can add antibody property into it, and may change the class to a genenral "**Protein**".

