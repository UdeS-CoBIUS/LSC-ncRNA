# LSC-ncRNA
 Large Scale Classification of non-coding RNA



# Experimentation Pipeline
## 1) get the dataset:
Download the ncRNA from somewhere, form example Rfam (https://rfam.xfam.org/)

## 2) Train Test split:
For each family (fasta file), split the sequnces into train and test dataset.
We use the code in ***constructTrainTestFiles***.
- dir input for fasta file, each file is a family.
- dir output.  the output will be in "dir output"/Train , and "dir output"/Test
- min_nb_seqs_allowed: minumum sequences in family, (let say 10).
- percentage_nb_seqs_train: let say 70%.

The output will be, 1) a folder Train with *percentage_nb_seqs_train* sequnces, 2) and the same goes for Test folder with *100%-percentage_nb_seqs_train*. 3) a csv file that conatins some useful information.

**note**, be aware of problem of file encoding, and endline. sometimse it causes problems. In my case, I have to change to encoding to UTF8 to avoid extra-space when reading a line (see my answer: https://stackoverflow.com/a/73952980/3429103)

## 3) Motifs computation and selection:
which generate a 2d vector (matrix) at the end, of motifs and their number of occurences.

## 4) Train, and test and generate the score values:



### supplementary sub-folder: 
contains supplementary information and data for our paper.