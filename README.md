# LSC-ncRNA
 Large Scale Classification of non-coding RNA



# Experimentation Pipeline
## A) dataset preparation
## A.1) get the dataset:
Download the ncRNA from somewhere, form example Rfam (https://rfam.xfam.org/)

## A.2) datset Train Test split:
For each family (fasta file), split the sequnces into train and test dataset.
We use the code in ***constructTrainTestFiles***.
- dir input for fasta file, each file is a family.
- dir output.  the output will be in "dir output"/Train , and "dir output"/Test
- min_nb_seqs_allowed: minumum sequences in family, (let say 10).
- percentage_nb_seqs_train: let say 70%.

The output will be, 1) a folder Train with *percentage_nb_seqs_train* sequnces, 2) and the same goes for Test folder with *100%-percentage_nb_seqs_train*. 3) a csv file that conatins some useful information.

**note**, be aware of problem of file encoding, and endline. sometimse it causes problems. In my case, I have to change to encoding to UTF8 to avoid extra-space when reading a line (see my answer: https://stackoverflow.com/a/73952980/3429103)

## B) Motifs computation and selection:
This is the first step in our method whoch is the computation and the selection of sequence motifs that allow defining a vectorial representation of ncRNAs sequences.
In this part as result we generate a 2d vector (matrix), that contain the number of occurences for each motif that exist in the ncRNA sequence. The matrix is saved as a single csv file.

The program can be used like this:

./MatrixCmsNStrNbOccs -in \<string\> [-nf \<integer\> -mins \<integer\> -maxs \<integer\> -minl \<integer\> -maxl \<integer\> -d \<integer\> -b \<integer\> -a \<integer\> -g \<integer\> -tn \<string\>]

- *MatrixCmsNStrNbOccs* : the program name
- *-in* : \<string\> a path directory for fasta files.
- *-tn* : \<string\> test name, give a specific name to you your excrement, (default "test")
- *-nf* : \<integer\> number of families, (**default 10**)
- *-mins* : \<integer\>, min number of sequences, (**default 4**)
- *-maxs* : \<integer\>, max number of sequences, (**default 1000**)
- *-minl* : \<integer\>, min length of motif, (**default 2**)
- *-maxl* : \<integer\>, max length of motif, (**default 10**)
- *-d* : \<integer\> (0: false, 1 or other: true), is delete sub-motifs, (**default 0**)
- *-b* : \<integer\> beta (between [0 and 100]), (**default 40**)
- *-a* : \<integer\>, alpha  (-1 no alpha, or: 0 equal number of occurrences, or 1,2,3,..., , (**default -1**)
- *-g* : \<integer\> ( >=1), gamma, number of occurrences allowed, (**default 1**)

All the parameters between [] are optional. The path to the directory of fast files is obligatory. The motifs length (minl, maxl) and Beta parameter and test name parameters are recommended to use.

Example:
```shell
nohup ./MatrixCmsNStrNbOccs -in "/data/ibra/Rfam_14.1_dataset/Rfam14.1_Sample_Train_Test/Rfam_600_Train_Test/Train" -minl 2 -maxl 8 -b 50 -g 1 -tn F_600 > out_F_600 &
```

The output csv file name is as follows: del_[No/Yes:-d]_nbF_[test_name:-tn]_min_[-minl]_max_[-maxl]_beta_[-b]_alpha_[-a]_nbOccrs_[-g].csv, the prvious test produce the following name: del_No_nbF_F_600_min_2_max_8_beta_50_alpha_-1_nbOccrs_1.csv


## C) Training and Testing experimental:

The step of the selection of supervised learning classification algorithms that allow achieving the most accurate classification of ncRNA sequences.

## C.1) selection of supervised learning classification algorithms:
To choose the best classification algorithm we use the following python program ***Classification/modelstest.py*** as follows:
`python3 modelstest.py path_motifs.csv`. The code test the following algorithm ['ext','knn','rdf','gnb', 'dt', 'nlp', 'svc'] using 10-fold cross-validation with a split of 0.3 for testing.

## C.2) Experimental with the chosen models.

To lucnh the program we use the main python script as follow:
*python3 Main.py mod "path/2d_matrix.csv" "path_test_folder"*, where:
- mod is "EXT" (Extra-tree), "RDF" (Random forest), "NLP" (MultiLayer perseptron) and "VOT" (voting model). 
- path/2d_matrix.csv: the path of the csv matrix generated using the previous program.
- path_test_folder: test folder, that contain a set of fasta files.

Example:

```shell
python3 Main.py EXT "/data/ibra/del_No_nbF_F_600_min_6_max_7_beta_50_alpha_-1_nbOccrs_1.csv" "/data/ibra/Rfam_14.1_dataset/Rfam14.1_Sample_Train_Test/Rfam_600_Train_Test/Test" > res_EXT_Single_del_No_nbF_Clans_min_6_max_7_beta_50_alpha_-1_nbOccrs_1
```

The training and test part, generate at the end different scores:
- Procissing time, of 1) trainig part and 2) testing part.
- Scores of training `score Train` and predection `pred Train`.
- pred Test (accuracy), Precision, Recall, and fbeta_score.



## 5) Blast based classification (blastn(av-s)):
Blast is search tool for the best hits. We use it and at to it some post processing to produce classification. We use a script to do the folloing steps:
- a) FRom the Train folder: gather all sequences in all files in fasta file X.
- b) use balst to creat a database on this one single file X (that contains all Train seqs)
- c) Gather all sequnces in Test folder in one single file Y.
- d) use blast to search sequnces in Y against database for X.
- e) Process balst result best hit to generate the classification.

To lunch the script `blast_classification.py`, we need to give as argument the path for Train and Test foder, and the test (experment) name as shown in the follwing example:
```shell
nohup python blast_classification.py "/data/chei2402/ibra/test_infernal/Clans_family_train_test/Train/" "/data/chei2402/ibra/test_infernal/Clans_family_train_test/Test/" Clans > res_blast_class_Clans_time_acc &
```


## 6) Infernal based classification:


### supplementary sub-folder: 
contains supplementary information and data for our paper.