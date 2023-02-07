# LSC-ncRNA
Large Scale Classification of non-coding RNA

# General use


# Experimentation Pipeline

<!-- Requirements -->
<h3 id="requirements"> :hammer_and_pick: Requirements</h3>

*   __`C++14`__
*   __`python3 (at least python 3.6)`__
*   __`Biopython`__
*   __`datatable`__
*   __`sklearn`__
*   __`Pandas`__
*   __`Numpy`__

## A) Dataset preparation

Download the ncRNA from somewhere, form example Rfam (https://rfam.xfam.org/). In our case we used the following 3 datasets

### Dataset 1: Rfam 14.1 seed sequences with 3016 ncRNAs families
Download the from the ftp server on https://ftp.ebi.ac.uk/pub/databases/Rfam/14.1/Rfam.seed.gz

After that, we have to exatrct each family in a separate file. Where the file name is Rfam family name as Rfxxxxx (xxxxx: is a number). Also we change the sequences identifier inside the fasta file to the family name >Rfamxxxxx , so we can get the sequence class in our future processing.

The first dataset is split to 30% Test and 70% Train.
The part in 70%, is used to do cross validation exeperemetation, where we generate a sample with 600 families, and a sample with 350 families. each of them is also split into 30% test and 70% train.

### Dataset 2: 88 Rfam with noise
This dataset is available at https://github.com/bioinformatics-sannio/ncrna-deep. This benchmark dataset is used in the method deepncRNA in the paper Deep learning predicts short non-coding RNA functions from only raw sequence data ([paper link](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1008415)).

The dataset includes 306,016 sequences distributed among 88 RFAM families of small ncRNA such that each family contains more than 400 sequences.  Each dataset class was split into three random subsets for train (84\%), validation (8\%), and test (8\%) in the deepncRNA method. We kept the initial partition of the dataset to be able to compare our results with those of deepncRNA.

***Noise datset*** : For each sequence in the test dataset, we added x\% noise at its start and its end, i.e., x/2\% at each extremity, such that the percentage x\% of noise compared to the length of the original sequence varied between 0\% and 200\% with steps of 25\%. The noise was generated by randomly shuffling the sequence while preserving the nucleotide and dinucleotide frequency of the original sequence. (The noise dataset is generated by deepncRNA method).


### Dataset 3: Clans from Rfam 14.8
- ***Low sequence similarity dataset*** We construct the low sequence similarity dataset included $36$ classes defined as the clans of RFAM 14.8 containing at least $4$ ncRNA families (https://rfam.xfam.org/clan/). This resulted in $36$ clans containing $4$ to $11$ families per clan and a total of $199$ families for all $36$ clans. For each clan, we selected $12$ sequences from each family.

***The high sequence similarity dataset*** included $199$ classes defined as the $199$ ncRNA families contained in the $36$ clans from the low sequence similarity dataset.

*** List of used clans: ***: *["CL00051", "CL00003", "CL00069", "CL00106", "CL00038", "CL00054", "CL00002", "CL00014", "CL00102", "CL00001", "CL00117", "CL00021", "CL00057", "CL00111", "CL00118", "CL00066", "CL00112", "CL00116", "CL00010", "CL00005", "CL00027", "CL00053", "CL00063", "CL00100", "CL00119", "CL00004", "CL00032", "CL00034", "CL00035", "CL00040", "CL00096", "CL00093", "CL00015", "CL00121", "CL00036", "CL00045"]*

Each class of each dataset was split into two random subsets for train (70\%) and test (30\%).

We use the script in ***dataset/clans_dataset.py*** to download and prepare the dataset.

### Datset Train Test split and dataset statistics:

We use ***constructTrainTestFiles*** to split the dataset into train and test, and to get the dataset statistics.
- To split the dataset, for each family (fasta file), we split the sequences in it, into train and test file.
- dataset statistics: for a given folder, for each family (in separate file) get the follwing informations: nb seqs, min seq len, max seq len, average seq len , and save to csv file.

The program can be use like this:
constructTrainTestFiles -in \<string\> [-out \<string\> -nf \<integer\> -mins \<integer\> -maxs \<integer\> -pt \<integer\> -m \<string\>]

- *constructTrainTestFiles* : the program name

- -in : \<string\> a path directory for fasta files.
- -out : \<string\> a main path directory the results out.
- -nf : \<integer\> number of families.
- -mins : \<integer\>, min number of sequences (default 4).
- -maxs : \<integer\>, max number of sequences.
- -pt : \<integer\>, percentage number sequences Test, defult value is 30.

- -m : \<string\> the used mode, several mode are avaialable: 
    - -m i: information, get all informations as nb seqs, min seq len, max seq len, average seq len , and save to csv file.
    - -m s: Sample dataset, get n random families that have nb seqs between min and max, and save them to dir_output
    - -m sttmm : Split Tarin Test Min Max, for a given nb of families, and min max number of seqs, split to train and test
    - -m sttm : Split Tarin Test Min, consider only min number of seqs, and split all files in input folder to train and test
    - -m stt : Split Tarin Test, for all files in input folder split to train and test

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


:busts_in_silhouette: __Authors__
`ibrahim Chegrane & Nabil Bendjaa & Aida Ouangraoua`, CoBIUS LAB, Department of Computer Science, Faculty of Science, Université de Sherbrooke, Sherbrooke, Canada

> :bulb: If you are using our algorithm in your research, please cite our recent paper: __Upcoming__ 

> :e-mail: Contact: ibrahim[dot]chegrane[at]usherbrooke[dot]ca
