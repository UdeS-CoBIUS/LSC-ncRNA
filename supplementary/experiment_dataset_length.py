# experiment_dataset_length:
#  generate a graph where on the X axis put family size (number of seqs),
#  and on Y axis the average len of seqs.
#  this for 3 dataset:
#  1) Rfam dataset,
#  2) Noise dataset for 0 noise,
#  3) Clan dataset (and for Rfam from clans)
#
# 	  for that we have compute for each family 1) number of seqs, and the average len seqs.


# count the number of sequences in each family
# with simply counting the sign of sequence id ">"
# we can use : `grep -c ">" file.fasta` but in shell, with `subprocess`
def count_number_of_seqs(file_family):
    nb_seqs = 0
    with open(file_family) as f:
        for line in f:
            if line.startswith('>'):
                nb_seqs += 1

    return nb_seqs

# this form : http://biopython.org/DIST/docs/tutorial/Tutorial.html#sec:low-level-fasta-fastq
def count_and_average_len_seqs(file_family):
    
    from Bio.SeqIO.FastaIO import SimpleFastaParser
    count = 0
    total_len = 0
    with open(file_family) as in_handle:
        for title, seq in SimpleFastaParser(in_handle):
            count += 1
            total_len += len(seq)
            #print(" seq num =", count, " | len =",len(seq))
            
    return count, total_len


def main():
    file_test = r"C:\Users\ibra\Desktop\Infernal\Clans ncRNA\CL00001\RF00005.fasta"
    nb_seqs, total_seqs_len = count_and_average_len_seqs(file_test)

    print(" nb seqs = ", nb_seqs)
    print(" total seqs len = ", total_seqs_len)
    print(" average len = ", total_seqs_len/nb_seqs)



if __name__ == '__main__':
    main()