
import itertools

def generate_kMer_MinMax(min=3, max=6):
    total = 0
    list_kMer = []
    list_i = None
    for i in range(min, max + 1):
        list_i = [''.join(i) for i in itertools.product('ACGU', repeat=i)]
        list_kMer += list_i
        print(i, " : ", len(list_i))
        total += len(list_i)

    print(" Total = ", total)
    return list_kMer

def get_k_mer(str, k):
    lis_k_mer = []

    for i in range(0, len(str) - k):
        print(str[i:i + k])
        lis_k_mer.append(str[i:i + k])

    return lis_k_mer

def get_all_k_mer(str, mink, maxk):
    lis_k_mer = []

    if (maxk > len(str)):
        maxk = len(str)

    for k in range(mink, maxk + 1):
        for i in range(0, len(str) - k):
            print(str[i:i + k])
            lis_k_mer.append(str[i:i + k])

    return lis_k_mer
