"""
Program wykonany na potrzeby pracy licencjackiej:
"Porównanie sekwencyjnych i równoległych algorytmów przyrównujących sekwencje genomowe"
Autorka: Hanna Pęciak
----------------------------------------------------
Uniwersytet Przyrodniczy we Wrocławiu
Wydział Biologii i Hodowlii Zwierząt
Bioinformatyka
Wrocław 2023
"""

import numpy as np
import psutil
import os
from datetime import datetime

def seq_smith_waterman(x, y, match=2, mismatch=-1, gap=-2):
    m = len(x)
    n = len(y)

    # macierz zer
    d = np.zeros((m+1, n+1), dtype='int32')

    # uzupełnienie macierzy zgodnie z algorytmem
    maximum = 0
    max_i, max_j = 0, 0

    for i in range(1, m+1):
        for j in range(1, n+1):
            if x[i-1] == y[j-1]:
                c = match
            else:
                c = mismatch

            d[i][j] = max(0,
                          d[i-1][j-1] + c,
                          d[i-1][j] + gap,
                          d[i][j-1] + gap)

            if d[i][j] > maximum:
                maximum = d[i][j]
                max_i = i
                max_j = j

    # poszukiwanie ścieżki
    dopasowana_1 = ""
    dopasowana_2 = ""
    i = max_i
    j = max_j

    while d[i][j] > 0:
        if x[i-1] == y[j-1]:
            dopasowana_1 = x[i-1] + dopasowana_1
            dopasowana_2 = y[j-1] + dopasowana_2
            i -= 1
            j -= 1

        elif d[i][j] == d[i-1][j-1] + mismatch:
            dopasowana_1 = x[i-1] + dopasowana_1
            dopasowana_2 = y[j-1] + dopasowana_2
            i -= 1
            j -= 1

        elif d[i][j] == d[i-1][j] + gap:
            dopasowana_1 = x[i-1] + dopasowana_1
            dopasowana_2 = "-" + dopasowana_2
            i -= 1

        else:
            dopasowana_1 = "-" + dopasowana_1
            dopasowana_2 = y[j-1] + dopasowana_2
            j -= 1

    return dopasowana_1, dopasowana_2


def read_fasta(filename):
    seq = ''
    file_path = filename
    with open(file_path, "r") as file:
        for line in file:
            line = line.strip()
            if not line.startswith(">"):
                seq += line
    return seq


if __name__ == "__main__":
    input_one = input("Enter the name of first file in FASTA format (e.g., AT2G39770.1): ")
    input_two = input("Enter the name of second file in FASTA format (e.g., AT2G39770.2): ")

    A = read_fasta(input_one)
    B = read_fasta(input_two)

    # Pomiary zasobów
    start = datetime.now()
    process = psutil.Process(os.getpid())

    dopasowana_1, dopasowana_2 = seq_smith_waterman(A, B)
    print(f"Sequenced Smith-Waterman algorithm took: {(datetime.now() - start).seconds} s")

    print("Smith-Waterman:")
    print("Aligned Sequence X:", dopasowana_1)
    print("Aligned Sequence Y:", dopasowana_2)

    print(f"Memory usage percent: {process.memory_info().rss} bytes")
