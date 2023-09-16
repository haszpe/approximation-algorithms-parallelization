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

def seq_levenshtein_distance(x, y):
    m = len(x)
    n = len(y)

    # macierz zer
    d = np.zeros((m+1, n+1), dtype='int32')

    # uzupełnienie pierwszej kolumny
    for i in range(m+1):
        d[i][0] = i

    # uzupełnienie pierwszego wiersza
    for j in range(n+1):
        d[0][j] = j

    # uzupełnienie macierzy zgodnie z algorytmem
    for j in range(1, n+1):
        for i in range(1, m+1):
            if x[i-1] == y[j-1]:
                c = 0
            else:
                c = 1

            d[i][j] = min(
                d[i-1][j] + 1,      # delecja
                d[i][j-1] + 1,      # insercja
                d[i-1][j-1] + c     # tranzycja
            )

    lev_dist = d[m][n]
    return lev_dist


def read_fasta(filename):
    seq = ''
    file_path = filename
    with open(file_path, "r") as file:
        for line in file:
            line = line.strip()
            if not line.startswith(">"):
                seq += line
    return seq


if __name__ == '__main__':
    input_one = input("Enter the name of first file in FASTA format (e.g., AT2G39770.1): ")
    input_two = input("Enter the name of second file in FASTA format (e.g., AT2G39770.2): ")

    A = read_fasta(input_one)
    B = read_fasta(input_two)

    start = datetime.now()
    lev = seq_levenshtein_distance(A, B)
    print(f"Classic Levenshtein edit distance took: {(datetime.now() - start).seconds} s")

    print(f"Classic Levenshtein edit distance between given sequences equals: {lev}")
    process = psutil.Process(os.getpid())
    print(f"Memory usage percent: {process.memory_info().rss} bytes")
