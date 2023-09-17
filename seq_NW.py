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

def seq_needleman_wunsch(x, y, match=2, mismatch=-1, gap=-2):
    m = len(x)
    n = len(y)

    # macierz zer
    d = np.zeros((m+1, n+1), dtype='int32')

    # uzupełnienie pierwszej kolumny
    for i in range(1, m+1):
        d[i][0] = d[i-1][0] + gap

    # uzupełnienie pierwszego wiersza
    for j in range(1, n+1):
        d[0][j] = d[0][j-1] + gap

    # uzupełnienie macierzy zgodnie z algorytmem
    for i in range(1, m+1):
        for j in range(1, n+1):
            if x[i-1] == y[j-1]:
                diag = d[i-1][j-1] + match
            else:
                diag = d[i-1][j-1] + mismatch

            d[i][j] = max(diag,
                          d[i-1][j] + gap,
                          d[i][j-1] + gap)

    # optymalne dopasowanie
    dopasowana_1 = ""
    dopasowana_2 = ""
    i = m
    j = n

    while i > 0 and j > 0:
        if x[i-1] == y[j-1]:
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

    while i > 0:
        dopasowana_1 = x[i-1] + dopasowana_1
        dopasowana_2 = "-" + dopasowana_2
        i -= 1

    while j > 0:
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

    dopasowana_1, dopasowana_2 = seq_needleman_wunsch(A, B)
    print(f"Sequenced Needleman-Wunsch algorithm took: {(datetime.now() - start).seconds} s")

    print("Needleman-Wunsch:")
    print("Aligned Sequence X:", dopasowana_1)
    print("Aligned Sequence Y:", dopasowana_2)

    print(f"Memory usage percent: {process.memory_info().rss} bytes")

