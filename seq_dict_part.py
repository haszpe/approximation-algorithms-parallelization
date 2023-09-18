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

#Funkcja aktualizująca słownik na podstawie odległości levenshteina
def levenshtein_distance(x, y, boundry_up, boundry_left, next_right=None, next_down=None, last = False):
    m = len(x)
    n = len(y)
    d = np.zeros((n+1, m+1), dtype='int32')
    assert(boundry_up[0] == boundry_left[0])

    #Uzupełnienie wierszy zerowych wartościami granicznymi z poprzednich submacierzy
    d[0, :] = boundry_up
    d[:, 0] = boundry_left

    #Edit distance
    for y_item in range(1, n+1):
        for x_item in range(1, m+1):
            if x[x_item-1] == y[y_item-1]:
                c = 0
            else:
                c = 1
            d[y_item, x_item] = np.min([d[y_item-1][x_item] + 1,
                                        d[y_item][x_item-1] + 1,
                                        d[y_item-1][x_item-1] + c])

    #Nadpisanie kolejnych słowników submacierzy (jeśli istnieją)
    if next_down:
        if next_down['granica_up'] == None:
            next_down['granica_up'] = list(d[n, :])
    if next_right:
        if next_right['granica_left'] == None:
            next_right['granica_left'] = list(d[:, m])

    #Zwraca ostatnią wartość w macierzy jeśli to ostatnia submacierz
    if last:
        return d[n, m]
    else:
        return


#Tworzenie słownika na podstawie tekstów i liczby dostępnych procesów
def dicts_creation(A, B, processes):
    x_parts = len(A)//processes
    y_parts = len(B)//processes

    pierwszy_wiersz = list(range(1, processes + 1))
    pierwsza_kolumna = list(range(1, (processes**2) + 1, processes))

    xs = []
    ys = []

    #Podział stringów na równe części
    for i in range(0, processes):
        if i == processes-1:
            xs.append(A[(i * x_parts):])
            ys.append(B[(i * y_parts):])
        else:
            xs.append(A[(i*x_parts):((i + 1) * x_parts)])
            ys.append(B[(i * y_parts):((i + 1) * y_parts)])

    #Tworzenie słowników
    dicts = {}
    num = 1
    for y in ys:            #po kolumnach
        for x in xs:        #po wierszach
            dict_name = f"dict_{num}"

            #Domyślne uzupełnienie granic górnych dla pierwszego wiersza submatryc
            if num in pierwszy_wiersz:
                if num == 1:
                    granica_up = list(range(0, len(x)+1))
                else:
                    granica_up = list(range(len(granica_up)-1, len(granica_up)+len(x)))
            else:
                granica_up = None

            # Domyślne uzupełnienie granic lewych dla pierwszej kolumny submatryc
            if num in pierwsza_kolumna:
                if num == 1:
                    granica_left = list(range(0, len(y)+1))
                    granica_left_save= granica_left
                else:
                    granica_left = list(range(len(granica_left_save)-1, len(granica_left_save)+len(y)))
                    granica_left_save = granica_left
            else:
                granica_left = None

            dicts[dict_name] = {'x': x,
                                'y': y,
                                'granica_up': granica_up,
                                'granica_left': granica_left}
            num += 1
    return dicts

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
    processes = 2

    #Pomiary zasobów
    process = psutil.Process(os.getpid())
    start = datetime.now()

    matrix = dicts_creation(A, B, processes)
    dict_keys = list(matrix.keys())

    for i, key in enumerate(dict_keys):

        current_dict = matrix[key]

        next_right_key = dict_keys[i + 1] if i + 1 < len(dict_keys) else None
        next_right_dict = matrix[next_right_key] if next_right_key else None

        next_down_key = dict_keys[i + processes] if i + processes < processes**2 else None
        next_down_dict = matrix[next_down_key] if next_down_key else None

        #Dla ostatniej matrycy - zmiana wartości last na True
        if i < processes**2-1:
            levenshtein_distance(current_dict['x'], current_dict['y'], current_dict['granica_up'],
                                 current_dict['granica_left'], next_right_dict, next_down_dict, False)
        elif i == processes**2-1:
            ld = levenshtein_distance(current_dict['x'], current_dict['y'], current_dict['granica_up'],
                                 current_dict['granica_left'], next_right_dict, next_down_dict, True)
            print(f"Levenshtein edit distance with dict between given sequences equals: {ld}")
            print(f"Levenshtein edit distance with dict took: {(datetime.now() - start).seconds} s")

    process = psutil.Process(os.getpid())
    print(f"Memory usage percent: {process.memory_info().rss} bytes")
