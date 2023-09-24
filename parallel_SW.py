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
from multiprocessing import Manager, Pool
import psutil
import os
from datetime import datetime

def SW_matrix(matrix, num_current, num_next_right, num_next_down, match=2, mismatch=-1, gap=-2):
    current_dict = matrix[f'dict_{num_current}']
    x = current_dict['x']
    y = current_dict['y']
    boundry_up = current_dict['granica_up']
    boundry_left = current_dict['granica_left']

    m = len(x)
    n = len(y)
    d = np.zeros((n+1, m+1), dtype='int32')
    assert(boundry_up[0] == boundry_left[0])

    #Uzupełnienie wierszy zerowych wartościami granicznymi z poprzednich submacierzy
    d[0, :] = boundry_up
    d[:, 0] = boundry_left

    for y_item in range(1, n+1):
        for x_item in range(1, m+1):
            if x[x_item-1] == y[y_item-1]:
                c = match
            else:
                c = mismatch
            d[y_item, x_item] = np.max([0,
                                        d[y_item-1][x_item] + gap,
                                        d[y_item][x_item-1] + gap,
                                        d[y_item-1][x_item-1] + c])

    #Nadpisanie kolejnych słowników submacierzy (jeśli istnieją)
    if num_next_right <= len(matrix):
        next_right = matrix[f'dict_{num_next_right}']
        if next_right['granica_left'] is None:
            next_right['granica_left'] = list(d[:, m])
            matrix[f'dict_{num_next_right}'] = next_right

    if num_next_down <= len(matrix):
        next_down = matrix[f'dict_{num_next_down}']
        if next_down['granica_up'] is None:
            next_down['granica_up'] = list(d[n, :])
            matrix[f'dict_{num_next_down}'] = next_down

    return (num_current, d)


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
                granica_up = [0 for _ in range(len(x)+1)]
            else:
                granica_up = None
            # Domyślne uzupełnienie granic lewych dla pierwszej kolumny submatryc
            if num in pierwsza_kolumna:
                granica_left = [0 for _ in range(len(y)+1)]
            else:
                granica_left = None

            dicts[dict_name] = {'x': x,
                                'y': y,
                                'granica_up': granica_up,
                                'granica_left': granica_left}
            num += 1
    return dicts


def create_acts(x):
    matrix = [[i + j + 1 for j in range(x)] for i in range(0, x * x, x)]
    num_rows = len(matrix)
    num_cols = len(matrix[0])
    result = []
    for i in range(num_rows + num_cols - 1):
        diagonal = []
        if i < num_rows:
            start_row, start_col = i, 0
        else:
            start_row, start_col = num_rows - 1, i - num_rows + 1
        while start_row >= 0 and start_col < num_cols:
            diagonal.append(matrix[start_row][start_col])
            start_row -= 1
            start_col += 1
        result.append(diagonal)
    return result


def SW_matrix_worker(matrix, num_current, num_next_right, num_next_down):
    return SW_matrix(matrix, num_current, num_next_right, num_next_down)


def smith_waterman_traceback(full_matrix, x, y, h, v, mismatch=-1, gap=-2):
    dopasowana_1 = ""
    dopasowana_2 = ""
    i = v
    j = h

    d = full_matrix.T
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

    processes = 2

    #Stworzenie wspólnego dla procesów źródła danych
    manager = Manager()
    matrix = manager.dict(dicts_creation(A, B, processes))

    akty = create_acts(processes)
    num_dicts = len(matrix)

    results = []

    for akt in akty:
        #Pojedyncze akty nie są zrównoleglane
        if len(akt) == 1:
            num = akt[0]
            num_current = num
            num_next_right = num + 1
            num_next_down = num + processes
            ld = SW_matrix(matrix,
                           num_current,
                           num_next_right,
                           num_next_down)
            results.append(ld)

        #Stworzenie rownoległych procesów dla danego aktu
        else:
            tasks = []
            for i in akt:
                num_current = i
                num_next_right = i + 1
                num_next_down = i + processes

                task_args = (matrix,
                             num_current,
                             num_next_right,
                             num_next_down)
                tasks.append(task_args)

            pool = Pool(processes)
            ld_results = pool.starmap(SW_matrix_worker, tasks)
            pool.close()
            pool.join()
            results.extend(ld_results)

    #Połaczenie submacierzy w jedną dużą macierz
    indicators = np.arange(1, processes * processes + 1).reshape((processes, processes))
    merged_rows = []
    for row_idx in range(processes):
        row_matrices = [matrix for number, matrix in results if number in indicators[row_idx]]
        merged_row = np.hstack(row_matrices)
        merged_rows.append(merged_row)
    result_matrix = np.vstack(merged_rows)

    #Odfiltrowanie powtórzonych wartości granicznych
    filtered = np.delete(result_matrix, len(B)//2, axis=0)
    filtered = np.delete(filtered, len(A)//2, axis=1)

    #Obliczenie maksymalnej wartości
    ind = np.unravel_index(np.argmax(filtered), filtered.shape)
    max_i = ind[0]
    max_j = ind[1]
    aligned_x, aligned_y = smith_waterman_traceback(filtered, A, B, max_i, max_j)

    print(f"Parallelized Smith-Waterman algorithm took: {(datetime.now() - start).seconds} s")

    print("Parallelized Smith-Waterman:")
    print("Aligned Sequence X:", aligned_x)
    print("Aligned Sequence Y:", aligned_y)

    print(f"Memory usage percent: {process.memory_info().rss} bytes")
