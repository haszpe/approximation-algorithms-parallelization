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

def levenshtein_distance(matrix, num_current, num_next_right, num_next_down):
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

    #edit distance
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
                if num == 1:
                    granica_up = list(range(0, len(x)+1))
                    granica_up_save = granica_up
                else:
                    granica_up = list(range(granica_up_save[-1], granica_up_save[-1]+len(x)+1))
                    granica_up_save = granica_up
            else:
                granica_up = None
            # Domyślne uzupełnienie granic lewych dla pierwszej kolumny submatryc
            if num in pierwsza_kolumna:
                if num == 1:
                    granica_left = list(range(0, len(y)+1))
                    granica_left_save= granica_left
                else:
                    granica_left = list(range(granica_left_save[-1], granica_left_save[-1]+len(y)+1))
                    granica_left_save = granica_left
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


def levenshtein_distance_worker(matrix, num_current, num_next_right, num_next_down):
    return levenshtein_distance(matrix, num_current, num_next_right, num_next_down)

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
    processes = 2

    #Pomiary zasobów
    start = datetime.now()
    process = psutil.Process(os.getpid())

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
            ld = levenshtein_distance(matrix,
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
            ld_results = pool.starmap(levenshtein_distance_worker, tasks)
            pool.close()
            pool.join()
            results.extend(ld_results)

    print(f"Memory usage percent before merging and filtering: {process.memory_percent()}")

    #Połaczenie submacierzy w jedną dużą macierz
    indicators = np.arange(1, processes * processes + 1).reshape((processes, processes))
    merged_rows = []
    for row_idx in range(processes):
        row_matrices = [matrix for number, matrix in results if number in indicators[row_idx]]
        merged_row = np.hstack(row_matrices)
        merged_rows.append(merged_row)
    result_matrix = np.vstack(merged_rows)

    #Odfiltrowanie powtórzonych wartości granicznych
    unique_rows, row_counts = np.unique(result_matrix, axis=0, return_counts=True)
    duplicated_row_indices = np.where(row_counts > 1)[0]
    unique_cols, cols_counts = np.unique(result_matrix, axis=1, return_counts=True)
    duplicated_cols_indices = np.where(cols_counts > 1)[0]
    filtered = np.delete(result_matrix, duplicated_row_indices, axis=0)
    filtered = np.delete(filtered, duplicated_cols_indices, axis=1)

    print(f"Levenshtein edit distance with dict parallelized with {processes} processes took: "
          f"{(datetime.now() - start).seconds} s")

    print("Levenshtein matrix:")
    print(filtered)

    print(f"Memory usage percent after merging and filtering: {process.memory_info().rss} bytes")
