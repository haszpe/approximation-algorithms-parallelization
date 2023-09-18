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

import pandas as pd
from datetime import datetime
from multiprocessing import Pool

zawodnicy = (['rat_one', 1.4],
             ['rat_two', 1.2],
             ['rat_three', 2.4],
             ['rat_four', 0.9])

def race(x):
    name = x[0]
    tempo = x[1]
    now = datetime.now()
    print(f"{name} started race at {now.time()}  ")
    progres = 0
    while progres < 100:
        progres += tempo
        result = 0
        for i in range(0, 1000000):
            result += i ** 2
    now = datetime.now()
    print(f"{name} finished race at {now.time()}")
    return

if __name__ == '__main__':

    start = datetime.now()
    test = []
    for x in range(0, 11):
        result = 0
        start = datetime.now()
        for i in range(0, 400000):
            result += i ** 2
        czas = datetime.now() - start
        td = pd.Timedelta(czas, unit='ns')
        test.append(td)
    avg = pd.to_timedelta(pd.Series(test)).mean()
    print(f"Średni czas wykonania jednego kroku: {avg} [CPU-bound task]")

    start = datetime.now()
    pool = Pool(processes=4)
    pool.map(race, zawodnicy)
    pool.close()
    pool.join()
    print(f"Zrównoleglone: {(datetime.now() - start).seconds} s")

    start = datetime.now()
    for x in zawodnicy:
        race(x)
    print(f" Sekwencyjnie: {(datetime.now() - start).seconds} s")
