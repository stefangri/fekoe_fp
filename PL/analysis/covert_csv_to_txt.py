import numpy as np
import os

listOfFiles  = os.listdir('csv_files')

for file in listOfFiles:
    with open('csv_files/' + file, 'r') as r:
        with open(file.split('.')[0] + '.txt', 'w') as w:
            for line in r:
                w.write(line.replace(',', '.').replace('E', 'e').replace('[', '#['))
