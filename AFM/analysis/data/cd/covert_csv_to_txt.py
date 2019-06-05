import numpy as np
import os

#listOfFiles  = os.listdir('csv_files')


with open('cd_profile_korrektur.txt', 'r') as r:
   with open('cd_profile_korrektur_new.txt', 'w') as w:
        for line in r:
            w.write(line.replace(',', '.').replace('E', 'e').replace('[', '#['))
