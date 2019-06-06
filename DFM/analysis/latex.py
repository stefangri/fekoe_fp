from pandas import Series, DataFrame
import pandas as pd
import collections
import numpy
import uncertainties
import pint
from uncertainties import ufloat
from pint import UnitRegistry
import os.path
ureg = UnitRegistry()
Q_ = ureg.Quantity



def return_int(num):
    num_str = str(num)
    
    num_str = num_str.split('.')[1]
    num_str = num_str[0:1]
    return int(num_str)

def abs_path(filename):
    return os.path.join(os.path.dirname(__file__), filename)


def sign_digits(num):
    chars = f'{num}'
    chars = chars.split('+/-')[-1]
    chars = list(chars)
    while '.' in chars: chars.remove('.')
    if 'e' in chars:
        return chars[0] + chars[1]
    else:
        return chars[-2] + chars[-1]

def best_value(num):
    value = f'{num}'.split('+/-')[0]
    return value




class Latexdocument(object):
    def __init__(self, filename):
        self.name = filename
        self.data = DataFrame(columns=(['tex', 'var']))
    

    def add_result(self, name, value):
            if (type(value.magnitude) == uncertainties.core.Variable or type(value.magnitude) == uncertainties.core.AffineScalarFunc):
                latex_value = f'{value.magnitude:+.2uS}'
                latex_unit = (f'{value:Lx}'.split('}{'))[1].split('}')[0]
                value = value.magnitude
                df = DataFrame({'var': pd.Series(value, index = [name]),
                'tex': '\SI{' + latex_value + '}{' + latex_unit + '}'})

            else:
                latex_unit = (f'{value:Lx}'.split('}{'))[1].split('}')[0]
                latex_value = str(value.magnitude)
                df = DataFrame({'var': pd.Series(value, index = [name] ),
                'tex': '\SI{' + latex_value + '}{' + latex_unit + '}'})

            self.data = self.data.append(df, sort = True)
            with open(abs_path('results/result_' + name.replace('\\', '') + '.tex'), 'w') as f:
                #f.write('\\begin{equation} \n')
                f.write(self.data['tex'][name] + '\n')
                #f.write('\label{eq: result_' +  name +  '}\n')
                #f.write(r'\end{equation}')


    
