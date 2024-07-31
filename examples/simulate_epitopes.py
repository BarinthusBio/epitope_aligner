import numpy as np
import pandas as pd
from epimap import map,stretch,utils

sequence = utils.random_seq(10)
epitopes = utils.random_epitopes(sequence, 3, (3,6), index=0, includeend=False)
stretch.stretch(epitopes)

epitopes = utils.random_epitopes(sequence, 3, (3,6), index=1, includeend=False)
stretch.stretch(epitopes)

epitopes = utils.random_epitopes(sequence, 3, (3,6), index=0, includeend=True)
stretch.stretch(epitopes)

epitopes = utils.random_epitopes(sequence, 3, (3,6), index=1, includeend=True)
stretch.stretch(epitopes)


sequence
epitopes
