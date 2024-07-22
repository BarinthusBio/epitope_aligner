import numpy as np
import pandas as pd
from epimap import map,stretch,utils

sequence = utils.random_seq(10)
epitopes = utils.random_epitopes(sequence, 9, (3,6), index=0, includeend=False)
sequence
epitopes

