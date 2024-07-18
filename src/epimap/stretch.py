"""Stretch epitopes to all the positions they cover
"""
import numpy as np
import pandas as pd

def stretch(epitopes, length_col="length", start_col="start", seq_col="seq"):
    stretched = []
    epitopes = epitopes.copy()
    epitopes["position"] = epitopes.start
    for i in range(epitopes.length.max()):
        updated_pos = epitopes.copy()
        mask = updated_pos.length >= i
        updated_pos = updated_pos[mask]
        updated_pos.position = updated_pos.position + i
        updated_pos['residue'] = updated_pos[seq_col].apply(lambda x: x[i])#
        stretched.append(updated_pos)
    stretched = pd.concat(stretched)
    stretched = stretched.sort_values([start_col, "position"])
    return stretched

def make_empty_grid(stretched, index, row_col, seq_length, default_value=0):
    rows = stretched[row_col].unique()
    cols = range(index, seq_length+index)
    grid = np.zeros((len(rows), len(cols)))
    grid = grid + default_value
    grid = pd.DataFrame(grid)
    grid.index = rows
    grid.columns = cols
    return grid

def make_grid(stretched, grid_values, index, row_col, seq_length):
    grid = make_empty_grid(stretched, index, row_col, seq_length)
    for i,count in grid_values.items():
        grid.loc[i[0], i[1]] = count
    return grid

