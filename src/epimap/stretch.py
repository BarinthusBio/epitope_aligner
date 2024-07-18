"""Stretch epitopes to all the positions they cover
"""
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import linkage, leaves_list

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

def add_empty_positions(series, seq_length, index, empty_value, position_col="position"):
    def get_func(series):
        if series.index.nlevels ==1:
        # isinstance(series.index, pd.core.indexes.base.Index):
            return add_empty_positions_index
        elif series.index.nlevels > 1:
        # isinstance(series.index, pd.core.indexes.multi.MultiIndex):
            return add_empty_positions_multiindex
        else:
            raise ValueError(type(series.index))
    def add_empty_positions_index(series, seq_length, index, empty_value, **kwargs):
        series = series.copy()
        for i in range(index, seq_length+index):
            if not i in series.index:
                series.loc[i] = empty_value
        series = series.sort_index()
        return series
    def add_empty_positions_multiindex(series, seq_length, index, empty_value, position_col="position"):
        series = series.copy()
        names = list(series.index.names)
        names.remove("position")
        for name in names:
            for level in series.index.get_level_values(name).unique():
                for i in range(index, seq_length+index):
                    if not i in series.loc[level].index:
                        series.loc[(level, i)] = empty_value
        series = series.sort_index()
        return series
    func = get_func(series)
    series = func(
        series=series,
        seq_length=seq_length,
        index=index,
        empty_value=empty_value,
        position_col=position_col
    )
    return series

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

def order_grid(grid):
    linkage_data = linkage(grid, method="ward", metric="euclidean")
    order = leaves_list(linkage_data)
    ordered_grid = grid.iloc[order]
    return ordered_grid
