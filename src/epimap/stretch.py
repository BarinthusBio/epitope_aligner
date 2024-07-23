"""Stretch epitopes to all the positions they cover
"""
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import linkage, leaves_list
import itertools

def stretch(epitopes, length_col="length", start_col="start", seq_col="seq"):
    stretched = []
    epitopes = epitopes.copy()
    epitopes["position"] = epitopes.start
    for i in range(epitopes[length_col].max()):
        updated_pos = epitopes.copy()
        mask = updated_pos.length > i
        updated_pos = updated_pos[mask]
        updated_pos.position = updated_pos.position + i
        updated_pos['residue'] = updated_pos[seq_col].apply(lambda x: x[i])#
        stretched.append(updated_pos)
    stretched = pd.concat(stretched)
    stretched = stretched.sort_values([start_col, "position"])
    return stretched

def add_empty_positions(series, index, seq_length, empty_value, position_name="position"):
    assert position_name in series.index.names, f"Expected {position_name} in series.index.names"
    series = series.copy()
    full_positions = pd.Series(range(index, seq_length+index))
    names = [name for name in series.index.names if name != position_name]
    index_levels = [series.index.unique(name) for name in names]
    for levels in itertools.product(*index_levels):
        series.sort_index(inplace=True)
        try:
            missing_positions = full_positions[~full_positions.isin(series.loc[levels].index)]
        except KeyError:
            missing_positions = full_positions.copy()
        for i in missing_positions:
            series.loc[levels + (i,)] = empty_value
    series = series.sort_index()
    return series


def _non_position_index(index, position_col):
    names = list(index.names)
    if len(names) != 2:
        raise AssertionError(f"Expected two index names, got {names}")
    names.remove(position_col)
    non_position_name = names[0]
    return non_position_name


def order_grid(grid):
    linkage_data = linkage(grid, method="ward", metric="euclidean")
    order = leaves_list(linkage_data)
    ordered_grid = grid.iloc[order]
    return ordered_grid


def make_grid(grid_values, index, seq_length, empty_value, position_col="position", row_col=None):
    grid_values = add_empty_positions(
        series=grid_values,
        seq_length=seq_length,
        index=index,
        empty_value=empty_value
    )
    if not row_col:
        row_col = _non_position_index(grid_values.index, position_col=position_col)
        grid_values = grid_values.reset_index()
    grid = grid_values.pivot(index=row_col, columns=position_col)
    grid = order_grid(grid)
    return grid


# def make_empty_grid(stretched, index, row_col, seq_length, default_value):
#     rows = stretched[row_col].unique()
#     cols = range(index, seq_length+index)
#     grid = np.zeros((len(rows), len(cols)))
#     grid = grid + default_value
#     grid = pd.DataFrame(grid)
#     grid.index = rows
#     grid.columns = cols
#     return grid

# def make_grid(stretched, grid_values, index, row_col, seq_length, default_value=0):
#     grid = make_empty_grid(stretched, index, row_col, seq_length, default_value)
#     for i,count in grid_values.items():
#         grid.loc[i[0], i[1]] = count
#     return grid

