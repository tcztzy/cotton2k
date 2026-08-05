//! Bounded iteration helpers for the model's fixed-size ndarray grids.
//!
//! These helpers centralize row, cell, layer-span, and fruiting-site traversal
//! so domain code respects the active dimensions stored in legacy state. They
//! borrow arrays for the duration of a callback and never change grid shape.

use ndarray::{Array1, Array2, ArrayView1, ArrayViewMut1, Axis};

/// Visits each active row of a two-dimensional array.
pub fn for_each_row<T>(
    grid: &Array2<T>,
    nrows: usize,
    mut f: impl FnMut(usize, ArrayView1<'_, T>),
) {
    for (row_idx, row) in grid.axis_iter(Axis(0)).take(nrows).enumerate() {
        f(row_idx, row);
    }
}

/// Visits each active row of a two-dimensional array mutably.
pub fn for_each_row_mut<T>(
    grid: &mut Array2<T>,
    nrows: usize,
    mut f: impl FnMut(usize, ArrayViewMut1<'_, T>),
) {
    for (row_idx, row) in grid.axis_iter_mut(Axis(0)).take(nrows).enumerate() {
        f(row_idx, row);
    }
}

/// Visits every cell in a bounded row-major rectangle.
pub fn for_each_cell(nrows: usize, ncols: usize, mut f: impl FnMut(usize, usize)) {
    for row in 0..nrows {
        for col in 0..ncols {
            f(row, col);
        }
    }
}

/// Visits every column in an inclusive-start, exclusive-end row range.
pub fn for_each_cell_in_rows(
    row_start: usize,
    row_end: usize,
    ncols: usize,
    mut f: impl FnMut(usize, usize),
) {
    for row in row_start..row_end {
        for col in 0..ncols {
            f(row, col);
        }
    }
}

/// Visits each cell between per-layer left and right column bounds.
pub fn for_each_layer_col_span(
    nrows: usize,
    ncols: usize,
    left_cols: &[i32],
    right_cols: &[i32],
    mut f: impl FnMut(usize, usize),
) {
    for row in 0..nrows {
        let left = left_cols.get(row).copied().unwrap_or_default().max(0) as usize;
        let right = right_cols
            .get(row)
            .copied()
            .unwrap_or(-1)
            .min((ncols as i32) - 1);
        if right < left as i32 {
            continue;
        }
        for col in left..=right as usize {
            f(row, col);
        }
    }
}

/// Visits active fruiting sites and passes each site's active node count.
pub fn for_each_fruiting_site(
    num_veg_branches: i32,
    num_fruit_branches: &Array1<i32>,
    num_nodes: &Array2<i32>,
    mut f: impl FnMut(usize, usize, usize),
) {
    for k in 0..num_veg_branches.max(0) as usize {
        for l in 0..num_fruit_branches[k].max(0) as usize {
            f(k, l, num_nodes[[k, l]].max(0) as usize);
        }
    }
}
