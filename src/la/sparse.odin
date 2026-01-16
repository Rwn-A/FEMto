// SPDX-FileCopyrightText: 2025 Rowan Apps
// SPDX-License-Identifier: MIT
/*
 Data & functions for dealing with sparse matrices.
*/
package la

// TODO: simd or even threaded maybe spmv.

// Sparsity pattern represents which columns, in which rows, might be non-zero.
Sparsity :: struct {
	row_ptrs: []int,
	columns:  []int,
}


Sparse_Matrix :: struct {
	using sparsity: Sparsity,
	values:         []f64,
}

// Creates a new sparse matrix from given sparsity, does not copy sparsity.
sparse_mat_from_sparsity :: proc(sparsity: Sparsity, allocator := context.allocator) -> Sparse_Matrix {
	return {sparsity = sparsity, values = make([]f64, len(sparsity.columns), allocator)}
}

// Does not destroy sparsity pattern.
sparse_mat_destroy :: #force_inline proc(mat: Sparse_Matrix, allocator := context.allocator) {delete(mat.values, allocator)}

sparse_mat_rows :: #force_inline proc(sparsity: Sparsity) -> int {return len(sparsity.row_ptrs) - 1}

// Specific value index
sparse_mat_idx :: proc(sp: Sparsity, row, col: int) -> int {
	for i in sp.row_ptrs[row] ..< sp.row_ptrs[row + 1] {
		if sp.columns[i] == col {return i}
	}
	panic("Queried mat entry outside of sparsity")
}

//--Sparse Matrix Operations--


//y = a * A * x + b * y
sparse_mv :: proc(A: Sparse_Matrix, x, y: Vec, a := 1.0, b := 0.0) #no_bounds_check {
	assert(len(y) == sparse_mat_rows(A))
	assert(len(x) == len(y))

	if b == 0.0 {
		for &elem in y {elem = 0}
	} else if b != 1.0 {
		scal(y, b)
	}

	for row := 0; row < sparse_mat_rows(A); row += 1 {
		for idx := A.row_ptrs[row]; idx < A.row_ptrs[row + 1]; idx += 1 {
			col := A.columns[idx]
			y[row] += a * A.values[idx] * x[col]
		}
	}
	return
}
