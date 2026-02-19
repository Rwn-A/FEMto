// SPDX-FileCopyrightText: 2026 Rowan Apps
// SPDX-License-Identifier: MIT
package la

import "core:slice"
import "core:mem"
import "core:math"

Vector :: []f64

// Row-major
Dense_Matrix :: struct {
    values: []f64,
    rows, columns: int
}

Sparsity :: struct {
	row_ptrs: []int,
	columns:  []int,
}

Sparse_Matrix :: struct {
	using sparsity: Sparsity,
	values:         []f64,
}

Block_Vector :: struct {
    values: Vector,
    offsets: []int,
}

// Row-major
// Does not own the inner sparse matrices
Block_Sparse_Matrix :: struct {
    values: []^Sparse_Matrix,
    block_rows, block_columns: int,
}

// Row-major

Block_Dense_Info :: struct {
    offset: int,
    rows, columns: int
}

Block_Dense_Matrix :: struct {
    values: []f64,
    blocks: []Block_Dense_Info,
    block_rows, block_columns: int,
}

// Creation

dense_matrix_create :: proc(rows, columns: int, allocator := context.allocator) -> Dense_Matrix {
    return {
        rows = rows,
        columns = columns,
        values = make([]f64, rows * columns, allocator),
    }
}

sparse_matrix_from_sparsity :: proc(sp: Sparsity, allocator := context.allocator) -> Sparse_Matrix {
    return {
        sparsity = sp,
        values = make([]f64, len(sp.columns), allocator)
    }
}

block_vector_create :: proc(block_sizes: []int, allocator := context.allocator) -> Block_Vector {
    n := len(block_sizes)

    offsets := make([]int, n + 1, allocator)
    for i in 0..<n { offsets[i + 1] = offsets[i] + block_sizes[i] }

    return {
        make(Vector, offsets[n], allocator),
        offsets,
    }
}

block_vector_zero_clone :: proc(v: Block_Vector, allocator := context.allocator) -> Block_Vector {
    return {
        make(Vector, v.offsets[len(v.offsets) - 1], allocator),
        slice.clone(v.offsets),
    }
}

block_sparse_create :: proc(block_rows, block_columns: int, allocator := context.allocator) -> Block_Sparse_Matrix {
    return {
        block_rows = block_rows,
        block_columns = block_columns,
        values = make([]^Sparse_Matrix, block_rows * block_columns, allocator)
    }
}

block_dense_create :: proc(row_sizes, col_sizes: []int, allocator := context.allocator) -> Block_Dense_Matrix {
    block_rows, block_cols := len(row_sizes), len(col_sizes)

    blocks := make([]Block_Dense_Info, block_rows * block_cols, allocator)

    offset: int
    for br in 0..<block_rows {
        for bc in 0..<block_cols {
            block_index := br * block_cols + bc

            blocks[block_index].offset = offset
            blocks[block_index].rows = row_sizes[br]
            blocks[block_index].columns = col_sizes[bc]

            offset += row_sizes[br] * col_sizes[bc]
        }
    }

    return {
        values = make([]f64, offset, allocator),
        blocks = blocks,
        block_rows = block_rows,
        block_columns = block_cols,
    }
}

// misc

sparse_mat_rows :: #force_inline proc(sparsity: Sparsity) -> int {
    return len(sparsity.row_ptrs) - 1
}

sparse_mat_zero_row :: proc(mat: Sparse_Matrix, row: int) {
    slice.zero(mat.values[mat.row_ptrs[row]:mat.row_ptrs[row + 1]])
}

block_sparse_set_block :: proc(bs: Block_Sparse_Matrix, block_row, block_col: int, sm: ^Sparse_Matrix) {
    bs.values[block_row * bs.block_columns + block_col] = sm
}

// View

block_dense_view :: proc(m: Block_Dense_Matrix, row, col: int) -> Dense_Matrix {
    block := m.blocks[(row * m.block_columns) + col]
    return Dense_Matrix {
        rows = block.rows,
        columns = block.columns,
        values = m.values[block.offset:block.offset + (block.rows * block.columns)]
    }
}

block_sparse_view :: proc(m: Block_Sparse_Matrix, row, col: int) -> ^Sparse_Matrix {
    return m.values[(row * m.block_columns) + col]
}

block_vector_view :: proc(v: Block_Vector, block: int) -> Vector {
    return v.values[v.offsets[block]:v.offsets[block + 1]]
}

// indexing

dense_idx :: proc(m: Dense_Matrix, row, col: int) -> int {
    return row * m.columns + col
}

sparse_idx :: proc(sp: Sparsity, row, col: int) -> int {
    for i in sp.row_ptrs[row] ..< sp.row_ptrs[row + 1] {
        if sp.columns[i] == col {return i}
	}
	panic("Sparse matrix did not contain requested entry")
}

vector_idx :: proc(v: Vector, idx: int) -> int {
    return idx
}

block_vector_idx :: proc(bv: Block_Vector, block, idx: int) -> int {
    return bv.offsets[block] + idx
}

// returns both the index of the sparse matrix, and the index of the individual entry within it
block_sparse_idx :: proc(bs: Block_Sparse_Matrix, block_row, block_col, row, col: int) -> (int, int) {
    sm_idx := block_row * bs.block_columns + block_col

    sm := bs.values[sm_idx]
    assert(sm != nil, "Requested sparse block is structurally 0")

    return sm_idx, sparse_idx(sm^, row, col)
}

block_dense_idx :: proc(bm: Block_Dense_Matrix, block_row, block_col, row, col: int) -> int {
    block := bm.blocks[block_row * bm.block_columns + block_col]
    return block.offset + (row * block.columns + col)
}

idx :: proc{
    dense_idx,
    sparse_idx,
    vector_idx,
    block_vector_idx,
    block_sparse_idx,
    block_dense_idx,
}

// Arithmetic

dot_vec :: proc(x, y: Vector) -> (result: f64) #no_bounds_check {
	assert(len(x) == len(y))
	for i := 0; i < len(x); i += 1 {result += x[i] * y[i]}
	return result
}

scal_vec :: proc(x: Vector, a: f64) #no_bounds_check {
	if a == 0 {mem.zero(raw_data(x), len(x) * size_of(f64)); return}
	for i := 0; i < len(x); i += 1 {x[i] *= a}
}

axpy_vec :: proc(x, y: Vector, a: f64) #no_bounds_check {
	assert(len(x) == len(y))
	for i := 0; i < len(x); i += 1 {y[i] += x[i] * a}
}

axpby_vec :: proc(x, y: Vector, a, b: f64) #no_bounds_check {
	assert(len(x) == len(y))
	for i := 0; i < len(x); i += 1 {y[i] = x[i] * a + (y[i] * b)}
}

nrm2_vec :: proc(x: Vector) -> f64 {
	return math.sqrt(dot(x, x))
}

swap_vec :: slice.swap_between

dot_block_vec :: proc(x, y: Block_Vector) -> (result: f64) #no_bounds_check {
    return dot_vec(x.values, y.values)
}

scal_block_vec :: proc(x: Block_Vector, a: f64) #no_bounds_check {
    scal_vec(x.values, a)
}

axpy_block_vec :: proc(x, y: Block_Vector, a: f64) #no_bounds_check {
    axpy_vec(x.values, y.values, a)
}

axpby_block_vec :: proc(x, y: Block_Vector, a, b: f64) #no_bounds_check {
    axpby_vec(x.values, y.values, a, b)
}

nrm2_block_vec :: proc(x: Block_Vector) -> f64 {
    return nrm2_vec(x.values)
}

swap_block_vec :: proc(a: Block_Vector, b: Block_Vector) {
    swap_vec(a.values, b.values)
}

scal_dense_mat :: proc(dm: Dense_Matrix, a: f64) { scal_vec(dm.values, a) }
scal_sparse_mat :: proc(sm: Sparse_Matrix, a: f64) {scal_vec(sm.values, a) }
scal_block_dense_mat :: proc(bm: Block_Dense_Matrix, a: f64) {scal_vec(bm.values, a)}

scal_block_sparse_mat :: proc(bs: Block_Sparse_Matrix, a: f64) {
    for sm in bs.values {
        if sm == nil {continue}
        scal_sparse_mat(sm^, a)
    }
}

dot :: proc{
    dot_vec,
    dot_block_vec,
}

scal :: proc{
    scal_vec,
    scal_block_vec,
    scal_dense_mat,
    scal_sparse_mat,
    scal_block_dense_mat,
    scal_block_sparse_mat,
}

axpy :: proc{
    axpy_vec,
    axpy_block_vec,
}

axpby :: proc{
    axpby_vec,
    axpby_block_vec,
}

nrm2 :: proc{
    nrm2_vec,
    nrm2_block_vec,
}

swap :: proc{
    swap_vec,
    swap_block_vec,
}

gemv_dense :: proc(A: Dense_Matrix, x, y: Vector, a := 1.0, b := 0.0) {
	assert(A.columns == len(x))
	assert(A.rows == len(y))

	scal(y, b)
	for i in 0..<A.rows {
		for j in 0..<A.columns {
		  y[i] += (a * x[j]) * A.values[idx(A, i, j)]
		}
	}
}

gemv_block_dense :: proc(A: Block_Dense_Matrix, x, y: Block_Vector, a := 1.0, b := 0.0) {
    scal(y, b)
    for br in 0..<A.block_rows {
        y_block := block_vector_view(y, br)
        for bc in 0..<A.block_columns {
            A_block := block_dense_view(A, br, bc)
            gemv_dense(A_block, block_vector_view(x, bc), y_block, a, 1)
        }
    }
}

gemv_sparse :: proc(A: Sparse_Matrix, x, y: Vector, a := 1.0, b := 0.0) {
    assert(len(y) == sparse_mat_rows(A))

	scal(y, b)
	for row := 0; row < sparse_mat_rows(A); row += 1 {
		for idx := A.row_ptrs[row]; idx < A.row_ptrs[row + 1]; idx += 1 {
			col := A.columns[idx]
			y[row] += a * A.values[idx] * x[col]
		}
	}
	return
}

gemv_block_sparse :: proc(A: Block_Sparse_Matrix, x, y: Block_Vector, a := 1.0, b := 0.0) {
   scal(y, b)
   for br in 0..<A.block_rows {
        y_block := block_vector_view(y, br)
        for bc in 0..<A.block_columns {
            A_block := block_sparse_view(A, br, bc)
            if A_block == nil{continue}
            gemv_sparse(A_block^, block_vector_view(x, bc), y_block, a, 1)
        }
    }
}

gemv :: proc {
    gemv_dense,
    gemv_sparse,
    gemv_block_dense,
    gemv_block_sparse,
}