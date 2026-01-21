// SPDX-FileCopyrightText: 2025 Rowan Apps
// SPDX-License-Identifier: MIT
/*
 Data & functions for dealing with dense matrices & dense vectors.
*/
package la

import "core:math"
import "core:mem"
import "core:slice"

// TODO: explicit SIMD for math ops.
// TODO: is gaussian elimination really the best bet (I have my doubts)

Vec :: []f64

Dense_Matrix :: struct {
	values:        [^]f64,
	rows, columns: i32, //i32 is the smallest integer we can do here because of padding.
}

Dense_LU_Factorized :: struct {
	mat:    Dense_Matrix,
	pivots: []i16,
}

dense_create :: proc(#any_int rows, cols: i32, allocator := context.allocator) -> Dense_Matrix {
	return {
		rows = rows,
		columns = cols,
		values = make([^]f64, rows * cols, allocator)
	}
}

dense_mat_idx :: #force_inline proc(mat: Dense_Matrix, #any_int row, col: i32) -> int {
	assert(row < mat.rows && col < mat.columns)
	return int(col * mat.rows + row)
}


//--Dense vector operations--

// x . y
dot :: proc(x, y: Vec) -> (result: f64) #no_bounds_check {
	assert(len(x) == len(y))
	for i := 0; i < len(x); i += 1 {result += x[i] * y[i]}
	return result
}

// x = a * x
scal :: proc(x: Vec, a: f64) #no_bounds_check {
	if a == 0 {mem.zero(raw_data(x), len(x) * size_of(f64)); return}
	for i := 0; i < len(x); i += 1 {x[i] *= a}
}

// y = a * x + y
axpy :: proc(x, y: Vec, a: f64) #no_bounds_check {
	assert(len(x) == len(y))
	for i := 0; i < len(x); i += 1 {y[i] += x[i] * a}
}

// y = a * x + b * y
axpby :: proc(x, y: Vec, a, b: f64) #no_bounds_check {
	assert(len(x) == len(y))
	for i := 0; i < len(x); i += 1 {y[i] = x[i] * a + (y[i] * b)}
}

// r = ||x||
nrm2 :: proc(x: Vec) -> f64 {
	return math.sqrt(dot(x, x))
}

// x <=> y
swap :: slice.swap_between


//--Dense matrix operations--

dense_scal :: proc(A: Dense_Matrix, a: f64) #no_bounds_check {
	if a == 0 {mem.zero(A.values, int(A.rows * A.columns) * size_of(f64)); return}
	for i: i32 = 0; i < A.rows * A.columns; i += 1 {A.values[i] *= a}
}

// y = a * A * x + y * b
dense_mv :: proc(A: Dense_Matrix, x, y: Vec, a := 1.0, b := 0.0) #no_bounds_check {
	assert(A.columns == i32(len(x)))
	assert(A.rows == i32(len(y)))

	scal(y, b)
	for j: i32 = 0; j < A.columns; j += 1 {
		for i: i32 = 0; i < A.rows; i += 1 {y[i] += (a * x[j]) * A.values[dense_mat_idx(A, i, j)]}
	}
}

// C = a * A * B + b * C
dense_mm :: proc(A, B, C: Dense_Matrix, a := 1.0, b := 0.0) #no_bounds_check {
	assert(A.columns == B.rows)
	assert(A.rows == C.rows)
	assert(B.columns == C.columns)

	dense_scal(C, b)
	for k: i32 = 0; k < B.columns; k += 1 {
		for j: i32 = 0; j < A.columns; j += 1 {
			b_jk := B.values[dense_mat_idx(B, j, k)]
			for i: i32 = 0; i < A.rows; i += 1 {
				C.values[dense_mat_idx(C, i, k)] += a * A.values[dense_mat_idx(A, i, j)] * b_jk
			}
		}
	}
}

// A = lu(A)
dense_lu_factorize :: proc(A: Dense_Matrix, pivots: []i16) -> Dense_LU_Factorized {
	assert(A.rows == A.columns)
	assert(len(pivots) == int(A.rows))

	lu := Dense_LU_Factorized{A, pivots}

	//guassian elimination with partial pivoting (not optimal)
	for k: i32 = 0; k < lu.mat.rows; k += 1 {
		pivot_row := k
		max_val := abs(lu.mat.values[dense_mat_idx(lu.mat, k, k)])
		for i: i32 = k + 1; i < lu.mat.rows; i += 1 {
			val := abs(lu.mat.values[dense_mat_idx(lu.mat, i, k)])
			if val > max_val {
				max_val = val
				pivot_row = i
			}
		}

		pivots[k] = i16(pivot_row)
		if pivot_row != k {
			for j: i32 = 0; j < lu.mat.rows; j += 1 {
				k_idx := dense_mat_idx(lu.mat, k, j)
				p_idx := dense_mat_idx(lu.mat, pivot_row, j)
				lu.mat.values[k_idx], lu.mat.values[p_idx] = lu.mat.values[p_idx], lu.mat.values[k_idx]
			}
		}

		pivot_val := lu.mat.values[dense_mat_idx(lu.mat, k, k)]
		assert(pivot_val != 0, "LU Failed: Encountered Singular Matrix.")
		for i: i32 = k + 1; i < lu.mat.rows; i += 1 {
			multiplier := lu.mat.values[dense_mat_idx(lu.mat, i, k)] / pivot_val
			lu.mat.values[dense_mat_idx(lu.mat, i, k)] = multiplier
			for j: i32 = k + 1; j < lu.mat.rows; j += 1 {
				lu.mat.values[dense_mat_idx(lu.mat, i, j)] -= multiplier * lu.mat.values[dense_mat_idx(lu.mat, k, j)]
			}
		}
	}

	return lu
}

// Solves A * x = b using LU factorization of A
// b is overwritten with the solution.
dense_lu_solve :: proc(lu: ^Dense_LU_Factorized, b: Vec) {
	assert(lu.mat.rows == i32(len(b)))

	for i: i32 = 0; i < lu.mat.rows; i += 1 {
		pivot_row := int(lu.pivots[i])
		if pivot_row != int(i) {b[i], b[pivot_row] = b[pivot_row], b[i]}
	}

	for i: i32 = 0; i < lu.mat.rows; i += 1 {
		for j: i32 = 0; j < i; j += 1 {b[i] -= lu.mat.values[dense_mat_idx(lu.mat, i, j)] * b[j]}
	}

	for i := lu.mat.rows - 1; i >= 0; i -= 1 {
		for j: i32 = i + 1; j < lu.mat.rows; j += 1 {b[i] -= lu.mat.values[dense_mat_idx(lu.mat, i, j)] * b[j]}
		b[i] /= lu.mat.values[dense_mat_idx(lu.mat, i, i)]
	}
}

// Solves A * X = B using LU factorization of A
// X and B are matrices, representing multiple right hand sides.
dense_lu_solve_multiple :: proc(lu: ^Dense_LU_Factorized, B: Dense_Matrix) {
	assert(lu.mat.rows == B.rows)
	for col: i32 = 0; col < B.columns; col += 1 {dense_lu_solve(lu, B.values[col * B.rows:][:B.rows])}
}
