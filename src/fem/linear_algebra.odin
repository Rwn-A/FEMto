// Linear algebra data structures and utilites.
// Not exhaustive, only operations needed up to now have been implemented.
// Includes voigt form structures for elasticity type weak forms.
package fem

import "core:math"
import "core:math/linalg"
import "core:slice"

Vec3 :: linalg.Vector3f64
Mat3 :: linalg.Matrix3x3f64

Voigt6 :: [6]f64
Voigt6x6 :: [36]f64

Vector :: []f64

Sparse_Matrix :: struct {
	using sp: Sparsity,
	values:   []f64,
}

Dense_Matrix :: struct {
	rows, columns: int,
	values:        []f64,
}

// i32 is more then enough, since we store lots of these, worth doing to lower memory.
Sparsity :: struct {
	row_ptrs: []i32,
	columns:  []i32,
}

dense_matrix_create :: proc(rows, columns: int, allocator := context.allocator) -> Dense_Matrix {
	return {rows = rows, columns = columns, values = make([]f64, rows * columns, allocator)}
}

dense_matrix_destroy :: proc(m: Dense_Matrix, allocator := context.allocator) {
	defer delete(m.values)
}

dense_matrix_idx :: proc(dm: Dense_Matrix, row, col: int) -> int {
	return row * dm.columns + col
}

dense_matrix_gemv :: proc(A: Dense_Matrix, x, y: Vector, a := 1.0, b := 0.0) {
	assert(A.columns == len(x))
	assert(A.rows == len(y))

	for i in 0 ..< A.rows {
		sum := 0.0
		for j in 0 ..< A.columns {sum += A.values[idx(A, i, j)] * x[j]}
		y[i] = a * sum + b * y[i]
	}
}


sparse_matrix_from_sparsity :: proc(sp: Sparsity, allocator := context.allocator) -> Sparse_Matrix {
	return {sp = sp, values = make([]f64, len(sp.columns), allocator)}
}

// does not destory sparsity
sparse_matrix_destroy :: proc(sm: Sparse_Matrix, allocator := context.allocator) {
	delete(sm.values, allocator)
}

sparse_matrix_rows :: proc(sparsity: Sparsity) -> int {
	return len(sparsity.row_ptrs) - 1
}

sparse_matrix_idx :: proc(sp: Sparsity, #any_int row, col: i32) -> int {
	row_slice := sp.columns[sp.row_ptrs[row]:sp.row_ptrs[row + 1]]
	idx, found := slice.binary_search(row_slice, col)
	if !found {panic("Sparse matrix dot have the (row, col) pair in sparsity pattern.")}
	return int(sp.row_ptrs[row]) + idx
}

vector_dot :: proc(a, b: Vector) -> (r: f64) {
	assert(len(a) == len(b))
	for i in 0 ..< len(a) {
		r += a[i] * b[i]
	}
	return r
}

vector_axpy :: proc(x, y: Vector, a: f64 = 1) {
	assert(len(x) == len(y))
	for i in 0 ..< len(x) {
		y[i] += a * x[i]
	}
}

vector_scal :: proc(x: Vector, a: f64) {
	for &e in x {e *= a}
}

vector_nrm2 :: proc(x: Vector) -> f64 {
	return math.sqrt(vector_dot(x, x))
}

voigt_gemv :: proc(A: Voigt6x6, x: Voigt6, y: ^Voigt6, a: f64 = 1.0, b: f64 = 0.0) {
	for i in 0 ..< 6 {
		sum: f64 = 0
		for j in 0 ..< 6 {
			sum += A[voigt6x6_idx(i, j)] * x[j]
		}
		y[i] = a * sum + b * y[i]
	}
}

voigt_axpy :: proc(x: Voigt6, y: ^Voigt6, a: f64 = 1) {
	for i in 0 ..< 6 {
		y[i] += a * x[i]
	}
}

voigt6x6_idx :: proc(row, col: int) -> int {
	return (row * 6) + col
}

voigt6x6_scal :: proc(x: ^Voigt6x6, a: f64) {
	for &e in x {e *= a}
}

voigt6_scal :: proc(x: ^Voigt6, a: f64) {
	for &e in x {e *= a}
}

frob :: proc(A, B: Mat3) -> (r: f64) {
	#unroll for i in 0 ..< 3 {r += dot(A[i], B[i])}
	return r
}

idx :: proc {
	dense_matrix_idx,
	sparse_matrix_idx,
	voigt6x6_idx,
}


dot :: proc {
	linalg.vector_dot, // for the vec3 and voigt6
	vector_dot,
}


axpy :: proc {
	vector_axpy,
	voigt_axpy,
}

gemv :: proc {
	voigt_gemv,
	dense_matrix_gemv,
}

scal :: proc {
	vector_scal,
	voigt6x6_scal,
	voigt6_scal,
}

nrm2 :: proc {
	vector_nrm2,
}
