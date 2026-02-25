// SPDX-FileCopyrightText: 2026 Rowan Apps
// SPDX-License-Identifier: MIT
/*
 seperate from `la` for clarity.

 Uses hardcoded arrays in hopes the compiler can better unroll and vectorize.
*/
package fem

Voigt6 :: [6]f64

// row major
Voigt6x6 :: [36]f64

VOIGT6x6_MAP: [3][3]int = {{0, 5, 4}, {5, 1, 3}, {4, 3, 2}}

voigt_dot :: proc(a, b: Voigt6) -> (r: f64) {
	for i in 0 ..< 6 {r += a[i] * b[i]}
	return r
}

voigt_gemv :: proc(A: Voigt6x6, x: Voigt6, y: ^Voigt6, a, b: f64) {
	for i in 0 ..< 6 {
		sum: f64 = 0
		for j in 0 ..< 6 {
			sum += A[voigt6x6_idx(i, j)] * x[j]
		}
		y[i] = a * sum + b * y[i]
	}
}

voigt6x6_tensor_index :: proc(i, j, k, l: int) -> int {
	I := VOIGT6x6_MAP[i][j]
	J := VOIGT6x6_MAP[k][l]
	return I * 6 + J
}

voigt6x6_idx :: proc(row, col: int) -> int {
	return (row * 6) + col
}
