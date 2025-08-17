import sys
import numpy as np
from prueba import qr_decomposition


def multiply_matrices(A, B):
    m = len(A)  # Filas de A
    n = len(B[0])  # Columnas de B
    l = len(A[0])  # Columnas de A / Filas de B

    C = [[0.0 for _ in range(n)] for _ in range(m)]

    for i in range(m):
        for j in range(n):
            for k in range(l):
                C[i][j] += A[i][k] * B[k][j]

    return C


def transpose_matrix(A):
    m = len(A)
    n = len(A[0])

    AT = [[0.0 for _ in range(m)] for _ in range(n)]

    for i in range(m):
        for j in range(n):
            AT[j][i] = A[i][j]

    return AT


def print_matrix(name, matrix):
    print(f"Matrix {name}:")
    for row in matrix:
        print(" ".join(f"{val:10.4f}" for val in row))
    print()


def diagonalize_matrix(A, num_iterations=10):
    AT = transpose_matrix(A)
    M = multiply_matrices(A, AT)

    print_matrix("Initial M", M)

    for iterat in range(num_iterations):
        Q, R = qr_decomposition(M)

        M = multiply_matrices(R, Q)

        print(f"\nIteration {iterat + 1}:")
        print_matrix(f"M (Step {iterat + 1})", M)

    print("\nDiagonal matrix approximation (final M):")
    print_matrix("Final M", M)

def main2():
    A = [
        [12.0, -51.0, 4.0],
        [6.0, 167.0, -68.0],
        [-4.0, 24.0, -41.0],
    ]

    A2 = [
         [3, 2, 2],
         [2, 3, -2],
    ]

    try:
        Q, R = np.linalg.qr(A2, mode='reduced')
        print(Q, R, "\n")


    except ValueError as e:
            print(f"Error durante la descomposición QR: {e}", file=sys.stderr)

def main():
    # Matriz de ejemplo 3x2
    A = np.array([[3, 2, 2],
                  [2, 3, -2],], dtype=float)

    # SVD
    U, sigma, Vt = np.linalg.svd(A, full_matrices=True)

    Sigma = np.zeros((A.shape[0], A.shape[1]))
    np.fill_diagonal(Sigma, sigma)

    print("\nMatriz U:")
    print(U)
    print("\nMatriz Sigma:")
    print(Sigma)
    print("\nMatriz V (no transpuesta):")
    print(Vt.T)

if __name__ == "__main__":
    main()