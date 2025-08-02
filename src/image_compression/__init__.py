import sys
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

    for iter in range(num_iterations):
        Q, R = qr_decomposition(M)

        M = multiply_matrices(R, Q)

        print(f"\nIteration {iter + 1}:")
        print_matrix(f"M (Step {iter + 1})", M)

    print("\nDiagonal matrix approximation (final M):")
    print_matrix("Final M", M)

def main():
    A = [
        [12.0, -51.0, 4.0],
        [6.0, 167.0, -68.0],
        [-4.0, 24.0, -41.0],
    ]

    try:
        diagonalize_matrix(A)

    except ValueError as e:
        print(f"Error durante la descomposición QR: {e}", file=sys.stderr)


if __name__ == "__main__":
    main()