import sys
from qr_decomposition import qr_decomposition

def print_matrix(name, matrix):
    print(f"Matrix {name}:")
    for row in matrix:
        print(" ".join(f"{val:10.4f}" for val in row))
    print()

def main():
    A = [
        [12.0, -51.0,  4.0],
        [6.0,  167.0, -68.0],
        [-4.0, 24.0,  -41.0],
    ]

    print("Input Matrix A:")
    print_matrix("A", A)

    try:
        Q, R = qr_decomposition(A)

        print_matrix("Q", Q)
        print_matrix("R", R)

    except ValueError as e:
        print(f"Error durante la descomposición QR: {e}", file=sys.stderr)

if __name__ == "__main__":
    main()