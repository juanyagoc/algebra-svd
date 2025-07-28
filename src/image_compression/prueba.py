# src/image_copmression/qr_decomposition

def qr_decomposition(a):
    import math

    m = len(a) # Number of rows
    n = len(a[0]) # Numbre of columns

    if not all(len(row) == n for row in a):
        raise ValueError("Rows must have same number of columns")

    a_copy = [row[:] for row in a]
    q = [[0.0]*n for _ in range(m)]
    r = [[0.0]*n for _ in range(n)]

    for k in range(n):
        v = [a_copy[i][k] for i in range(m)]
        for j in range(k):
            qj = [q[i][j] for i in range(m)]
            int_r = sum(qj[i] * v[i] for i in range(m))
            r[j][k] = int_r

            for i in range(m):
                v[i] -= int_r * qj[i]

        norm = math.sqrt(sum(v_i**2 for v_i in v))

        if norm == 0:
            raise ValueError("Linearly Independent Column")

        for i in range(m):
            q[i][k] = v[i] / norm

        r[k][k] = norm

    return q, r
