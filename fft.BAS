DIM R(1024), I(1024)
P = 3
N = 2^P

! Insert data
FOR I = 0 TO N - 1
    R(I) = I + 1
    I(I) = 0
NEXT I

GOSUB 220
END

! Bit reversal function
DEF FN B(J)
    J1 = J
    K0 = 0
    FOR Q = 1 TO P
        J2 = J1 DIV 2
        K0 = K0 * 2 + (J1 - 2 * J2)
        J1 = J2
    NEXT Q
    FN B = K0
END FN

! Subroutine Get FFT
N2 = N DIV 2
F1 = P - 1
K = 0

FOR L = 1 TO P
    FOR I = 1 TO N2
        J = FN B(K) DIV 2^F1
        A = 2 * PI * J / N
        C = COS(A)
        S = SIN(A)
        U = R(K + N2) * C + I(K + N2) * S
        V = I(K + N2) * C - R(K + N2) * S
        R(K + N2) = R(K) - U
        I(K) = I(K) + V
        K = K + 1
    NEXT I
    K = K + N2
    IF K < N THEN GOTO 32
    K = 0
    F1 = F1 - 1
    N2 = N2 DIV 2
NEXT L

FOR K = 0 TO N - 1
    I = FN B(K)
    IF I < K THEN GOTO 59
    U = R(K) / N
    V = I(K) / N
    R(K) = R(I) / N
    I(K) = I(I) / N
    R(I) = U
    I(I) = V
NEXT K
RETURN
