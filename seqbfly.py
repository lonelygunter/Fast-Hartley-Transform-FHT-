import math
import numpy as np
import sys

def cosN(n, l):
    return math.cos((2 * math.pi * n) / (2 ** l))

def sinN(n, l):
    return math.sin((2 * math.pi * n) / (2 ** l))

def transform_array(data):
    n = len(data)
    
    if n == 2:
        return data
    
    temp = []
    
    # Separate even and odd index elements
    even_indices = data[0:n-1:2]  # elements at even positions (0, 2, 4, ...)
    odd_indices = data[1:n:2]  # elements at odd positions (1, 3, 5, ...)
    
    temp += transform_array(even_indices)
    temp += transform_array(odd_indices)
    
    return temp

# Check if the number of command-line arguments is correct
if len(sys.argv) != 2:
    print("Usage: python script.py <number>")


# Take the number from the command line and convert it to an integer
N = int(sys.argv[1])

# Create a list of that length with default values (e.g., 0)
data = [0] * N

for i in range(N):
    data[i] = i+1

print(f"Data: \t{np.round(data, 2)}")

logN = int(math.log2(N))


# Stage 0
# for k in range(1, N, 2):
#     k2 = int(k+((N/2)-1))
#     if k2 < N:
#         temp = data[k]
#         data[k] = data[k2]
#         data[k2] = temp
data = transform_array(data)

print(f"f0: \t{np.round(data, 2)}")


# Stages
temp = data.copy()
for s in range(logN):
    m = 2 ** s
    u = m
    for k in range(0, N, 2 * m):

        for j in range(m):
            t_d = k+j
            u1_d = k+j+m
            u2_d = u+k
            t = data[k + j]
            u1 = data[k + j + m]
            u2 = data[u+k]

            # print(f"data[{k + j}] \t= [{k+j}] \t+ ([{k+j+m}] * cosN({j}, {s+1})) \t+ ([{u+k}] * sinN({j}, {s+1}))")
            # print(f"data[{k + j + m}] \t= [{k+j}] \t+ ([{k+j+m}] * cosN({j+m}, {s+1})) \t+ ([{u+k}] * sinN({j+m}, {s+1}))")
            
            # print(f"data[{k + j}] \t= {round(t, 2)} \t+ ({round(u1, 2)} * cosN({j}, {s+1})) \t+ ({round(u2, 2)} * sinN({j}, {s+1}))")
            # print(f"data[{k + j + m}] \t= {round(t, 2)} \t+ ({round(u1, 2)} * cosN({j+m}, {s+1})) \t+ ({round(u2, 2)} * sinN({j+m}, {s+1}))")

            temp[k + j] = t + (u1 * cosN(j, s+1)) + (u2 * sinN(j, s+1))
            temp[k + j + m] = t + (u1 * cosN(j+m, s+1)) + (u2 * sinN(j+m, s+1))

            if u == m:
                u = pow(2, s+1)-1
            else:
                u -= 1
    data = temp.copy()

    print(f"f{s+1}: \t{np.round(data, 2)}")


# H(v)
for k in range(N):
    data[k] = data[k] / N

print(f"H(v): \t{np.round(data, 2)}")