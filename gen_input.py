import random

# Function to generate a file with a specified number of random numbers
def gen_input(filename, count):
    # Generate the random numbers
    numbers = [random.randint(1, 1000) for _ in range(count)]

    # Write the numbers to the file
    with open(filename, 'w') as f:
        f.write(f"{count}\n")
        f.write(" ".join(map(str, numbers)) + "\n")

# List of file sizes
sizes = [64, 128, 256, 512, 1024]

# Create files for each size
for size in sizes:
    gen_input(f"input{size}", size)
