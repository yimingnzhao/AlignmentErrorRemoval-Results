import sys
import os



alignment_file = sys.argv[1]

f1 = open(alignment_file, "r")

# Ignores the first line
f1.readline();

# Gets the number of nongap characters of the first sequence
line = f1.readline();
total_chars = 0
non_gap_chars = 0;
for char in line:
    total_chars += 1;
    if (char != '-'):
        non_gap_chars += 1

f1.close()
#print('Total Characters in Line: ' + str(total_chars))
#print('Non-Gap Characters in Line: ' + str(non_gap_chars))
print(non_gap_chars)
