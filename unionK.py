import sys
import os


def getXCount( small, pos, limit ):
    char = small[pos]
    x_count = 1
    pos_count = 1
    while pos < len(small):
        if small[pos] == 'X':
            x_count += 1
            pos_count += 1
            pos += 1
        elif small[pos] == '-':
            pos_count += 1
            pos += 1
        else:
            break
    return (x_count, pos_count)
    

def getIntersection( small, large, small_limit ):
    
    # Generates the intersection of the small and large counts
    result = ""
    i = 0
    while i < (len(large)):
        if i < len(small):
            small_char = small[i]
        else:
            small_char = '.'
        large_char = large[i]

        if small_char == 'X':
            cnt = getXCount( small, i, small_limit )
            if cnt[0] < small_limit:
                result += small[i: i + cnt[1] - 1]
                i += cnt[1] - 1 
            else:
                result += large[i : i + cnt[1] - 1]
                i += cnt[1] - 1
        else:
            result += (large_char)
            i += 1
    return result

    





file1 = sys.argv[1]
file2 = sys.argv[2]
small_max = int(sys.argv[3])
small_max_mult = int(sys.argv[4]);

f1 = open(file1, "r")
f2 = open(file2, "r")
output = open("OUTPUT", "a")

lines1 = f1.readlines()
lines2 = f2.readlines()

if len(lines1) != len(lines2):
    print("Error: Alignment length is not the same")
    sys.exit()

for i in range(len(lines1)):
    line1 = lines1[i]
    line2 = lines2[i]
    if ">" in line1:
        output.write(line1)
        continue
    if "X" in line1 or "X" in line2:
        output.write(getIntersection(line1, line2, small_max_mult*small_max))
        continue
    output.write(line1)



f1.close()
f2.close()
output.close()
