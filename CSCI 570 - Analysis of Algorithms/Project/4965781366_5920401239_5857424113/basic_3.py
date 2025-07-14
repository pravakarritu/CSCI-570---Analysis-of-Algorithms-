"""
Basic sequence alignment algorithm
"""
import sys
import os
import numpy as np
import time
import psutil

# To run:
# python3 basic_3.py SampleTestCases/input1.txt out/output1.txt

# DNA bases to index
Base = {'A': 0, 'C': 1, 'G': 2, 'T': 3}


def runBasic(x, y, gap, alpha):
    start_time = time.time()

    cost, str1, str2 = seqAlign(x, y, gap, alpha)
    end_time = time.time()
    total_time = (end_time - start_time) * 1000

    mem = process_memory()
    return cost, str1, str2, total_time, mem

def seqAlign(x, y, gap, alpha):
    m, n = len(x), len(y)
    dp = np.empty((m + 1, n + 1))
    # Init col 0 and row 0
    for i in range(m + 1):
        dp[i][0] = gap * i
    for j in range(n + 1):
        dp[0][j] = gap * j

    # Find min cost buttom-up
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            if x[i - 1] != y[j - 1]:
                match = alpha[Base[x[i - 1]]][Base[y[j - 1]]] + dp[i - 1][j - 1]
                x_not = gap + dp[i - 1][j]
                y_not = gap + dp[i][j - 1]
                dp[i][j] = min(match, x_not, y_not)
            else:
                dp[i][j] = dp[i - 1][j - 1]

    # Find string alignments top-down
    m_id, n_id = m, n
    align1, align2 = "", ""
    while m_id > 0 and n_id > 0:
        match = dp[m_id - 1][n_id - 1] + alpha[Base[x[m_id - 1]]][Base[y[n_id - 1]]]
        x_not = gap + dp[m_id - 1][n_id]
        y_not = gap + dp[m_id][n_id - 1]
        if dp[m_id][n_id] == match:
            align1 = x[m_id - 1] + align1
            align2 = y[n_id - 1] + align2
            m_id -= 1
            n_id -= 1
        elif dp[m_id][n_id] == x_not:
            align1 = x[m_id - 1] + align1
            align2 = "_" + align2
            m_id -= 1
        elif dp[m_id][n_id] == y_not:
            align1 = "_" + align1
            align2 = y[n_id - 1] + align2
            n_id -= 1
            
    while m_id > 0:
        align1 = x[m_id - 1] + align1
        align2 = "_" + align2
        m_id -= 1
    while n_id > 0:
        align2 = y[n_id - 1] + align2
        align1 = "_" + align1
        n_id -= 1
    return dp[m][n], align1, align2


# Sample code for time and memory calculation
def process_memory():
    process = psutil.Process()
    memory_info = process.memory_info()
    memory_consumed = int(memory_info.rss / 1024)
    return memory_consumed
    
# Generate two input strings
def generateInput(input_file):
    lines = input_file.readlines()
    s1, s2 = "", ""
    final_s1, final_s2 = None, None
    seen = False
    j, k = 0, 0
    for line in lines:
        line = line.rstrip("\n")
        if not line.isdigit():
            if seen:
                final_s2 = line
                s2 = line
                break
            s1 = line
            final_s1 = line
            seen = True
        else:
            if len(final_s1) >= int(line) + 1:
                final_s1 = final_s1[:int(line)+1] + final_s1 + final_s1[int(line)+1:]
            else:
                final_s1 = final_s1 + s1
            j += 1
    
    for line in lines[j + 2 :]:
        line = line.rstrip("\n")
        if line.isdigit():
            if len(final_s2) >= int(line) + 1:
                final_s2 = final_s2[:int(line)+1] + final_s2 + final_s2[int(line)+1:]
            else:
                final_s2 = final_s2 + s2
            k += 1

    return final_s1, final_s2



def main(argv):
    if len(argv) > 2:
        return "Invalid arguments"
    input_file = open(argv[0], "r")

    dir = "/".join(argv[1].split('/')[:-1])
    if not os.path.exists(dir):
        os.makedirs(dir)
    output_file = open(argv[1], "w+")

    # Generate two input strings
    final_s1, final_s2 = generateInput(input_file)
    
    # Hardcoded costs
    gap = 30
    alpha = np.array(
        [[0, 110, 48, 94], [110, 0, 118, 48], [48, 118, 0, 110], [94, 48, 110, 0]]
    )

    cost, align1, align2, time, mem = runBasic(final_s1, final_s2, gap, alpha)
    lines = [str(int(cost))+"\n", align1+"\n", align2+"\n", str(time)+"\n", str(mem)+"\n"]
    output_file.writelines(lines)

    input_file.close()
    output_file.close()


if __name__ == "__main__":
    main(sys.argv[1:])
