import time as t
import matplotlib.pyplot as plt

# Values to be used for score calculation
penalty_gap = -1
award_match = 1
penalty_mismatch = -1


# The helper function to print out the score matrices
def print_scores_matrix(matrix):
    for x in matrix:
        print(*x, sep="\t")


# The helper function to print out the matrices
def print_matrix(ref_matrix):
    # Looping over all the rows
    for i in range(0, len(ref_matrix)):
        print("[")
        # Looping over each column present in row 'i'
        for j in range(0, len(ref_matrix[i])):
            # Printing out the values contained in (row,column) = (i,j)
            print(ref_matrix[i][j])
            # Adding a tab if not in the last column
            if j != len(ref_matrix[i]) - 1:
                print("\t")
        print("]\n")


# Function to develop a matrix of zeros
def zeros(rows, cols):
    # Definition of an empty list
    retrieve_value = []
    # Filling up the rows of the matrix
    for i in range(rows):
        # Appending an empty list for each row
        retrieve_value.append([])
        # Filling up the columns for each row
        for j in range(cols):
            # Appending zero to each column in each row
            retrieve_value[-1].append(0)
    # Returning the matrix of zeros
    return retrieve_value


# Function for obtaining the score between any two bases in the given alignment
def score_match(input1, input2):
    if input1 == input2:
        return award_match
    elif input1 == '-' or input2 == '-':
        return penalty_gap
    else:
        return penalty_mismatch


def needleman_wunsch(testDNA, refDNA):
    # Storing the length of two input sequences
    len_test = len(testDNA)
    len_ref = len(refDNA)

    # Generating matrix of zeros to store scores
    score = zeros(len_ref + 1, len_test + 1)

    # Calculate score table

    # Filling out the first column
    for i in range(0, len_ref + 1):
        score[i][0] = penalty_gap * i

    # Filling out the first row
    for j in range(0, len_test + 1):
        score[0][j] = penalty_gap * j

    # Filling out all other values for the score matrix
    for i in range(1, len_ref + 1):
        for j in range(1, len_test + 1):
            # Calculating the score by comparing the top, left, and diagonal cells
            matches = score[i - 1][j - 1] + score_match(testDNA[j - 1], refDNA[i - 1])
            deleted = score[i - 1][j] + penalty_gap
            inserted = score[i][j - 1] + penalty_gap
            # Record the maximum score from the three possible scores calculated above
            score[i][j] = max(matches, deleted, inserted)

    # Tracebacking and computing the alignment

    # Defining the variables to store alignment
    align_storage1 = ""
    align_storage2 = ""

    # Starting from the bottom right cell in matrix
    i = len_ref
    j = len_test

    # Using 'i' and 'j' for keeping a track of where we are in the matrix
    while i > 0 and j > 0:
        score_current = score[i][j]
        score_diagonal = score[i - 1][j - 1]
        score_up = score[i][j - 1]
        score_left = score[i - 1][j]

        # Check to find out which cell the current score was calculated from
        # Eventually update i and j corresponding to that cell.
        if score_current == score_diagonal + score_match(testDNA[j - 1], refDNA[i - 1]):
            align_storage1 += testDNA[j - 1]
            align_storage2 += refDNA[i - 1]
            i -= 1
            j -= 1
        elif score_current == score_up + penalty_gap:
            align_storage1 += testDNA[j - 1]
            align_storage2 += '-'
            j -= 1
        elif score_current == score_left + penalty_gap:
            align_storage1 += '-'
            align_storage2 += refDNA[i - 1]
            i -= 1

    # Finish the tracing up to the top left cell
    while j > 0:
        align_storage1 += testDNA[j - 1]
        align_storage2 += '-'
        j -= 1
    while i > 0:
        align_storage1 += '-'
        align_storage2 += refDNA[i - 1]
        i -= 1

    # Since the score matrix has been traversed from the bottom right, the two sequences will be reversed.
    # Reversing the order of the characters in each sequence
    align_storage1 = align_storage1[::-1]
    align_storage2 = align_storage2[::-1]

    return (align_storage1, align_storage2, score)


# Find longest common subsequence (LCS) of two strings
# Dynamic programming with top-down approach
def lcs(seq1, seq2):
    # Return c where c[i][j] contains length of LCS of u[i:] and v[j:]
    table = [[-1] * (len(seq2) + 1) for _ in range(len(seq1) + 1)]
    lcs_helper(seq1, seq2, table, 0, 0)
    return table


def lcs_helper(seq1, seq2, table, i, j):
    # Return the length of LCS of u[i:] and v[j:] and fill in table
    # table[i][j] contains the length of LCS of seq1[i:] and seq2[j:].
    # The function fills in table
    if table[i][j] >= 0:
        return table[i][j]

    if i == len(seq1) or j == len(seq2):
        q = 0
    else:
        if seq1[i] == seq2[j]:
            q = 1 + lcs_helper(seq1, seq2, table, i + 1, j + 1)
        else:
            q = max(lcs_helper(seq1, seq2, table, i + 1, j),
                    lcs_helper(seq1, seq2, table, i, j + 1))
    table[i][j] = q
    return q


def print_lcs(seq1, seq2, table):
    # Print one LCS of se1 and seq2 seq1sing table
    i = j = 0
    tolist = []
    while not (i == len(seq1) or j == len(seq2)):
        if seq1[i] == seq2[j]:
            tolist.append(seq1[i].capitalize())
            i += 1
            j += 1
        elif table[i][j + 1] > table[i + 1][j]:
            j += 1
        else:
            i += 1
    return tolist


def string_compare(s1, s2):
    count = 0
    gap_count = 0
    for i, j in zip(s1, s2):
        if i == '-' or j == '-':
            gap_count += 1
        if i == j:
            count += 1
    return count, gap_count

# Reference for Hamming and Edit Distance: https://nbviewer.jupyter.org/github/BenLangmead/comp-genomics-class/blob/master/notebooks/CG_DP_EditDist.ipynb 
def hamming_distance(x, y):
    ''' Return Hamming distance between x and y '''
    assert len(x) == len(y)
    nmm = 0
    for i in range(0, len(x)):
        if x[i] != y[i]:
            nmm += 1
    return nmm


def bound_edit_distance(x, y):
    ''' Return loose lower and upper bounds on the edit distance
        between x and y in O(|x| + |y|) time. '''
    if x == y: return 0, 0
    if len(x) == 0: return len(y), len(y)
    if len(y) == 0: return len(x), len(x)
    diff = abs(len(x) - len(y))
    lower = diff
    if lower == 0 and x != y:
        lower = 1
    minlen = min(len(x), len(y))
    upper = hamming_distance(x[:minlen], y[:minlen]) + diff
    return lower, upper

if __name__ == '__main__':
    list_test = []
    with open("testDNA.txt", 'r') as f1:
        for line in f1.readlines():
            list_test.append(line)

    list_ref = []
    with open("refDNA.txt", 'r') as f2:
        for line in f2.readlines():
            list_ref.append(line)

    end_time = []
    gap_data = []
    similarity_data = []

    for test, ref in zip(list_test, list_ref):
        testDNA = test
        refDNA = ref
        print("#" + "=" * 80 + "#")
        print("Sequence 1 : ", testDNA)
        print("Sequence 2 : ", refDNA)
        print("#" + "=" * 80 + "#")
        start_time = t.time()
        output1, output2, score = needleman_wunsch(testDNA, refDNA)
        end_time.append((t.time() - start_time) * 1000)
        print("Globally aligned ref-DNA sequence : " + str(output1))
        print("Globally aligned test-DNA sequence : " + str(output2))
        print("#" + "=" * 80 + "#")
        # print_scores_matrix(score)

        # print similarity and gap
        cnt, gap_cnt = string_compare(output1, output2)
        similarity = (cnt / len(output1)) * 100
        similarity_data.append(similarity)
        print("#" + "=" * 80 + "#")
        print("Similarity : " + str(similarity) + "%")
        gap = (gap_cnt / len(output1)) * 100
        gap_data.append(gap)
        print("Gap : " + str(gap) + "%")
        print("#" + "=" * 80 + "#")

        print('Longest Common Subsequence : ', end="")
        print(''.join(print_lcs(testDNA, refDNA, lcs(testDNA, refDNA))))
        print("#" + "=" * 80 + "#")

        #Print Hamming and Edit Distance
        # print("#" + "=" * 80 + "#")
        print("Hamming Distance : " + str(hamming_distance(output1, output2)))
        print("Edit Distance : " + str(bound_edit_distance(output1, output2)))
        print("#" + "=" * 80 + "#")

    name = [2, 4, 8, 16, 32, 64, 128, 256, 512]
    for n, end in zip(name, end_time):
        print(str(n) + " : " + str(end))

    # import the data from emboss into a list
    emboss_gap = []
    emboss_similarity = []
    with open("result.txt", 'r') as f:
        for line in f.readlines():
            emboss_similarity.append(float(line.split()[0].replace("%", "")))
            emboss_gap.append(float(line.split()[1].replace("%", "")))

    # plot functions
    plt.plot(name, emboss_similarity, color='red', label='EmbossSimilarity', marker='o')
    plt.plot(name, similarity_data, color='blue', label='Similarity', marker='*')
    plt.legend()
    plt.figure()
    plt.plot(name, emboss_gap, color='red', label='EmbossGap', marker='o')
    plt.plot(name, gap_data, color='blue', label='Gap', marker='*')
    plt.legend(loc='best')
    plt.show()
