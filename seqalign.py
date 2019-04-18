# Values to be used for score calculation
penalty_gap = -1
award_match = 1
penalty_mismatch = -1

# Make a score matrix with these two sequences
testDNA = "ccggggatgcgtgatgtttaagagttaaagcctggcctcccaaccgacat"
refDNA = "agagcccatgttgtgcggtgtcgccatttcgggggccaaa"
len_test = len(testDNA)
len_ref = len(refDNA)


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
    retrive_value = []
    # Filling up the rows of the matrix
    for i in range(rows):
        # Appending an empty list for each row
        retrive_value.append([])
        # Filling up the columns for each row
        for j in range(cols):
            # Appending zero to each column in each row
            retrive_value[-1].append(0)
    # Returning the matrix of zeros
    return retrive_value


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
    list = []
    while not (i == len(seq1) or j == len(seq2)):
        if seq1[i] == seq2[j]:
            list.append(seq1[i].capitalize())
            i += 1
            j += 1
        elif table[i][j + 1] > table[i + 1][j]:
            j += 1
        else:
            i += 1
    return list


def string_compare(s1, s2):
    count = 0
    gap_count = 0
    for i, j in zip(s1, s2):
        if i == '-' or j == '-':
            gap_count += 1
        if i == j:
            count += 1
    return count, gap_count


# testDNA = ""
print("#==============================================================================#")
print("Sequence 1 : ", testDNA)
print("Sequence 2 : ", refDNA)
print("#==============================================================================#")
output1, output2, score = needleman_wunsch(testDNA, refDNA)
print("Globally aligned ref-DNA sequence : " + str(output1))
print("Globally aligned test-DNA sequence : " + str(output2))
print("#==============================================================================#")
print_scores_matrix(score)
cnt, gap_cnt = string_compare(output1, output2)
similarity = cnt / len(output1)
print("#==============================================================================#")
print("Similarity : " + str(similarity * 100) + "%")
gap = gap_cnt / len(output1)
print("Gap : " + str(gap * 100) + "%")
print("#==============================================================================#")

c = lcs(testDNA, refDNA)
print('Longest Common Subsequence: ', end="")
list = print_lcs(testDNA, refDNA, c)
str_list = ''.join(list)
print(str_list)
print("#==============================================================================#")

# print(output1 + "\n" + output2)
