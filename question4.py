
# Global alignment from previous task
def global_alignment(seq1, seq2, scoring_function):
    # Gap penalty
    gap = -1

    n = len(seq1)
    m = len(seq2)

    # ---------------------------------------------------------
    # 1. Create and initialise the dynamic programming matrix
    # ---------------------------------------------------------

    # score[i][j] = best alignment score for
    # seq1[:i] and seq2[:j]
    score = [[0] * (m + 1) for _ in range(n + 1)]

    # Aligning a sequence with an empty sequence requires
    # inserting a gap for every character.
    for i in range(1, n + 1):
        score[i][0] = i * gap

    for j in range(1, m + 1):
        score[0][j] = j * gap

    # ---------------------------------------------------------
    # 2. Fill the dynamic programming matrix
    # ---------------------------------------------------------

    for i in range(1, n + 1):
        for j in range(1, m + 1):

            # Option 1: align the two characters
            diagonal = (
                score[i - 1][j - 1]
                + scoring_function(seq1[i - 1], seq2[j - 1])
            )

            # Option 2: align seq1[i-1] with a gap
            up = score[i - 1][j] + gap

            # Option 3: align seq2[j-1] with a gap
            left = score[i][j - 1] + gap

            # Keep the best possible score
            score[i][j] = max(diagonal, up, left)

    # ---------------------------------------------------------
    # 3. Traceback
    # ---------------------------------------------------------

    aligned_seq1 = []
    aligned_seq2 = []

    i = n
    j = m

    while i > 0 or j > 0:

        # Diagonal: seq1[i-1] aligned with seq2[j-1]
        if (
            i > 0
            and j > 0
            and score[i][j]
            == score[i - 1][j - 1]
            + scoring_function(seq1[i - 1], seq2[j - 1])
        ):
            aligned_seq1.append(seq1[i - 1])
            aligned_seq2.append(seq2[j - 1])

            i -= 1
            j -= 1

        # Up: seq1[i-1] aligned with a gap
        elif (
            i > 0
            and score[i][j] == score[i - 1][j] + gap
        ):
            aligned_seq1.append(seq1[i - 1])
            aligned_seq2.append("-")

            i -= 1

        # Left: seq2[j-1] aligned with a gap
        else:
            aligned_seq1.append("-")
            aligned_seq2.append(seq2[j - 1])

            j -= 1

    aligned_seq1 = "".join(reversed(aligned_seq1))
    aligned_seq2 = "".join(reversed(aligned_seq2))

    return aligned_seq1, aligned_seq2, float(score[n][m])