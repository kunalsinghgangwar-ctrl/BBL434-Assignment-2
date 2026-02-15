import numpy as np
import argparse

INF = 10**9

# ---------- FASTA READER ----------
def read_fasta(file):
    seq = ""
    with open(file) as f:
        for line in f:
            if not line.startswith(">"):
                seq += line.strip().upper()
    return seq


# ---------- GOTOh AFFINE GAP GLOBAL ALIGNMENT ----------
def affine_global_alignment(seq1, seq2, match, mismatch, gap_open, gap_extend):

    n, m = len(seq1), len(seq2)

    # Three DP matrices (Gotoh algorithm)
    M = np.full((n+1, m+1), -INF)   # match/mismatch
    X = np.full((n+1, m+1), -INF)   # gap in seq2
    Y = np.full((n+1, m+1), -INF)   # gap in seq1

    # Traceback matrices
    tM = np.zeros((n+1, m+1), dtype=int)
    tX = np.zeros((n+1, m+1), dtype=int)
    tY = np.zeros((n+1, m+1), dtype=int)

    # ---------- INITIALIZATION ----------
    M[0,0] = 0

    for i in range(1, n+1):
        X[i,0] = gap_open + (i-1)*gap_extend
        tX[i,0] = 1

    for j in range(1, m+1):
        Y[0,j] = gap_open + (j-1)*gap_extend
        tY[0,j] = 2

    # ---------- FILL MATRICES ----------
    for i in range(1, n+1):
        for j in range(1, m+1):

            score = match if seq1[i-1] == seq2[j-1] else mismatch

            # M matrix
            prev = [M[i-1,j-1], X[i-1,j-1], Y[i-1,j-1]]
            M[i,j] = max(prev) + score
            tM[i,j] = np.argmax(prev)

            # X matrix (gap in seq2)
            open_gap = M[i-1,j] + gap_open + gap_extend
            extend_gap = X[i-1,j] + gap_extend
            if open_gap >= extend_gap:
                X[i,j] = open_gap
                tX[i,j] = 0
            else:
                X[i,j] = extend_gap
                tX[i,j] = 1

            # Y matrix (gap in seq1)
            open_gap = M[i,j-1] + gap_open + gap_extend
            extend_gap = Y[i,j-1] + gap_extend
            if open_gap >= extend_gap:
                Y[i,j] = open_gap
                tY[i,j] = 0
            else:
                Y[i,j] = extend_gap
                tY[i,j] = 2

    # ---------- CHOOSE BEST FINAL SCORE ----------
    matrices = [M[n,m], X[n,m], Y[n,m]]
    matrix = np.argmax(matrices)
    score = matrices[matrix]

    # ---------- TRACEBACK ----------
    align1, align2 = "", ""
    i, j = n, m

    while i > 0 or j > 0:

        if matrix == 0:   # from M
            prev = tM[i,j]
            align1 = seq1[i-1] + align1
            align2 = seq2[j-1] + align2
            i -= 1; j -= 1
            matrix = prev

        elif matrix == 1: # from X (gap in seq2)
            prev = tX[i,j]
            align1 = seq1[i-1] + align1
            align2 = "-" + align2
            i -= 1
            matrix = prev

        else:             # from Y (gap in seq1)
            prev = tY[i,j]
            align1 = "-" + align1
            align2 = seq2[j-1] + align2
            j -= 1
            matrix = prev

    return score, align1, align2


# ---------- MAIN ----------
def main():
    parser = argparse.ArgumentParser(description="Global DNA Alignment with Affine Gaps")
    parser.add_argument("fasta1")
    parser.add_argument("fasta2")
    parser.add_argument("--match", type=int, default=2)
    parser.add_argument("--mismatch", type=int, default=-1)
    parser.add_argument("--gap_open", type=int, default=-5)
    parser.add_argument("--gap_extend", type=int, default=-1)
    args = parser.parse_args()

    seq1 = read_fasta(args.fasta1)
    seq2 = read_fasta(args.fasta2)

    score, a1, a2 = affine_global_alignment(
        seq1, seq2,
        args.match, args.mismatch,
        args.gap_open, args.gap_extend
    )

    print("\n===== BEST GLOBAL ALIGNMENT =====")
    print("Score:", score)
    print(a1)
    print(a2)


if __name__ == "__main__":
    main()
