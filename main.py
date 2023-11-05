from Bio import SeqIO
import argparse


def needleman_wunsch(sequence1, sequence2, match=1, mismatch=-1, gap_penalty=-2):
    # inicjalizacja
    rows, cols = len(sequence1) + 1, len(sequence2) + 1
    score_matrix = [[0 for _ in range(cols)] for _ in range(rows)]

    # wypełnianie macierzy
    for i in range(1, rows):
        for j in range(1, cols):
            match_mismatch = score_matrix[i - 1][j - 1] + (match if sequence1[i - 1] == sequence2[j - 1] else mismatch)
            gap_in_seq1 = score_matrix[i - 1][j] + gap_penalty
            gap_in_seq2 = score_matrix[i][j - 1] + gap_penalty

            score_matrix[i][j] = max(match_mismatch, gap_in_seq1, gap_in_seq2)

    aligned_seq1, aligned_seq2 = "", ""
    i, j = rows - 1, cols - 1

    while i > 0 or j > 0:
        if i > 0 and score_matrix[i][j] == score_matrix[i - 1][j] + gap_penalty:
            aligned_seq1 = sequence1[i - 1] + aligned_seq1
            aligned_seq2 = "-" + aligned_seq2
            i -= 1
        elif j > 0 and score_matrix[i][j] == score_matrix[i][j - 1] + gap_penalty:
            aligned_seq1 = "-" + aligned_seq1
            aligned_seq2 = sequence2[j - 1] + aligned_seq2
            j -= 1
        else:
            aligned_seq1 = sequence1[i - 1] + aligned_seq1
            aligned_seq2 = sequence2[j - 1] + aligned_seq2
            i -= 1
            j -= 1

    return aligned_seq1, aligned_seq2, score_matrix[rows - 1][cols - 1]


def read_fasta_file(file_path):
    sequences = []
    with open(file_path, 'r') as file:
        for record in SeqIO.parse(file, 'fasta'):
            sequences.append(str(record.seq))
    return sequences


def main():
    parser = argparse.ArgumentParser(
        description='Perform Needleman-Wunsch alignment on two amino acid sequences from a FASTA file.')
    parser.add_argument('file_path', metavar='file_path', type=str,
                        help='Path to the FASTA file containing the sequences')
    args = parser.parse_args()

    sequences = read_fasta_file(args.file_path)
    if len(sequences) != 2:
        print("Error: Input file should contain exactly two sequences.")
        return

    sequence1, sequence2 = sequences
    aligned_seq1, aligned_seq2, alignment_score = needleman_wunsch(sequence1, sequence2)

    print("Optymalne dopasowanie 1:", aligned_seq1)
    print("Optymalne dopasowanie 2:", aligned_seq2)
    print("Alignment Score:", alignment_score)


if __name__ == "__main__":
    main()
