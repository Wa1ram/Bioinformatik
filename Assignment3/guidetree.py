import numpy as np

GAP = '*'

def read_fasta(filename):
    sequences = []
    metadata = []
    with open(filename, 'r') as f:
        seq = ''
        for line in f:
            if line.startswith('>'):
                metadata.append(line.strip("> \t\n\r"))
                if seq:
                    sequences.append(seq)
                    seq = ''
                continue
            seq += line.strip()
        if seq:
            sequences.append(seq)
    return sequences, metadata

def load_scoring_matrix(filename):
    matrix = {}
    alphabet = []
    with open(filename, 'r') as f:
        # Skip comment lines
        for line in f:
            if not line.startswith('#'):
                break
        # Parse column headers
        col_headers = line.split()
        alphabet = col_headers
        # Parse matrix rows
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            parts = line.split()
            row_header = parts[0]
            values = list(map(int, parts[1:]))
            matrix[row_header] = {}
            for col_header, val in zip(col_headers, values):
                matrix[row_header][col_header] = val
    return matrix, alphabet


def get_alignment_score(seq1, seq2, scoring_matrix):
    rows = len(seq1) + 1
    cols = len(seq2) + 1
    matrix = np.zeros((rows, cols), dtype=int)

    # Initialize first row and column
    for i in range(1, cols):
        matrix[0][i] = matrix[0][i-1] + scoring_matrix[GAP][seq2[i-1]]
    for i in range(1, rows):
        matrix[i][0] = matrix[i-1][0] + scoring_matrix[seq1[i-1]][GAP]

    # Fill the rest of the matrix
    for i in range(1, rows):
        for j in range(1, cols):
            left = matrix[i][j-1] + scoring_matrix[GAP][seq2[j-1]]
            up = matrix[i-1][j] + scoring_matrix[seq1[i-1]][GAP]
            diagonal = matrix[i-1][j-1] + scoring_matrix[seq1[i-1]][seq2[j-1]]
            matrix[i][j] = max(left, up, diagonal)
    score = matrix[rows-1][cols-1]
    return score


def create_similarity_matrix(sequences, metadata, scoring_matrix):
    sim = {}

    for i in range(len(sequences)):
        sim[metadata[i]] = {}
        for j in range(i+1, len(sequences)):
            sim[metadata[i]][metadata[j]] = get_alignment_score(sequences[i], sequences[j], scoring_matrix)
            
    return sim
            
    

def print_matrix(matrix, labels):
    col_width = max(max(len(str(c)) for c in labels), 4) + 1  # at least 4, plus space

    print(" " * (col_width), end="")
    for c in labels:
        print(f"{c:>{col_width}}", end="")
        
    print()
    
    for row_c in labels:
        print(f"{row_c:>{col_width}}", end="")
        for col_c in labels:
            value = matrix.get(row_c, {}).get(col_c)
            if value is not None:
                print(f"{value:>{col_width}d}", end="")
            else:
                print(" " * col_width, end="")
        print()
        


def print_alignment_matrix(matrix, seq1, seq2):
    print("\nAlignment Matrix:")
    print("     " + " ".join(f"{c:>4}" for c in seq2))
    for i in range(len(seq1) + 1):
        if i == 0:
            row_label = " "
        else:
            row_label = seq1[i-1]
        print(f"{row_label:>4}", end=" ")
        for j in range(len(seq2) + 1):
            print(f"{matrix[i][j]:4d}", end=" ")
        print()
        

def main():
    sequences, metadata = read_fasta("sequences.fasta")
    scoring_matrix, alphabet = load_scoring_matrix("blosum62.txt")
    
    sim = create_similarity_matrix(sequences, metadata, scoring_matrix)
    print_matrix(sim, metadata)

if __name__ == "__main__":
    main()