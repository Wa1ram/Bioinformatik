import numpy as np
import heapq

GAP = '*'

def read_fasta(filename):
    sequences = []
    seq_names = []
    with open(filename, 'r') as f:
        seq = ''
        for line in f:
            if line.startswith('>'):
                seq_names.append(line.strip("> \t\n\r"))
                if seq:
                    sequences.append(seq)
                    seq = ''
                continue
            seq += line.strip()
        if seq:
            sequences.append(seq)
    return sequences, seq_names



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

            for col_header, val in zip(col_headers, values):
                matrix[(row_header, col_header)] = val
    return matrix, alphabet


def get_alignment_score(seq1, seq2, scoring_matrix):
    rows = len(seq1) + 1
    cols = len(seq2) + 1
    matrix = np.zeros((rows, cols), dtype=int)

    # Initialize first row and column
    for i in range(1, cols):
        matrix[0][i] = matrix[0][i-1] + scoring_matrix[(GAP, seq2[i-1])]
    for i in range(1, rows):
        matrix[i][0] = matrix[i-1][0] + scoring_matrix[(seq1[i-1], GAP)]

    # Fill the rest of the matrix
    for i in range(1, rows):
        for j in range(1, cols):
            left = matrix[i][j-1] + scoring_matrix[(GAP, seq2[j-1])]
            up = matrix[i-1][j] + scoring_matrix[(seq1[i-1], GAP)]
            diagonal = matrix[i-1][j-1] + scoring_matrix[(seq1[i-1], seq2[j-1])]
            matrix[i][j] = max(left, up, diagonal)
    score = matrix[rows-1][cols-1]
    return score



def create_similarity_matrix(sequences, scoring_matrix):
    sim = {}
    max_heap = []
    
    for i in range(len(sequences)):
        for j in range(i+1, len(sequences)):
            score = get_alignment_score(sequences[i], sequences[j], scoring_matrix)
            sim[(i, j)] = score
            sim[(j, i)] = score
            max_heap.append((-score, frozenset([i]), frozenset([j])))
            
    heapq.heapify(max_heap)  # O(k^2) Heap-construction
    return sim, max_heap
            
    
    
def compute_guide_tree(sequences, sim_mat: dict, sim_max_heap: heapq):
    active = {frozenset([i]) for i in range(len(sequences))}
    
    while active:
        # pop pair with lowest score from heap
        # lazy delete irrelevant pairs (i or j already picked)
        while True:
            _, i, j = heapq.heappop(sim_max_heap)
            if i in active and j in active:
                break
        
        active.remove(i)
        active.remove(j)
        print_inner_node(i, j)
        pair = i | j
        
        # fill new column with scores
        max_score = float('-inf')
        max_cluster = None
        for cluster in active:
            score = sum(sim_mat[(seq1, seq2)] for seq1 in cluster for seq2 in pair)
            score /= len(cluster) * len(pair)
            if score > max_score:
                max_score = score
                max_cluster = cluster
        
        # add only max score to heap - disregard all other 
        heapq.heappush(sim_max_heap, (-max_score, max_cluster, pair))
        
        # only add pair if there is another cluster left to join with
        if active:
            active.add(pair)


def print_inner_node(i, j):
    i = [x + 1 for x in i]
    j = [x + 1 for x in j]
    if max(i) > max(j):
        i, j = j, i
    print(f"({'+'.join(map(str, sorted(i)))}, {'+'.join(map(str, sorted(j)))})")
    

def print_matrix(matrix, labels):
    col_width = max(max(len(str(c)) for c in labels), 4) + 1  # at least 4, plus space

    print(" " * (col_width), end="")
    for label in labels:
        print(f"{label:>{col_width}}", end="")
        
    print()
    
    for i, row in enumerate(labels):
        print(f"{row:>{col_width}}", end="")
        for j, col in enumerate(labels):
            value = matrix.get((i, j))
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
    sequences, seq_names = read_fasta("sequences.fasta")
    scoring_matrix, alphabet = load_scoring_matrix("blosum62.txt")
    
    sim_mat, sim_max_heap = create_similarity_matrix(sequences, scoring_matrix)
    print_matrix(sim_mat, seq_names)
    
    compute_guide_tree(sequences, sim_mat, sim_max_heap)

if __name__ == "__main__":
    main()