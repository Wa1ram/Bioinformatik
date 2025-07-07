import gzip
import sys


def max_cliques(edges: set):
    clique_size = 2
    current_cliques = edges

    while True:
        clique_size += 1
        prev_cliques = list(current_cliques)
        current_cliques = set()

        for i in range(len(prev_cliques)):
            for j in range(i + 1, len(prev_cliques)):
                c1 = prev_cliques[i]
                c2 = prev_cliques[j]

                candidate = c1 | c2
                if len(candidate) != clique_size:
                    continue

                unique_in_c1 = c1 - c2
                unique_in_c2 = c2 - c1

                if (unique_in_c1 | unique_in_c2 in edges):
                    current_cliques.add(candidate)

        if not current_cliques:
            break

    return prev_cliques, clique_size - 1

    
    
def get_all_edges_from_ppi_file(filepath):
    protein_to_index = {}
    index_to_protein = {}
    edges = set()
    line_count = 0

    p_idx = 0

    open_func = gzip.open if filepath.endswith('.gz') else open
    with open_func(filepath, 'rt', encoding='utf-8') as f:
        next(f)
        for line in f:
            line_count += 1
            if line_count == 100000:
                break
            p1, p2, _ = line.split()

            for p in (p1, p2):
                if p not in protein_to_index:
                    protein_to_index[p] = p_idx
                    index_to_protein[p_idx] = p
                    p_idx += 1
            
            p1_idx = protein_to_index[p1]
            p2_idx = protein_to_index[p2]

            edge = frozenset([p1_idx, p2_idx])
            if edge not in edges:
                edges.add(edge)

    
    return edges, protein_to_index, index_to_protein


def print_cliques(cliques, clique_size, idx_to_protein):
    print(f"Size of maximal clique(s): {clique_size}")
    for nr, clique in enumerate(cliques, 1):
        print(f"Clique nr. {nr}")
        print(", ".join(idx_to_protein[idx] for idx in clique))



def main(filepath_ppi_edges):
    edges, protein_to_idx, idx_to_protein = get_all_edges_from_ppi_file(filepath_ppi_edges)
    cliques, clique_size = max_cliques(edges)
    print_cliques(cliques, clique_size, idx_to_protein)


if __name__ == "__main__":
    main(sys.argv[1])