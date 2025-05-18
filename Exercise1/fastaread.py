import sys

def fasta_sequence_lengths(file):
    """
    Returns the length of each sequence in a FASTA file.
    """
    seq_lengths = []
    seq_count = -1
    
    for line in file:
        line = line.strip()
        if not line:
            continue
        
        # description line -> new sequence
        if line.startswith(">"):
            seq_count += 1
            seq_lengths.append(0)
        else:
            seq_lengths[seq_count] += len(line)
        
    return seq_lengths


def main(args):
    for file_name in args:
        with open(file_name) as fasta_file:
            lengths = fasta_sequence_lengths(fasta_file)
        for l in lengths:
            print(l)

if __name__ == '__main__':
    # pass all arguments (not including the invoked .py file)
    main(sys.argv[1:])