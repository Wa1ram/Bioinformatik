import sys

def boyer_moore(search_file, pattern, bcr_pre, gsr_pre):
    """
    Searches for occurrences of a pattern in a sequence using the Boyer-Moore algorithm
    with a sliding window approach. Returns the number of bad character rule uses,
    the number of occurrences found, and the positions of the first 10 matches.

    Args:
        search_file: File object to search in (FASTA format, sequence lines).
        pattern: String, the pattern to search for.
        bcr_pre: Dictionary, bad character rule preprocessing table.
        gsr_pre: Dictionary, good suffix rule preprocessing table.

    Returns:
        Tuple: (num_bcr_uses, num_occurences, positions)
    """
    
    pattern_length = len(pattern)
    window_size = 10 * pattern_length
    window = fill_window(search_file, "", window_size)
    t = 0
    
    num_bcr_uses = 0
    num_occurences = 0
    positions = []

    while len(window) >= pattern_length:
        p = pattern_length
        match = True
        while match and p >= 1:
            if window[p-1] == pattern[p-1]:
                p -= 1
            else:
                match = False

        if match:
            if num_occurences < 10:
                positions.append(t)
            num_occurences += 1
            shift = 1
        else:
            match_length = pattern_length - p
            bcr_shift = bcr_pre.get((window[p-1], p-1), p)
            gsr_shift = gsr_pre.get(match_length, 1)
            shift = max(bcr_shift, gsr_shift)
            if bcr_shift == shift:
                num_bcr_uses += 1

        t += shift
        
        # slide window on search file to the right; fill up window if necessary
        window = window[shift:]
        if len(window) < pattern_length:
            window = fill_window(search_file, window, window_size)

    return num_bcr_uses, num_occurences, positions


def fill_window(search_file, window, window_size):
    """
    Fills the window string with characters from the search_file until the window
    reaches the desired window_size. Removes all newline and carriage return characters.

    Args:
        search_file: File object to read from.
        window: Current window string.
        window_size: Desired length of the window.

    Returns:
        String: The filled window of length up to window_size.
    """
    while len(window) < window_size:
        chunk = search_file.read(window_size - len(window))
        if not chunk:
            break
        window += chunk.replace('\n', '').replace('\r', '')
        
    return window

            
def bcr_preprocessing(seq):
    """
    Preprocesses the pattern for the bad character rule.
    Builds a dictionary mapping (character, position) to the shift value.

    Args:
        seq: String, the pattern to preprocess.

    Returns:
        Dictionary: Mapping (character, position) -> shift value.
    """
    bcr = {}
    
    for i, c in enumerate(seq):        
        for j in range(i+1, len(seq)):
            # shift the search file by j-i to the right 
            bcr[(c,j)] = j-i
    
    return bcr
            

def gsr_preprocessing(seq):
    """
    Preprocesses the pattern for the good suffix rule.
    For each possible suffix, computes the shift needed if a mismatch occurs.

    Args:
        seq: String, the pattern to preprocess.

    Returns:
        Dictionary: Mapping suffix length -> shift value.
    """
    seq_len = len(seq)
    gsr = {}
    
    for suffix_size in range(1, seq_len-1):
        suffix = seq[-suffix_size:]
        bad_character = seq[-suffix_size-1]
        
        # atgca̲g̲g̲aaaaa  ->  a̲g̲g̲aaaaa    (suffix not found)
        # atgta̲g̲g̲            atgtagg
        shift = seq_len - suffix_size+1
        
        for index in range(1, seq_len - suffix_size-1):
           
            # taggca̲g̲g̲aaaaa  ->     taggca̲g̲g̲aaaaa    (suffix found and leading literal different from bad character)
            # taggta̲g̲g̲                  taggtagg
            if seq[-suffix_size - index : -index] == suffix and bad_character != seq[-suffix_size - index - 1]:
                shift = index
                break
            
        gsr[suffix_size] = shift
    
    return gsr
                

def get_all_sequences(file):
    """
    Reads all sequences from a FASTA file and returns them as a list of strings.
    Ignores header lines (starting with '>') and concatenates sequence lines.

    Args:
        file: File object in FASTA format.

    Returns:
        List of strings: Each string is a sequence from the file.
    """
    sequences = []
    current_seq = []
    
    for line in file:
        line = line.strip()
        
        if not line:
            continue
        
        # description line -> new sequence
        if line.startswith(">"):
            if current_seq:
                sequences.append(''.join(current_seq))
                current_seq = []
        else:
            current_seq.append(line)
            
    if current_seq:
        sequences.append(''.join(current_seq))
        
    return sequences


def print_results(result):
    output = " / ".join([str(result[0])] + [str(result[1])] + [str(x) for x in result[2]])
    print(output)
    

def main(args):
    with open(args[0]) as search_file, open(args[1]) as pattern_file:
        pattern_sequences = get_all_sequences(pattern_file)
          
        for seq in pattern_sequences:
            bcr_pre = bcr_preprocessing(seq)
            gsr_pre = gsr_preprocessing(seq)
            
            search_file.seek(0)     # reset file pointer
            search_file.readline()  # skip description line
            print_results(boyer_moore(search_file, seq, bcr_pre, gsr_pre))             
          

if __name__ == '__main__':
    # pass all arguments (not including the invoked .py file)
    #main(sys.argv[1:])
    main(["test1.fasta", "test2.fasta"])