#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#define MAX_SEQUENCES 2
#define ASCII_SIZE 256
#define GAP '*'
#define max(a,b) (((a) > (b)) ? (a) : (b))

typedef struct {
    int **matrix;
    int score;
    size_t rows;
    size_t cols;
} AlignmentMatrix;


int max3(int a, int b, int c)
{
    return max(max(a, b), c);
}


void print_matrix(int **matrix, char* alphabet)
{
    printf("\n");
    for (int i = 0; alphabet[i]; i++)
    {
        for (int j = 0; alphabet[j]; j++)
        {
            printf("%4d ", matrix[(unsigned char)alphabet[i]][(unsigned char)alphabet[j]]); // %4d nice format
        }
        printf("\n");
    }
}


char **read_fasta(const char *filename)
{
    FILE *file = fopen(filename, "r");

    if (!file)
    {
        perror("File couldn't be opened.");
        return NULL;
    }

    char **sequences = malloc(MAX_SEQUENCES * sizeof(char *));
    size_t *capacities = malloc(MAX_SEQUENCES * sizeof(size_t));
    size_t *lengths = calloc(MAX_SEQUENCES, sizeof(size_t));

    for (int i = 0; i < MAX_SEQUENCES; i++)
    {
        capacities[i] = 128; // start size
        sequences[i] = malloc(capacities[i]);
    }

    char line[1024];
    int seq_index = -1;

    while (fgets(line, sizeof(line), file))
    {
        if (line[0] == '>')
        {
            seq_index++;
            if (seq_index >= MAX_SEQUENCES)
                break;
            continue;
        }

        // strcspn returns index of first occurence
        line[strcspn(line, "\r\n")] = '\0';
        size_t line_len = strlen(line);

        // extend seq len if necessary
        while (lengths[seq_index] + line_len + 1 > capacities[seq_index])
        {
            capacities[seq_index] *= 2;
            sequences[seq_index] = realloc(sequences[seq_index], capacities[seq_index]);
        }

        // append line to seq
        strcpy(sequences[seq_index] + lengths[seq_index], line);
        lengths[seq_index] += line_len;
    }

    fclose(file);
    free(capacities);
    free(lengths);
    return sequences;
}


int** load_scoring_matrix(const char *filename, char **alphabet_out) {
    FILE *f = fopen(filename, "r");
    if (!f) {
        perror("Failed to open file");
        return NULL;
    }

    int** matrix = malloc(ASCII_SIZE * sizeof(int *));

    for (int i = 0; i < ASCII_SIZE; i++){
        matrix[i] = malloc(ASCII_SIZE * sizeof(int));
    }

    char line[1024];
    char col_headers[ASCII_SIZE];
    int num_cols = 0;

    // Skip comment lines starting with '#'
    while (fgets(line, sizeof(line), f)) {
        if (line[0] == '#') {
            continue;
        } else {
            break;
        }
    }


    // Parse column headers (e.g. read all single characters)
    char *token = strtok(line, " \t\r\n");
    while (token != NULL && num_cols < ASCII_SIZE) {
        col_headers[num_cols++] = token[0];
        token = strtok(NULL, " \t\r\n");
    }

    // Extract used alphabet
    *alphabet_out = malloc((num_cols + 1) * sizeof(char));
    memcpy(*alphabet_out, col_headers, num_cols);
    (*alphabet_out)[num_cols] = '\0';


    // Now read the matrix values
    while (fgets(line, sizeof(line), f)) {
        // Skip comment or empty lines
        if (line[0] == '#' || line[0] == '\n' || line[0] == '\r')
            continue;

        char *ptr = line;
        // Trim leading spaces and tabs
        while (*ptr == ' ' || *ptr == '\t') ptr++;

        if (*ptr == '\0') continue; // Empty line, skip

        // First char of the line is the row header
        char row_header = *ptr;
        ptr++;

        int col_idx = 0;
        // Parse integer values for each column
        while (col_idx < num_cols) {
            // Skip whitespace before number
            while (*ptr == ' ' || *ptr == '\t') ptr++;

            if (*ptr == '\0' || *ptr == '\n') break;

            char *endptr;
            int val = (int)strtol(ptr, &endptr, 10);

            if (ptr == endptr) break; // No number found, break

            ptr = endptr;

            // Store the value in the matrix using ASCII indexes
            matrix[(unsigned char)row_header][(unsigned char)col_headers[col_idx]] = val;
            col_idx++;
        }
    }

    fclose(f);
    return matrix;
}


int score(char c1, char c2, int** scoring_matrix){
    return scoring_matrix[(unsigned char) c1][(unsigned char) c2];
}


AlignmentMatrix create_matrix(const char *seq1, const char *seq2, int **scoring_matrix) {
    AlignmentMatrix m;
    m.rows = strlen(seq1) + 1;
    m.cols = strlen(seq2) + 1;

    // extra row and column allow for gap before first character
    m.matrix = malloc(m.rows * sizeof(int *));
    for (size_t i = 0; i < m.rows; i++)
        m.matrix[i] = calloc(m.cols, sizeof(int));


    // initialise first row & first column
    for (size_t i = 1; i < m.cols; i++)
        m.matrix[0][i] = m.matrix[0][i - 1] + score(GAP, seq2[i - 1], scoring_matrix);

    for (size_t i = 1; i < m.rows; i++)
        m.matrix[i][0] = m.matrix[i - 1][0] + score(seq1[i - 1], GAP, scoring_matrix);


    for (size_t i = 1; i < m.rows; i++) {
        for (size_t j = 1; j < m.cols; j++) {
            int left = m.matrix[i][j - 1] + score(GAP, seq2[j - 1], scoring_matrix);
            int up = m.matrix[i - 1][j] + score(seq1[i - 1], GAP, scoring_matrix);
            int diagonal = m.matrix[i - 1][j - 1] + score(seq1[i - 1], seq2[j - 1], scoring_matrix);

            m.matrix[i][j] = max3(left, up, diagonal);
        }
    }
    
    m.score = m.matrix[m.rows-1][m.cols-1];
    return m;
}

void free_alignment_matrix(AlignmentMatrix *m) {
    for (size_t i = 0; i < m->rows; i++)
        free(m->matrix[i]);
    free(m->matrix);
}

void free_scoring_matrix(int **matrix) {
    for (int i = 0; i < ASCII_SIZE; i++)
        free(matrix[i]);
    free(matrix);
}

void free_sequences(char **sequences) {
    for (int i = 0; i < MAX_SEQUENCES; i++)
        free(sequences[i]);
    free(sequences);
}

int main(void) {
    char **sequences = read_fasta("sequence_pair.fasta");
    if (!sequences) return 1;

    char *alphabet = NULL;
    int **scoring_matrix = load_scoring_matrix("blosum62.txt", &alphabet);
    if (!scoring_matrix) {
        free_sequences(sequences);
        return 1;
    }

    AlignmentMatrix m = create_matrix(sequences[0], sequences[1], scoring_matrix);

    printf("%d\n", m.score);


    free_alignment_matrix(&m);
    free_scoring_matrix(scoring_matrix);
    free_sequences(sequences);
    free(alphabet);

    return 0;
}
