#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAX_SEQUENCES 2
#define MATCH 1
#define GAP -1
#define MISMATCH -1
#define max(a,b) (((a) > (b)) ? (a) : (b))

static char **sequences;
static int **mat;
static size_t m;
static size_t n;
static int max_score;


int max4(int a, int b, int c, int d)
{
    return max(max(a, b), max(c, d));
}


void print_matrix()
{
    printf("\n");
    for (int i = 0; i < m + 1; i++)
    {
        for (int j = 0; j < n + 1; j++)
        {
            printf("%4d ", mat[i][j]); // %4d nice format
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


void create_matrix(char *seq1, char *seq2)
{
    m = strlen(seq1);
    n = strlen(seq2);

    mat = malloc((m + 1) * sizeof(int *));

    max_score = 0;

    for (int i = 0; i < m + 1; i++)
    {
        mat[i] = calloc(n + 1, sizeof(int)); // every entry initialised with 0
    }

    for (int i = 1; i < m + 1; i++)
    {
        for (int j = 1; j < n + 1; j++)
        {
            int left = mat[i][j - 1] + GAP;
            int up = mat[i - 1][j] + GAP;
            int diagonal = mat[i - 1][j - 1] + (seq1[i-1] == seq2[j-1] ? MATCH : MISMATCH);

            mat[i][j] = max4(left, up, diagonal, 0);

            if (mat[i][j] > max_score)
            {
                max_score = mat[i][j];
            }
        }
    }
}


void backtrack_core(int i, int j){
    int max_len = max(m, n);
    char *s1 = malloc(max_len + 1);
    char *s2 = malloc(max_len + 1);

    //initialise *-String
    memset(s1, '*', max_len);
    memset(s2, '*', max_len);
    s1[max_len] = '\0';
    s2[max_len] = '\0';

    int score = max_score;
    int index = max(i, j) - 1;  

    while(score > 0){
        // left -> score increases
        if(score == mat[i][j-1] + GAP){
            score -= GAP;
            j--;
            s1[index] = '_';
            s2[index] = sequences[1][j];
        }
        // up -> score increases
        else if(score == mat[i-1][j] + GAP){
            score -= GAP;
            i--;
            s1[index] = sequences[0][i];
            s2[index] = '_';
        }
        // diagonal
        else{
            score -= (sequences[0][i-1] == sequences[1][j-1] ? MATCH : MISMATCH);
            i--;
            j--;
            s1[index] = sequences[0][i];
            s2[index] = sequences[1][j];
        }
        index--;
    }
    
    printf("\n%s\n%s", s1, s2);
    free(s1);
    free(s2);
}


void find_alignment_cores()
{
    for (int i = 1; i < m + 1; i++)
    {
        for (int j = 1; j < n + 1; j++)
        {
            if(mat[i][j] == max_score){
                backtrack_core(i, j);
                break;  //no need to look for other max_scores in this row
            }
        }
    }
}


void free_resources() {
    for (int i = 0; i < MAX_SEQUENCES; i++)
        free(sequences[i]);
    free(sequences);

    for (int i = 0; i < m + 1; i++)
        free(mat[i]);
    free(mat);
}


int main(int argc, char *argv[])
{
    sequences = read_fasta("very_short_test.fasta");
    create_matrix(sequences[0], sequences[1]);

    print_matrix();
    printf("%d", max_score);
    find_alignment_cores();

    free_resources();
    return 0;
}
