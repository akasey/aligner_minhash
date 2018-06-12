//
// Created by Akash Shrestha on 6/11/18.
//


#define TYPE_MATCH      0
#define TYPE_MISMATCH   1
#define TYPE_INSERTION  2
#define TYPE_DELETION   3

#define PENALTY_MATCH      1
#define PENALTY_MISMATCH  -1
#define PENALTY_INSERTION -1
#define PENALTY_DELETION  -1

#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include "ssw.h"

const char kTypeToSymbol[] = {'M', 'X', 'I', 'D'};
#define TypeToSymbol(x) ((x < 4) ? (kTypeToSymbol[x]) : 'X')

void
Verbose(FILE *fp, char *query, uint32_t query_length, char *reference, uint32_t reference_length, std::string &cigar,
        std::vector<std::vector<int32_t> > &dp_matrix, std::vector<std::vector<uint8_t> > &dp_traceback);

int32_t
GenerateCigar(char *query, uint32_t query_length, char *reference, uint32_t reference_length, std::string *ret_cigar,
              uint32_t *ret_alignment_length, std::string *ret_alignment) {
    if (query == NULL || reference == NULL || query_length == 0 || reference_length == 0)
        return -1;

    // Preallocate the memory for the DP matrix and traceback, and initialize all values to 0.
    std::vector<std::vector<int32_t> > dp_matrix((query_length + 1), std::vector<int32_t>((reference_length + 1), 0));
    std::vector<std::vector<uint8_t> > dp_traceback((query_length + 1),
                                                    std::vector<uint8_t>((reference_length + 1), 0));
    // Temporary variables for storing the options for the current step. (Could have been done differently, but increases readability.)
    int32_t up = 0, left = 0, diagonal = 0;

    // Initialize the first row and column according to Needleman-Wunsch.
    for (uint32_t i = 0; i <= query_length; i++)
        dp_matrix[i][0] = i * PENALTY_INSERTION;
    for (uint32_t i = 0; i <= reference_length; i++)
        dp_matrix[0][i] = i * PENALTY_DELETION;

    int32_t max_score_index[] = {0,0};
    int32_t max_score = 0;

    for (uint32_t i = 1; i <= query_length; i++) {
        for (uint32_t j = 1; j <= reference_length; j++) {
            // Calculate the options for the current value of the DP matrix.
            up = dp_matrix[i - 1][j] + PENALTY_INSERTION;
            left = dp_matrix[i][j - 1] + PENALTY_DELETION;
            diagonal =
                    dp_matrix[i - 1][j - 1] + ((query[i - 1] == reference[j - 1]) ? PENALTY_MATCH : PENALTY_MISMATCH);

            // Find the maximum of the three values.
            dp_matrix[i][j] = diagonal;
            dp_traceback[i][j] = ((query[i - 1] == reference[j - 1]) ? 0 : 1);
            if (up > dp_matrix[i][j]) {
                dp_matrix[i][j] = up;
                dp_traceback[i][j] = 2;
            }
            if (left > dp_matrix[i][j]) {
                dp_matrix[i][j] = left;
                dp_traceback[i][j] = 3;
            }
            if (dp_matrix[i][j]<0) {
                dp_matrix[i][j] = 0;
            }

            // If this is the last row of the matrix, we need to look up for the maximum score. This will be our starting point for the
            // CIGAR traceback calculation.
            if (dp_matrix[i][j] > max_score) {
                max_score_index[0] = i, max_score_index[1] = j;
                max_score = dp_matrix[i][j];
            }
        }
    }

    // Initialize the traceback column to the index of the highest score, and the traceback row to the last row of the matrix.
    int32_t current_traceback_row = max_score_index[0], current_traceback_column = max_score_index[1];

    std::stringstream cigar_stream;
    std::string cigar = "";
    std::string alignment = "";
    uint32_t edit_distance = 0;

    // Traceback: Find the CIGAR string.
    uint32_t iterations = 0, num_insertions = 0, type = 0, last_type = 0, num_same_types = 0;
    int32_t last_score = dp_matrix[current_traceback_row][current_traceback_column];
    while (last_score != 0) {//(dp_matrix[current_traceback_row][current_traceback_column] != 0 ) {
        type = dp_traceback[current_traceback_row][current_traceback_column];
        last_score = last_score = dp_matrix[current_traceback_row][current_traceback_column];
        // Count the number of same transitions in the DP matrix (whether if it was to the left, right or up.)
        if (iterations == 0) {
            // If this is the first traceback step, initialize the last type for counting.
            last_type = type;
            num_same_types = 1;
        } else {
            if (type == last_type) {
                // If types of transition are the same, just count them.
                num_same_types += 1;
            } else {
                // If the type has changed, stream the CIGAR count, and reset the type to the current one.
                cigar_stream << TypeToSymbol(last_type) << num_same_types;
                last_type = type;
                num_same_types = 1;
            }
        }

        if (type != TYPE_MATCH)
            edit_distance += 1;

        // Accumulate the characters for the alignment, but only if this was asked for (that's why there is
        // the NULL testing). We insert a special character '+' in this implementation, which represents that
        // there was an insertion on the reference (one character was removed from the query).
        // Character '-' represents the deletion in the reference (or an insertion on the query).
        if (ret_alignment != NULL)
            alignment += ((type == TYPE_DELETION) ? '-' :
                          ((type == TYPE_MATCH) ? query[current_traceback_row - 1] :
                           ((type == TYPE_INSERTION) ? '+' :
                            'x')));

        if (type == TYPE_INSERTION)
            num_insertions += 1;

        // Row changes in all cases except when a deletion occurs (type == 3).
        current_traceback_row -= (type != TYPE_DELETION) ? 1 : 0;
        // Column changes in all cases except when an insertion occurs (type == 2).
        current_traceback_column -= (type != TYPE_INSERTION) ? 1 : 0;

        // Count the traceback length.
        iterations += 1;
    }

    // Stream the remaining count to the CIGAR stream.
    cigar_stream << TypeToSymbol(last_type) << num_same_types;
    cigar = cigar_stream.str();
    std::reverse(cigar.begin(), cigar.end());
    std::reverse(alignment.begin(), alignment.end());

    if (ret_cigar != NULL)
        *ret_cigar = cigar;
    else
        Verbose(stdout, query, query_length, reference, reference_length, cigar, dp_matrix, dp_traceback);

    // Return the final alignment length. We need to subtract the number of insertions on the reference because the final
    // alignment should be shorter than the traceback for this amount.
    if (ret_alignment_length != NULL)
        *ret_alignment_length = (iterations - num_insertions);

    // Return the final alignment if required. The alignment.size() value is larger than ret_alignment_length, because
    // query deletions have not been completely removed, but changed with a special character.
    if (ret_alignment != NULL)
        *ret_alignment = alignment;

    return edit_distance;
}

/*int main() {
    std::string reference = "ACTGCTGCCTGCAAAAAAAAAAA";
    std::string query = "TGCCTGCAA";

    std::string cigar;
    std::string alignment = "";
    uint32_t alignment_length = 0;
    int32_t edit_distance = GenerateCigar((char *) query.c_str(), query.size(), (char *) reference.c_str(), reference.size(), &cigar, &alignment_length, &alignment);

    printf ("Reference sequence:\t%s\n", reference.c_str());
    printf ("Query sequence:\t\t%s\n", query.c_str());
    printf ("Returned CIGAR: %s\n", cigar.c_str());
    printf ("Alignment: %s\n", alignment.c_str());
    printf ("Alignment length: %d\n", alignment_length);
    printf ("Edit distance: %d\n", edit_distance);

    return 0;
}*/

void
Verbose(FILE *fp, char *query, uint32_t query_length, char *reference, uint32_t reference_length, std::string &cigar,
        std::vector<std::vector<int32_t> > &dp_matrix, std::vector<std::vector<uint8_t> > &dp_traceback) {
    fprintf(fp, "%3c", ' ');
    for (uint32_t i = 0; i <= reference_length; i++) {
        if (i > 0)
            fprintf(fp, "%3c", reference[i - 1]);
        else
            fprintf(fp, "%3c", '-');
    }
    fprintf(fp, "\n");

    for (uint32_t i = 0; i <= query_length; i++) {
        if (i > 0)
            fprintf(fp, "%3c", query[i - 1]);
        else
            fprintf(fp, "%3c", '-');

        for (uint32_t j = 0; j <= reference_length; j++) {
            fprintf(fp, "%3d", dp_matrix[i][j]);
        }
        fprintf(fp, "\n");
    }

    fprintf(fp, "\n");
    fprintf(fp, "Traceback:\n");

    fprintf(fp, "%3c", ' ');
    for (uint32_t i = 0; i <= reference_length; i++) {
        if (i > 0)
            fprintf(fp, "%3c", reference[i - 1]);
        else
            fprintf(fp, "%3c", '-');
    }
    fprintf(fp, "\n");

    for (uint32_t i = 0; i <= query_length; i++) {
        if (i > 0)
            fprintf(fp, "%3c", query[i - 1]);
        else
            fprintf(fp, "%3c", '-');

        for (uint32_t j = 0; j <= reference_length; j++) {
            fprintf(fp, "%3d", dp_traceback[i][j]);
        }
        fprintf(fp, "\n");
    }

    fprintf(fp, "\n");

    fprintf(fp, "CIGAR: %s\n", cigar.c_str());
    fprintf(fp, "\n");
}



//	Print the BLAST like output.
static void ssw_write (const s_align* a,
                       const char* ref_seq,
                       const char* read_seq,
                       const int8_t* table) {

    fprintf(stdout, "optimal_alignment_score: %d\tsub-optimal_alignment_score: %d\t", a->score1, a->score2);
    if (a->ref_begin1 + 1) fprintf(stdout, "target_begin: %d\t", a->ref_begin1 + 1);
    fprintf(stdout, "target_end: %d\t", a->ref_end1 + 1);
    if (a->read_begin1 + 1) fprintf(stdout, "query_begin: %d\t", a->read_begin1 + 1);
    fprintf(stdout, "query_end: %d\n\n", a->read_end1 + 1);
    if (a->cigar) {
        int32_t c = 0, left = 0, e = 0, qb = a->ref_begin1, pb = a->read_begin1;
        uint32_t i;
        while (e < a->cigarLen || left > 0) {
            int32_t count = 0;
            int32_t q = qb;
            int32_t p = pb;
            fprintf(stdout, "Target: %8d    ", q + 1);
            for (c = e; c < a->cigarLen; ++c) {
                char letter = cigar_int_to_op(a->cigar[c]);
                uint32_t length = cigar_int_to_len(a->cigar[c]);
                uint32_t l = (count == 0 && left > 0) ? left: length;
                for (i = 0; i < l; ++i) {
                    if (letter == 'I') fprintf(stdout, "-");
                    else {
                        fprintf(stdout, "%c", *(ref_seq + q));
                        ++ q;
                    }
                    ++ count;
                    if (count == 60) goto step2;
                }
            }
            step2:
            fprintf(stdout, "    %d\n                    ", q);
            q = qb;
            count = 0;
            for (c = e; c < a->cigarLen; ++c) {
                char letter = cigar_int_to_op(a->cigar[c]);
                uint32_t length = cigar_int_to_len(a->cigar[c]);
                uint32_t l = (count == 0 && left > 0) ? left: length;
                for (i = 0; i < l; ++i){
                    if (letter == 'M') {
                        if (table[(int)*(ref_seq + q)] == table[(int)*(read_seq + p)])fprintf(stdout, "|");
                        else fprintf(stdout, "*");
                        ++q;
                        ++p;
                    } else {
                        fprintf(stdout, "*");
                        if (letter == 'I') ++p;
                        else ++q;
                    }
                    ++ count;
                    if (count == 60) {
                        qb = q;
                        goto step3;
                    }
                }
            }
            step3:
            p = pb;
            fprintf(stdout, "\nQuery:  %8d    ", p + 1);
            count = 0;
            for (c = e; c < a->cigarLen; ++c) {
                char letter = cigar_int_to_op(a->cigar[c]);
                uint32_t length = cigar_int_to_len(a->cigar[c]);
                uint32_t l = (count == 0 && left > 0) ? left: length;
                for (i = 0; i < l; ++i) {
                    if (letter == 'D') fprintf(stdout, "-");
                    else {
                        fprintf(stdout, "%c", *(read_seq + p));
                        ++p;
                    }
                    ++ count;
                    if (count == 60) {
                        pb = p;
                        left = l - i - 1;
                        e = (left == 0) ? (c + 1) : c;
                        goto end;
                    }
                }
            }
            e = c;
            left = 0;
            end:
            fprintf(stdout, "    %d\n\n", p);
        }
    }
}

//	Align a pair of genome sequences.
int main (int argc, char * const argv[]) {
    int32_t l, m, k, match = 2, mismatch = 2, gap_open = 3, gap_extension = 1;	// default parameters for genome sequence alignment
    // reference sequence
    static const char ref_seq[40] = {'C', 'A', 'G', 'C', 'C', 'T', 'T', 'T', 'C', 'T', 'G', 'A', 'C', 'C', 'C', 'G', 'G', 'A', 'A', 'A', 'T',
                                     'C', 'A', 'A', 'A', 'A', 'T', 'A', 'G', 'G', 'C', 'A', 'C', 'A', 'A', 'C', 'A', 'A', 'A', '\0'};
    static const char read_seq[16] = {'C', 'T', 'G', 'A', 'G', 'C', 'C', 'G', 'G', 'T', 'A', 'A', 'A', 'T', 'C', '\0'};	// read sequence
    s_profile* profile;
    int8_t* num = (int8_t*)malloc(16);	// the read sequence represented in numbers
    int8_t* ref_num = (int8_t*)malloc(64);	// the read sequence represented in numbers
    s_align* result;

    /* This table is used to transform nucleotide letters into numbers. */
    static const int8_t nt_table[128] = {
            4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
            4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
            4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
            4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
            4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
            4, 4, 4, 4,  3, 0, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
            4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
            4, 4, 4, 4,  3, 0, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
    };

    // initialize scoring matrix for genome sequences
    //  A  C  G  T	N (or other ambiguous code)
    //  2 -2 -2 -2 	0	A
    // -2  2 -2 -2 	0	C
    // -2 -2  2 -2 	0	G
    // -2 -2 -2  2 	0	T
    //	0  0  0  0  0	N (or other ambiguous code)
    int8_t* mat = (int8_t*)calloc(25, sizeof(int8_t));
    for (l = k = 0; l < 4; ++l) {
        for (m = 0; m < 4; ++m) mat[k++] = l == m ? match : - mismatch;	/* weight_match : -weight_mismatch */
        mat[k++] = 0; // ambiguous base: no penalty
    }
    for (m = 0; m < 5; ++m) mat[k++] = 0;

    for (m = 0; m < 15; ++m) num[m] = nt_table[(int)read_seq[m]];
    profile = ssw_init(num, 15, mat, 5, 2);
    for (m = 0; m < 39; ++m) ref_num[m] = nt_table[(int)ref_seq[m]];

    // Only the 8 bit of the flag is setted. ssw_align will always return the best alignment beginning position and cigar.
    result = ssw_align (profile, ref_num, 39, gap_open, gap_extension, 1, 0, 0, 15);
    ssw_write(result, ref_seq, read_seq, nt_table);

    align_destroy(result);
    init_destroy(profile);
    free(mat);
    free(ref_num);
    free(num);
    return(0);
}