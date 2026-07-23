// utils.h

#include <stdlib.h> // For free()

#define MAXSTR 900*2 // Needs MAXSTR from tremer.h, but hardcoding for simplicity in this utility module.

// Safe memory deallocation macro
#define free_null(p)    {if(p)free(p);p=0;}

/**
 * Parses a line of text containing STR values and stores them in a character array.
 * @param line The input string (line from file).
 * @param str The output character array to store integer values.
 * @return The number of STR values successfully read.
 */
int getval(char *line, char *str);

/**
 * Finds the unique values within a list of values.
 * @param v The input list of values.
 * @param uv The output list of unique values.
 * @param n The number of elements in the input list.
 * @param nz If non-zero, only consider non-zero values.
 * @return The number of unique values found.
 */
int unique(char *v, char *uv, int n, int nz);

/**
 * Counts the occurrences of a specific value in a list.
 * @param v The input list of values.
 * @param n The number of elements in the input list.
 * @param u The value to count.
 * @return The count of the value 'u'.
 */
int count_val(char *v, int n, int u);
