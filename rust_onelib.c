// Rust-friendly wrappers for ONElib functions
// These provide clean accessors without macro complexity

#include <stdint.h>
#include "ONElib.h"
#include "alncode.h"

// Field accessors (reading)
int64_t one_int(OneFile *of, int index) {
    return oneInt(of, index);
}

double one_real(OneFile *of, int index) {
    return oneReal(of, index);
}

char one_char(OneFile *of, int index) {
    return oneChar(of, index);
}

// Get current line type
char one_line_type(OneFile *of) {
    return of->lineType;
}

// Get line count for this line type
int64_t one_line_count(OneFile *of) {
    return of->line;
}

// Field setters (writing)
void one_int_set(OneFile *of, int index, int64_t value) {
    oneInt(of, index) = value;
}

void one_real_set(OneFile *of, int index, double value) {
    oneReal(of, index) = value;
}

void one_char_set(OneFile *of, int index, char value) {
    oneChar(of, index) = value;
}
