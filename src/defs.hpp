#ifndef TYPE_DEFS_HPP
#define TYPE_DEFS_HPP

#include <cstdint>
#include <vector>

// CONSTANTS ---------------------------------------------------------
const unsigned DNA_ALPHABET_SIZE           = 4;
const char DNA_ALPHABET[DNA_ALPHABET_SIZE] = {'A', 'C', 'G', 'T'};
const unsigned ST_ALPHABET_SIZE            = DNA_ALPHABET_SIZE + 1;
const char ST_ALPHABET[ST_ALPHABET_SIZE]   = {'A', 'C', 'G', 'T', '$'};

const unsigned PROT_ALPHABET_SIZE              = 23;
const char PROT_ALPHABET[PROT_ALPHABET_SIZE]   = {'A', 'C', 'D', 'E', 'F', 'G',
                                                  'H', 'I', 'K', 'L', 'M', 'N',
                                                  'P', 'Q', 'R', 'S', 'T', 'V',
                                                  'W', 'Y', 'B', 'Z', 'X'};

// int vectors
typedef std::vector<int32_t> ivec_t;
typedef std::vector<int64_t> ivec64_t;

#endif /* TYPE_DEFS_HPP */
