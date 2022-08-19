#include <stdio.h>
#include <vector>
#include <iostream>
#include <unordered_map>
#include <queue>          // std::priority_queue
#include <string>
#include <assert.h>
#include <cstring>
#include <inttypes.h>
#include <stddef.h>
#include <stdint.h>
#include <atomic>
#include <string.h>
#include <algorithm>
#include <numeric>
#include <iostream>
#include <sstream>
#include <random>
#include <utility>
#include <stdlib.h>
#include <assert.h>
#include <errno.h>
#include <stdbool.h>
#include <math.h>
#include <time.h>
#include <omp.h>
#include <mpi.h>
#include "distribute_kmers.h"
#include "bucketing.h"

extern int rank, size;
kmer_t pred_mask=(((kmer_t)1 << (2))-1) << (KMER_LENGTH*2);
kmer_t k_mask=((kmer_t)1 << (KMER_LENGTH*2))-1;

bool comparebykmer(const KmerPair &a, const KmerPair &b)
{
    return extract_kmer(a.bkmer) < extract_kmer(b.bkmer) ||
           ((extract_kmer(a.bkmer) == extract_kmer(b.bkmer)) && (extract_pred(a.bkmer) < extract_pred(b.bkmer)));
}
