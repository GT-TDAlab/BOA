#ifndef BUCKETING_H
#define BUCKETING_H

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

extern kmer_t pred_mask;
extern kmer_t k_mask;

inline kmer_t kmerplusone_shift(kmer_t kmer_in,
                                ElType el)
{   
    return (kmer_t)((kmer_in<<2) | (kmer_t)el) & (kmer_t)KPLUS1_MASK;
}

inline kmer_t kmer_shift_rc(kmer_t kmer_in,
                                ElType el) {

      return (kmer_t)((kmer_in>>2) | ((kmer_t)el << (2*(KMER_LENGTH-1)))) & (kmer_t)KMER_MASK;
}

inline kmer_t extract_pred(kmer_t kmer_in) {

    return ((pred_mask & kmer_in) >> (KMER_LENGTH*2));
}

inline kmer_t extract_kmer(kmer_t kmer_in) {

    return (k_mask & kmer_in);
}

struct KmerPair
{ 
  kmer_t bkmer;
  read_t read_id;

  bool operator==(const KmerPair& a) const
  {
    return bkmer == a.bkmer;
  }

  auto operator<(const KmerPair& r) const
  {

    return std::make_tuple(extract_kmer(bkmer), extract_pred(bkmer), read_id) < std::make_tuple(extract_kmer(r.bkmer), extract_pred(r.bkmer), r.read_id);

   }

  KmerPair()=default;

};

struct RevKmerPair
{ 
  kmer_t fkmer;
  lmer_t min_lmer;
  bool is_blank=false;

  RevKmerPair()=default;

};


struct KmerBucketPair
{ 
  kmer_t seq;
  std::vector<read_t> ReadIds;
  ElType pred_bp;

  bool operator<(const KmerBucketPair& r) const
  {

     return std::make_pair(extract_kmer(seq), pred_bp) < std::make_pair(extract_kmer(r.seq), r.pred_bp);

  }

  auto operator==(const KmerBucketPair& b) const
  {

     return (extract_kmer(seq) == extract_kmer(b.seq)) && (pred_bp == b.pred_bp);

  }

  KmerBucketPair()=default;

};


bool comparebykmer(const KmerPair &a, const KmerPair &b);

void print_kmer_buckets_pair (int klen,
                              std::vector < std::vector<KmerBucket> > &kmer_send_buf);

void print_kmer_buckets_rank (int klen,
                              std::vector <KmerPair> &kmer_send_buf);

void create_readID_list(std::vector<KmerBucketPair> &vec);

size_t Sliding_window (const char *ptr, size_t length, int *n_kmers, 
                     std::vector<std::vector<KmerPair>> &partial_kmer_counts,
                     std::vector<std::vector<KmerPair>> &partial_blank_kmer_counts,
                     uint64_t reads_disp, size_t read_start_idx);
                  
void transfer_kmers (std::vector<int> &scounts_kmer, 
                     std::vector < std::vector<KmerBucket> > &kmer_send_buf,
                     std::vector<int> &send_num_blanks);

void process_remaining_kmers(
                     std::vector<std::vector<KmerPair>> &partial_kmer_counts,
                     std::vector<std::vector<KmerPair>> &partial_blank_kmer_counts);
#endif
