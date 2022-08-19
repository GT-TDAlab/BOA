#ifndef DISTRIBUTE_KMERS_H
#define DISTRIBUTE_KMERS_H

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
#include <atomic>
#include <string.h>
#include <algorithm>
#include <numeric>
#include <stdlib.h>
#include <assert.h>
#include <errno.h>
#include <stdbool.h>
#include <math.h>
#include <time.h>
#include <omp.h>
#include <mpi.h>
#include <cstdint>
#include "zoltan_partition.h"
#include "defs.h"

#define KMER_LENGTH                    (WINDW_SIZE+1)
#define KPLUSONE_LENGTH                (KMER_LENGTH+1)
#define LMER_SIZE                      (pow(4, LMER_LENGTH))
#define MN_LENGTH                      (KMER_LENGTH-1)

#define MOD 2147483647
#define HL 31


#ifdef EXTEND_KMER
typedef __uint128_t kmer_t;
#else
typedef uint64_t kmer_t;
#endif

enum PartitionerType {Zoltan=1, ParMetis=2};
using BasePair = uint8_t;
typedef BasePair ElType;
typedef uint64_t lmer_t;

static const size_t bsize_w = 3;
static const size_t nels_per_value_w = (sizeof(read_t)*8)/bsize_w;

#ifdef EXTEND_KMER
#define KMER_MASK            ((__uint128_t)(~0L)) >> ((sizeof(kmer_t)*8) - (2*KMER_LENGTH))
#define BASE_MASK            ((__uint128_t)(~0L)) >> ((sizeof(kmer_t)*8) - (2*64))
#define LMER_MASK            ((__uint128_t)(~0L)) >> ((sizeof(lmer_t)*8) - (2*LMER_LENGTH))
#define MN_MASK              ((__uint128_t)(~0L)) >> ((sizeof(kmer_t)*8) - (2*MN_LENGTH))
#define SUCC_MASK            ((__uint128_t)(~0L)) >> ((sizeof(kmer_t)*8) - (2*(KMER_LENGTH-1)))
#else
#define KMER_MASK           ((~0UL)) >> ((sizeof(kmer_t)*8) - (2*KMER_LENGTH))
#define BASE_MASK           ((~0UL)) >> ((sizeof(kmer_t)*8) - (bsize_w*nels_per_value_w))
#define LMER_MASK           ((~0UL)) >> ((sizeof(lmer_t)*8) - (2*LMER_LENGTH))
#define MN_MASK             ((~0UL)) >> ((sizeof(kmer_t)*8) - (2*MN_LENGTH))
#define SUCC_MASK           ((~0UL)) >> ((sizeof(kmer_t)*8) - (2*(KMER_LENGTH-1)))
#define KPLUS1_MASK         ((~0UL)) >> ((sizeof(kmer_t)*8) - (2*KPLUSONE_LENGTH))
#endif


enum KTerminal {N, Y};
enum KdirType {P, S};
extern int rank, size;


template<typename T>
inline bool
operator == (const std::vector<T>& v1,
             const std::vector<T>& v2) {
  if(v1.size() != v2.size()) {
    return false;
  }
  return std::equal(v1.begin(), v1.end(), v2.begin());
}


inline char el_to_char(unsigned x) {
      static char symbols[] = {'A', 'C', 'T', 'G', '\0', '\0', 'L', 'N'};

        return symbols[(ElType)x];
}

inline kmer_t mnmer_shift(kmer_t kmer_in, 
                        ElType el) {

  
  return (kmer_t)((kmer_in<<2) | (kmer_t)el) & (kmer_t)MN_MASK;
}


class Comp_rev{
    const std::vector<std::pair<int, int>> & _v;
  public:
    Comp_rev(const std::vector<std::pair<int, int>> & v) : _v(v) {}
    bool operator()(size_t i, size_t j){
         return ((_v[i].second > _v[j].second) ||
                 ((_v[i].second == _v[j].second) &&
                  ( _v[i].first > _v[j].first))
                );

   }
};

struct KmerBucket
{ 
  kmer_t seq;
  std::vector<read_t> ReadIds;

  KmerBucket()=default;

};


typedef struct begin_kmer_id
{
    kmer_t node;
    int terminal_prefix_id;

} BeginMN;



typedef struct __attribute__ ((__packed__)) kmer_pairs
{ 
  kmer_t seq;
  int k_count;

} KmerPairs;
static_assert(sizeof(KmerPairs) == (sizeof(kmer_t)+sizeof(int)), "struct KmerPairs shouldn't be padded");

//data structure for storing Reads
typedef struct Rd_Sequence
{
	// Read data
	size_t num_of_reads;
    size_t *read_data_ptr;
	char *read_data;
    size_t read_data_size;

    Rd_Sequence() : num_of_reads(0), read_data_ptr(nullptr), read_data(nullptr), read_data_size(0) { }

} input_read_data;


class Comp{
    const std::vector<kmer_t> & _v;
  public:
    Comp(const std::vector<kmer_t> & v) : _v(v) {}
    bool operator()(size_t i, size_t j){
         return _v[i] < _v[j];
   }
};

inline ElType kmerel(kmer_t kmer, unsigned pos) {
  assert(pos < KMER_LENGTH);

  return ElType(((kmer) >> ((KMER_LENGTH-(pos)-1)*2)) & (0x3));
}

inline ElType lmerel(lmer_t kmer, unsigned pos) {
  assert(pos < LMER_LENGTH);

  return ElType(((kmer) >> ((LMER_LENGTH-(pos)-1)*2)) & (0x3));
}

inline ElType kmerel_mn(kmer_t kmer, unsigned pos) {
      assert(pos < MN_LENGTH);

        return ElType(((kmer) >> ((MN_LENGTH-(pos)-1)*2)) & (0x3));
}


inline lmer_t kmer_to_lmer(kmer_t kmer_in, unsigned pos, lmer_t kmer_out) {
  assert(pos < KMER_LENGTH);

  return (lmer_t)((kmer_out<<2) | (lmer_t)(ElType(((kmer_in) >> ((KMER_LENGTH-(pos)-1)*2)) & (0x3)))) & (lmer_t)LMER_MASK;
}

// A/a: 0
// T/t: 2
// G/g: 3
// C/c: 1
// L/l: 6
// N/n: 7
inline ElType char_to_el(char ch) {
              return (ElType)((((ElType)ch)>>1) & 0x7);
}

inline ElType rev_comp(char ch) {
       char rch;

       switch (ch) {
       case 'A':
             rch='T';
             break;
       case 'a':
             rch='t';
             break;
       case 'C':
             rch='G';
             break;
       case 'c':
             rch='g';
             break;
       case 'G':
             rch='C';
             break;
       case 'g':
             rch='c';
             break;
       case 'T':
             rch='A';
             break;
       case 't':
             rch='a';
             break;
       default:
             printf("Error: bp:%c not valid \n", ch);
             assert(0 && "bp must be either A/C/G/T, Abort!!");
             break;
      }

      return char_to_el(rch);
}

inline kmer_t kmer_cons(kmer_t kmer_in, 
                        unsigned pos,
                        ElType el) {
  assert(pos < KMER_LENGTH);
 
  return kmer_in | ((kmer_t)el << ((KMER_LENGTH-pos-1)*2));
}

inline kmer_t kmer_shift(kmer_t kmer_in, 
                        ElType el) {

  return (kmer_t)((kmer_in<<2) | (kmer_t)el) & (kmer_t)KMER_MASK;
}

inline lmer_t lmer_shift(lmer_t kmer_in, 
                        ElType el) {

  return (lmer_t)((kmer_in<<2) | (lmer_t)el) & (lmer_t)LMER_MASK;
}

inline kmer_t mn_shift(kmer_t kmer_in)
{
    return (kmer_t)((kmer_in>>2)) & (kmer_t)MN_MASK;
}

inline kmer_t tokmer(const char *kmer_str,
                     int kmer_len) {
  int i;
  assert(kmer_len == KMER_LENGTH);
  assert(kmer_str != NULL);
  assert(kmer_str[kmer_len] == '\0');

  kmer_t km = 0;
  for(i=0; i<kmer_len; i++) {
    ElType el = char_to_el(kmer_str[i]);
    km = kmer_cons(km, i, el);
  }  
  return km;
}

// This function only is used for generating the read_t, so we should the wide version
template<typename T>
void set_bp(BasePair val, size_t pos, std::vector<T>& vec_) {
    assert(val < (1<<bsize_w));
    size_t idx = pos / nels_per_value_w;

    vec_[idx] = ((vec_[idx]<<bsize_w) | val) & static_cast<T>(BASE_MASK);
};


// This function only is used for generating the read_t, so we should the wide version
template<typename T>
void push_back_bp_w(BasePair val, int *size_, std::vector<T>& vec_)
{
    assert(val < (1<<bsize_w));
    *size_ += 1;
    if((*size_%nels_per_value_w) == 1) {
      vec_.push_back(T{0});
    }
    set_bp(val, *size_-1, vec_);
}

template<typename T>
BasePair get_bp(size_t pos, std::vector<T>& vec_, size_t off, int read_length) {
    size_t idx = pos / nels_per_value_w;
    size_t idx_p = read_length / nels_per_value_w;
    size_t off_size = idx < idx_p ? nels_per_value_w : read_length % nels_per_value_w;
    size_t off_pos = pos % nels_per_value_w;

    return ((vec_[idx+off]) >> ((off_size-(off_pos)-1)*bsize_w)) & ((1<<bsize_w)-1);

}


template <typename T> 
inline long uhash31( uint64_t a, uint64_t b, T x)
{

  T result;
  long lresult;  

  result=(a * x) + b;
  result = ((result >> HL) + result) & MOD;
  lresult=(long) result; 
  
  return(lresult);
}

inline long hash31(long long a, long long b, long long x)
{

  long long result;
  long lresult;  

  result=(a * x) + b;
  result = ((result >> HL) + result) & MOD;
  lresult=(long) result; 
  
  return(lresult);
}
 
int hash_fcn (const char* word, unsigned M);
void parse_alphabet_file (FILE * fp);
void free_reads (int num_reads);
bool is_sigma (char c);
void change_to_num (char *pred_kmer, int *num);
void change_to_char (char *pred_kmer, int num);
int find_max (int array[]);
int find_min (int array[], int *min);
int compare (const void * a, const void * b);
int FindIndex( const int a[], int size, int value);
char* DivideReads(MPI_File *in, const int rank, const int size, 
        const int overlap, uint64_t *nlines, size_t *data_size);
input_read_data perform_input_reading (const int rank, const int size,
                                       std::string &fileName,
                                       uint64_t *reads_disp, uint64_t *nlines);
void perform_kmer_counting(input_read_data &rdata_fwd, input_read_data &rdata_rev, uint64_t reads_disp, uint64_t local_nlines);
void print_kmers (int klen, int min_t);
void SortAndAggregate(std::vector<KmerPairs>& arr);
void sort_recv_buffer(std::vector<KmerPairs>& kmer_recv_buf, 
                      std::vector<int>& rcounts_kmer, 
                      std::vector<int>& rdisp_kmer);
void Sliding_window_l (const char *ptr, size_t length);
void Sliding_window_l_variableReadLength (input_read_data & reads);

void free_kmer_count_buffers();
int convert_hash_fcn (const char* word);
unsigned int hash_str(const char *str);

int retrieve_proc_id (lmer_t min_lmer);
int retrieve_proc_id (kmer_t min_lmer);
HGRAPH_DATA populate_hgraph_data(uint64_t reads_disp, uint64_t num_reads_per_p);
GRAPH_DATA populate_graph_data(uint64_t reads_disp, uint64_t num_reads_per_p);
void print_hgraph_bucket_data (HGRAPH_DATA &hg);
void print_hgraph_read_data (HGRAPH_DATA &hg);
void print_graph_read_data (GRAPH_DATA &g);

//////////////////////////////////////////////////////////////
#endif
