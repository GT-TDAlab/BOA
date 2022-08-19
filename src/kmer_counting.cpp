#include <stdio.h>
#include <mpi.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <inttypes.h>
#include <vector>
#include <set>
#include <algorithm>
#include <unordered_map>
#include <climits>
#include <parallel/algorithm>
#include <numeric>
#include <omp.h>
#include "distribute_kmers.h"
#include "timers.h"
#include "bucketing.h"
#include "serialize.h"
#include "zoltan_partition.h"
#include <filesystem>
#include <chrono>
#include <thread>

extern long int MAX_KMER_COUNT;
extern int rank, size;
extern int coverage;
extern int num_threads;
extern int num_batch_transfers;
extern bool ispaired;
extern int read_length;
extern std::string outputDirName;
uint64_t num_recalculate_lmer, global_num_recalculate_lmer;
extern uint64_t global_nlines;
extern uint64_t global_nvertices;
int max_reads_per_bucket=201;
extern PartitionerType ptype;
extern int read_size;

extern std::vector<lmer_t> lmer_frequency;
extern std::vector<lmer_t> global_lmer_frequency;

std::vector<KmerBucketPair> kmer_proc_buf;
extern std::vector<read_t> col_ind_readid;
extern std::vector<kmer_t> bucket_list;
extern std::vector<int> row_bucket_ptr;

extern std::vector<size_t> vtxDataPtrs;
extern std::vector<read_t> reads_data;

HGRAPH_DATA populate_hgraph_data(uint64_t reads_disp, uint64_t num_reads_per_p)
{
     return HGRAPH_DATA(num_reads_per_p, reads_disp, vtxDataPtrs, reads_data,
             bucket_list, row_bucket_ptr, col_ind_readid);
}

GRAPH_DATA populate_graph_data(uint64_t reads_disp, uint64_t num_reads_per_p)
{
     return GRAPH_DATA(num_reads_per_p, reads_disp, vtxDataPtrs, reads_data,
             row_bucket_ptr, col_ind_readid);
}

void sort_recv_buffer(std::vector<KmerPairs>& kmer_recv_buf, 
                      std::vector<int>& rcounts_kmer, 
                      std::vector<int>& rdisp_kmer)
{
     size_t rsize=0;
     std::vector<int> offsets(rdisp_kmer);
     offsets.push_back(rcounts_kmer[size-1]+rdisp_kmer[size-1]);
     for (int t=0; t<rcounts_kmer.size(); t++) rsize += rcounts_kmer[t];
     assert(kmer_recv_buf.size() == rsize);
     assert(offsets.size() == (size+1));

     while(offsets.size()>2) {
            assert(offsets.front() == 0);
            std::vector<int> new_offsets;
            int x = 0;
            while(x+2 < offsets.size()) {
                    // mergesort (offsets[x],offsets[x+1]) and (offsets[x+1],offsets[x+2])
                    std::inplace_merge(kmer_recv_buf.begin()+offsets[x]
                                 ,kmer_recv_buf.begin()+offsets[x+1]
                                 ,kmer_recv_buf.begin()+offsets[x+2] // this *might* be at the end
                                 ,[](const auto& i, const auto& j) {return i.seq < j.seq;} 
                                 );
                    // now they are sorted, we just put offsets[x] and offsets[x+2] into the new offsets.
                    // offsets[x+1] is not relevant any more
                    new_offsets.push_back(offsets[x]);
                    new_offsets.push_back(offsets[x+2]);
                    x += 2;
            }
            // if the number of offsets was odd, there might be a dangling offset
            // which we must remember to include in the new_offsets
            if(x+2==offsets.size()) {
                    new_offsets.push_back(offsets[x+1]);
            }
            offsets.swap(new_offsets);

    }
    offsets.clear();
    offsets.shrink_to_fit();
}  

void SortAndAggregate(std::vector<KmerPairs>& arr)
{

    std::vector<KmerPairs>::iterator low,up, it;
    std::vector<KmerPairs> new_arr;

    for( it = arr.begin(); it != arr.end(); )
          {
              kmer_t key = (*it).seq;
              low=std::lower_bound (arr.begin(), arr.end(), key,
                                   [] (const KmerPairs& lhs, kmer_t rhs) {
                             return (lhs.seq < rhs);
                             });
 
              up= std::upper_bound (arr.begin(), arr.end(), key,
                                   [] (kmer_t rhs, const KmerPairs& lhs) {
                             return (rhs < lhs.seq);
                             });

              int sum=0;
              for (auto itr=low; itr!= up; itr++) {
                        sum += (*itr).k_count;
              }
              new_arr.push_back(KmerPairs{key, sum});

              it = up;
          }

     arr=new_arr;
     new_arr.clear();
     new_arr.shrink_to_fit();
}

int SortAndAggregatePair(std::vector<KmerPair>& arr,
                          std::vector<KmerBucket>& new_arr)
{

    std::vector<KmerPair>::iterator low,up, it;
    int sum=0;
    for( it = arr.begin(); it != arr.end(); )
          {
              kmer_t kmerplusone = (*it).bkmer;
              std::pair<kmer_t, kmer_t> key = std::make_pair(extract_kmer(kmerplusone), extract_pred(kmerplusone));
              low=std::lower_bound (arr.begin(), arr.end(), key,
                                   [] (const KmerPair& lhs, std::pair<kmer_t, kmer_t> rhs) {
                             return (std::make_pair(extract_kmer(lhs.bkmer), extract_pred(lhs.bkmer)) < rhs);
                             });
 
              up= std::upper_bound (arr.begin(), arr.end(), key,
                                   [] (std::pair<kmer_t, kmer_t> rhs, const KmerPair& lhs) {
                             return (rhs < std::make_pair(extract_kmer(lhs.bkmer), extract_pred(lhs.bkmer)));
                             });

               std::vector<read_t> read_list;
              for (auto itr=low; itr!= up; itr++) {
                  read_list.push_back((*itr).read_id);
              }
              
              // performing the sort and deletion of temp duplicates for now
              std::sort(read_list.begin(), read_list.end());
              read_list.erase(std::unique(read_list.begin(), read_list.end()), read_list.end());
              new_arr.push_back(KmerBucket{kmerplusone, read_list});
              sum++;

              it = up;
          }

     return sum;
}

template <class It>
auto compute_duplicates(It begin, It end) -> std::pair<std::vector<read_t>, It>
{
  auto it = begin + 1;
  std::vector<read_t> read_ids;
  read_ids.push_back(begin->read_id);
  for (; it != end && *begin == *it; ++it) {
       read_ids.push_back(it->read_id);
  }

  std::sort(read_ids.begin(), read_ids.end());
  read_ids.erase(std::unique(read_ids.begin(), read_ids.end()), read_ids.end());
  
  return std::make_pair(read_ids, it);
}

template <class It>
auto aggregate_readIds(It start, It end) -> std::pair<std::vector<read_t>, It>
{ 
  auto it = start + 1;
  std::vector<read_t> read_id;
  read_id.insert(std::end(read_id), std::begin(start->ReadIds), std::end(start->ReadIds));
  for (; it != end && *start == *it; ++it) {
       read_id.insert(std::end(read_id), std::begin(it->ReadIds), std::end(it->ReadIds));
  }
  
  std::sort(read_id.begin(), read_id.end());
  read_id.erase(std::unique(read_id.begin(), read_id.end()), read_id.end());
  return std::make_pair(read_id, it);
}

void create_readID_list(std::vector<KmerBucketPair> &vec)
{
  std::vector<KmerBucketPair> vec_combine;

  for (std::vector<KmerBucketPair>::iterator i = std::begin(vec); i != std::end(vec); ++i) {
       decltype(i) j;
       std::vector<read_t> readId_list;
       std::tie(readId_list, j) = aggregate_readIds<std::vector<KmerBucketPair>::iterator>(i, std::end(vec));

       vec_combine.push_back(KmerBucketPair{i->seq, readId_list, i->pred_bp});
       vec.erase(i + 1, j);
  }
  vec=vec_combine;
  vec_combine.clear();
  vec_combine.shrink_to_fit();
}

void SortAndAggregatePairBucket(std::vector<KmerBucketPair>& arr)
{

    std::vector<KmerBucketPair>::iterator low,up, it;
    std::vector<KmerBucketPair> new_arr;
    
    for( it = arr.begin(); it != arr.end(); )
          {
              kmer_t kmerplusone = (*it).seq;
              ElType prec_bp = (*it).pred_bp;
              std::pair<kmer_t, ElType> key = std::make_pair(extract_kmer(kmerplusone), prec_bp);
              low=std::lower_bound (arr.begin(), arr.end(), key,
                                   [] (const KmerBucketPair& lhs, std::pair<kmer_t, ElType> rhs) {
                             return (std::make_pair(extract_kmer(lhs.seq), lhs.pred_bp) < rhs);
                             });
 
              up= std::upper_bound (arr.begin(), arr.end(), key,
                                   [] (std::pair<kmer_t, ElType> rhs, const KmerBucketPair& lhs) {
                             return (rhs < std::make_pair(extract_kmer(lhs.seq), lhs.pred_bp));
                             });

              std::vector<read_t> read_list;
              for (auto itr=low; itr!= up; itr++) {
                  read_list.insert(std::end(read_list), std::begin((*itr).ReadIds), std::end((*itr).ReadIds));
              }
              
              // performing the sort and deletion of temp duplicates for now
              std::sort(read_list.begin(), read_list.end());
              read_list.erase(std::unique(read_list.begin(), read_list.end()), read_list.end());
              new_arr.push_back(KmerBucketPair{kmerplusone, read_list, prec_bp});

              it = up;
          }

     
     arr=new_arr;
     new_arr.clear();
     new_arr.shrink_to_fit();
}

void GenerateSortedReadListPerBucket(std::vector<KmerBucketPair>& arr,
                                     std::vector<KmerBucket>& new_arr)
{

    std::vector<KmerBucketPair>::iterator low,up, it;

    for( it = arr.begin(); it != arr.end(); )
          {
              kmer_t kmerplusone = (*it).seq;
              kmer_t key = extract_kmer(kmerplusone);
              low=std::lower_bound (arr.begin(), arr.end(), key,
                                   [] (const KmerBucketPair& lhs, kmer_t rhs) {
                             return (extract_kmer(lhs.seq) < rhs);
                             });
 
              up= std::upper_bound (arr.begin(), arr.end(), key,
                                   [] (kmer_t rhs, const KmerBucketPair& lhs) {
                             return (rhs < extract_kmer(lhs.seq));
                             });

               std::vector<read_t> read_list;
              for (auto itr=low; itr!= up; itr++) {
                  read_list.insert(std::end(read_list), std::begin((*itr).ReadIds), std::end((*itr).ReadIds));
              }
              
              // performing the sort and deletion of temp duplicates for now
              std::sort(read_list.begin(), read_list.end());
              read_list.erase(std::unique(read_list.begin(), read_list.end()), read_list.end());
              new_arr.push_back(KmerBucket{key, read_list});

              it = up;
          }

}

void GenerateSortedCSRListPerBucket(std::vector<KmerBucketPair>& arr,
                                    std::vector<read_t>& col_ind_readid,
                                    std::vector<kmer_t>& bucket_list,
                                    std::vector<int>& row_bucket_ptr,
                                    std::vector<int>& hist_bucket,
                                    uint64_t *num_outliers,
                                    uint64_t *num_buckets)
{

    std::vector<KmerBucketPair>::iterator low,up, it;
    int offset=0;
    uint64_t numout=0;
    uint64_t nbuckets=0;

#ifdef DEGUG_ISFORWARD
    std::string output_filename("readlist_p" + std::to_string(rank) + "_k" + std::to_string(KMER_LENGTH) + ".csv");
    FILE *f = fopen(output_filename.c_str(), "w");
    if (f == NULL)
    {
        printf("Error opening read list file!\n");
        exit(1);
    }
#endif

    for( it = arr.begin(); it != arr.end(); )
          {
              kmer_t kmerplusone = (*it).seq;
              kmer_t key = extract_kmer(kmerplusone);
              nbuckets++;
              low=std::lower_bound (arr.begin(), arr.end(), key,
                                   [] (const KmerBucketPair& lhs, kmer_t rhs) {
                             return (extract_kmer(lhs.seq) < rhs);
                             });
 
              up= std::upper_bound (arr.begin(), arr.end(), key,
                                   [] (kmer_t rhs, const KmerBucketPair& lhs) {
                             return (rhs < extract_kmer(lhs.seq));
                             });

              std::vector<read_t> read_list;
              for (auto itr=low; itr!= up; itr++) {
                  read_list.insert(std::end(read_list), std::begin((*itr).ReadIds), std::end((*itr).ReadIds));
              }
              
              // performing the sort and deletion of temp duplicates for now
              std::sort(read_list.begin(), read_list.end());
              read_list.erase(std::unique(read_list.begin(), read_list.end()), read_list.end());

              size_t nreads=read_list.size();

#ifdef PRINT_BUCKET_HISTOGRAM
              if (nreads<max_reads_per_bucket)
                  hist_bucket[nreads]++;
#endif
            
              if (nreads>1 && nreads<max_reads_per_bucket) 
              {
                  col_ind_readid.insert(std::end(col_ind_readid), std::begin(read_list), std::end(read_list));
                  row_bucket_ptr.push_back(offset);
                  bucket_list.push_back(key);
                  offset = (int)col_ind_readid.size();
              }
              else
                  numout++;

              it = up;
          }
    row_bucket_ptr.push_back(offset);

#ifdef DEGUG_ISFORWARD
    fclose(f);
#endif

    assert((bucket_list.size()+1) == row_bucket_ptr.size());
    *num_outliers = numout;
    *num_buckets = nbuckets;     

}

uint64_t CartProdforBlankList(std::vector<read_t> &a, std::vector<read_t> &b, 
                          std::vector< std::vector< std::pair<read_t, read_t>> >& output)
{
   uint64_t vecsize = a.size()*b.size();
   uint64_t num_pairs=0;

  for(size_t i = 0, s = b.size(); i < vecsize; ++i) {
      std::pair<read_t, read_t> l_pair = std::make_pair(a[i/s],b[i%s]);
      output[retrieve_read_proc_owner(l_pair.first)].push_back(l_pair);
      num_pairs++;
  }

  return num_pairs;

}

uint64_t CartProd(std::vector<read_t> &a, std::vector<read_t> &b, 
              std::vector< std::vector< std::pair<read_t, read_t>> >& output)
{
   size_t vecsize=a.size()*b.size();
   uint64_t num_pairs=0;
   
  for(size_t i = 0, s = b.size(); i < vecsize; ++i) {
      std::pair<read_t, read_t> l_pair = std::make_pair(a[i/s], b[i%s]);
      std::pair<read_t, read_t> r_pair = std::make_pair(b[i%s], a[i/s]);
      output[retrieve_read_proc_owner(l_pair.first)].push_back(l_pair);
      output[retrieve_read_proc_owner(r_pair.first)].push_back(r_pair);
      num_pairs+=2;
  }

  return num_pairs;

}


std::vector <std::vector< std::pair<read_t, read_t>> > GenerateEdgeListPerProcess(
                                      std::vector<KmerBucketPair>& arr,
                                      uint64_t *num_outliers,
                                      uint64_t *num_buckets)
{

    std::vector<KmerBucketPair>::iterator low,up, it;
    std::vector< std::vector< std::pair<read_t, read_t>> > edgeListperBucket(size);
    uint64_t num_pairs=0;
    uint64_t numout=0, nbuckets=0;

    for( it = arr.begin(); it != arr.end(); )
          {
              kmer_t kmerplusone = (*it).seq;
              kmer_t key = extract_kmer(kmerplusone);
              nbuckets++;
              low=std::lower_bound (arr.begin(), arr.end(), key,
                                   [] (const KmerBucketPair& lhs, kmer_t rhs) {
                             return (extract_kmer(lhs.seq) < rhs);
                             });
 
              up= std::upper_bound (arr.begin(), arr.end(), key,
                                   [] (kmer_t rhs, const KmerBucketPair& lhs) {
                             return (rhs < extract_kmer(lhs.seq));
                             });

              size_t nreads=0;
              for (auto itr=low; itr!= up; itr++) 
                   nreads += (*itr).ReadIds.size();

              if (nreads<max_reads_per_bucket) {          
 
              for (auto itr=low; itr!= up; itr++) {
                  for (auto itr2=itr+1; itr2!= up; itr2++) {
                       num_pairs += CartProd((*itr).ReadIds, (*itr2).ReadIds, edgeListperBucket);
                  }
                  if((*itr).pred_bp == 4) {
                      num_pairs += CartProdforBlankList((*itr).ReadIds, (*itr).ReadIds, edgeListperBucket);
                  }

                }
              }
              else {
                   numout++;
              }

              it = up;
          }
    
    uint64_t all_pairs=0;
    for (size_t i=0; i<size; i++)
         all_pairs += edgeListperBucket[i].size();

    assert(num_pairs == all_pairs);

    *num_outliers = numout;
    *num_buckets = nbuckets;
    return edgeListperBucket;


}

bool IsEqual (std::pair<read_t, read_t> &left, std::pair<read_t, read_t> &right)
     { return ((left.first == right.first) && (left.second == right.second)); }



void construct_sgraph_buffers(std::vector< std::pair<read_t, read_t> >& arr, uint64_t reads_disp, uint64_t num_reads_per_p)
{
    std::vector<std::pair<read_t, read_t>>::iterator low,up,it;
    int offset=0;

    row_bucket_ptr.resize(num_reads_per_p + 1);
    row_bucket_ptr[0]=0;

    for(read_t vid=0; vid<num_reads_per_p; vid++) {
        read_t gvid=reads_disp + vid;
        low=std::lower_bound (arr.begin(), arr.end(), gvid,
                [] (const std::pair<read_t, read_t>& lhs, read_t rhs) {
                return (lhs.first < rhs);
                });

        up= std::upper_bound (arr.begin(), arr.end(), gvid,
                [] (read_t rhs, const std::pair<read_t, read_t>& lhs) {
                return (rhs < lhs.first);
                });

        int num_neighbors = std::distance(low, up);
        offset += num_neighbors;
        row_bucket_ptr[vid+1] = offset;
    }
   assert(size_t(row_bucket_ptr[num_reads_per_p]) == arr.size());

    for( it = arr.begin(); it != arr.end(); it++)
        col_ind_readid.push_back((*it).second);

}

#ifdef DEBUG_BUCKET_PARMETIS
void print_local_edgelist_entries(std::vector< std::pair<read_t, read_t> >& arr)
{
    std::string output_filename(outputDirName + "/localreadpairlist_p" + std::to_string(rank) + ".csv");

    FILE *f = fopen(output_filename.c_str(), "w");
    if (f == NULL)
    {
        printf("rank: %d, Error opening local edgelist readpairs file!\n", rank);
        exit(1);
    }

    for (size_t i=0; i<arr.size(); i++)
         fprintf(f, "(%lu, %lu)\n", arr[i].first, arr[i].second);

    fclose(f);

}

void print_gathered_edgelist()
{
    std::string output_filename1(outputDirName + "/colidxlist_p" + std::to_string(rank) + ".csv");
    FILE *f1 = fopen(output_filename1.c_str(), "w");
    if (f1== NULL)
    {
        printf("rank: %d, Error opening colidxlist readpairs file!\n", rank);
        exit(1);
    }

    std::string output_filename2(outputDirName + "/rowbist_p" + std::to_string(rank) + ".csv");
    FILE *f2 = fopen(output_filename2.c_str(), "w");
    if (f2== NULL)
    {
        printf("rank: %d, Error opening rowb readpairs file!\n", rank);
        exit(1);
    }

    for (size_t i=0; i<col_ind_readid.size(); i++)
         fprintf(f1, "%lu\n", col_ind_readid[i]);

    for (size_t i=0; i<row_bucket_ptr.size(); i++)
         fprintf(f2, "%d\n", row_bucket_ptr[i]);


    fclose(f1);
    fclose(f2);

}
#endif

void transfer_edge_list (std::vector < std::vector< std::pair<read_t, read_t>> >& vec, uint64_t reads_disp, uint64_t num_reads_per_p)
{
    uint64_t ssize=0, rsize=0;
    std::vector<int> scounts(size,0);
    std::vector<int> rcounts(size,0);
    std::vector<int> rdisp (size,0);
    std::vector<int> sdisp (size,0);
    std::vector<std::pair<read_t, read_t>> pair_send_buf;

    for (int t=0; t<size; t++)
    {
        std::sort(vec[t].begin(), vec[t].end(), [](const std::pair<read_t, read_t> &left, const std::pair<read_t, read_t> &right) {
                                         return (left.first < right.first) || ( (left.first == right.first) && (left.second < right.second) );
                                         });

        vec[t].erase(std::unique(vec[t].begin(), vec[t].end(), IsEqual), vec[t].end());
        scounts[t]=vec[t].size();
        ssize += vec[t].size();
        pair_send_buf.insert(std::end(pair_send_buf), std::begin(vec[t]), std::end(vec[t]));
    }

    if (ssize > INT_MAX)
        fprintf(stderr, "Error!!! MPI cannot handle count sizes greater than INT_MAX as ssize: %lu. Please increase the total number of cores and try again!\n", ssize);
    assert(ssize <= INT_MAX);
    //TODO:
    //In case ssize > INT_MAX
    //Serialize the send buffer
    //For send counts, perform multiple alltoall calls to communicate the total number of read-pairs

    assert(ssize==pair_send_buf.size());

    for (int i=0; i<size; i++) {
         vec[i].clear();
         vec[i].shrink_to_fit();
    }

    // create contiguous derived data type
    MPI_Datatype rowtype;
    MPI_Type_contiguous(sizeof(std::pair<read_t, read_t>), MPI_BYTE, &rowtype);
    MPI_Type_commit(&rowtype);

    sdisp[0] = 0;
    for (int i=1; i<size; i++) sdisp[i] = scounts[i-1] + sdisp[i-1];

    MPI_Barrier(MPI_COMM_WORLD);

    MPI_Alltoall (scounts.data(), 1, MPI_INT, rcounts.data(), 1, MPI_INT, MPI_COMM_WORLD);

    for (int t=0; t<size; t++) rsize += rcounts[t];
    rdisp[0] = 0;
    for (int i=1; i<size; i++) rdisp[i] = rcounts[i-1] + rdisp[i-1];

    std::vector<std::pair<read_t, read_t>> pair_recv_buf (rsize);

    MPI_Barrier(MPI_COMM_WORLD);

    MPI_Alltoallv(pair_send_buf.data(), scounts.data(), sdisp.data(), rowtype,
                   pair_recv_buf.data(), rcounts.data(), rdisp.data(), rowtype, 
                   MPI_COMM_WORLD);

    pair_send_buf.clear();
    pair_send_buf.shrink_to_fit();
    scounts.clear();
    sdisp.clear();
    rcounts.clear();
    rdisp.clear();

    // free datatype
    MPI_Type_free(&rowtype);

    //Sort and erase duplicate read-pairs
    std::sort(pair_recv_buf.begin(), pair_recv_buf.end(), [](const std::pair<read_t, read_t> &left, const std::pair<read_t, read_t> &right) {
                                         return (left.first < right.first) || ( (left.first == right.first) && (left.second < right.second) );
                                         });

    pair_recv_buf.erase(std::unique(pair_recv_buf.begin(), pair_recv_buf.end(), IsEqual), pair_recv_buf.end());

#ifdef DEBUG_BUCKET_PARMETIS
    print_local_edgelist_entries(pair_recv_buf);
#endif

    construct_sgraph_buffers(pair_recv_buf, reads_disp, num_reads_per_p);

#ifdef DEBUG_BUCKET_PARMETIS
    print_gathered_edgelist();
#endif

    pair_recv_buf.clear();
    pair_recv_buf.shrink_to_fit();

}

void process_inplace_merge(std::vector<KmerBucketPair> &kmer_proc_buf,
                           std::vector<KmerBucketPair> &kmer_recv_buf)
{
     if (kmer_proc_buf.size())
    {
       size_t offset = kmer_proc_buf.size();
       kmer_proc_buf.insert(kmer_proc_buf.end(), kmer_recv_buf.begin(), kmer_recv_buf.end());

       std::inplace_merge(kmer_proc_buf.begin(),
                  kmer_proc_buf.begin()+offset,
                  kmer_proc_buf.end());
    }
    else
      kmer_proc_buf.insert(kmer_proc_buf.end(), kmer_recv_buf.begin(), kmer_recv_buf.end());

}


void transfer_kmers (std::vector<int> &scounts_kmer, 
                     std::vector < std::vector<KmerBucket> > &kmer_send_buf,
                     std::vector<int> &send_num_blanks) 
{

    std::vector<int> rcounts_kmer (size,0);
    std::vector<int> recv_num_blanks (size,0);

    /* Alltoall to communicate the number of KmerBuckets with blank preceeding base-pair */
    MPI_Barrier(MPI_COMM_WORLD);
    double tblank1  = MPI_Wtime ();
    MPI_Alltoall (send_num_blanks.data(), 1, MPI_INT, recv_num_blanks.data(), 1, MPI_INT, MPI_COMM_WORLD);
    double tblank2  = MPI_Wtime ();
    alltoall_time += (tblank2 - tblank1);

    /*Alltoall to calculate the total number of KmerBuckets being communicated */
    MPI_Barrier(MPI_COMM_WORLD);
    double t7 = MPI_Wtime ();
    MPI_Alltoall (scounts_kmer.data(), 1, MPI_INT, rcounts_kmer.data(), 1, MPI_INT, MPI_COMM_WORLD);
    double t8 = MPI_Wtime ();
    alltoall_time += (t8 - t7);

    /*buffers to communicate serialized data */
    std::string send_buffer;
    std::string tmp_buffer;
    std::string recv_buffer;
    char* recv_ptr=nullptr;

    /* send and recv buffers for obtaining the serialized data */
     std::vector<uint64_t> scounts(size,0); //sending serialized data in bytes
     std::vector<uint64_t> rcounts (size,0);
     std::vector<uint64_t> rdisp (size,0);
     std::vector<int> scounts_dd(size,0);
     std::vector<int> rcounts_dd(size,0);
     std::vector<int> sdisp_dd(size,0);
     std::vector<int> rdisp_dd(size,0);

     //serialize
     std::ostringstream os_cnt(std::ios::binary | std::ios::out | std::ios::in);
     uint64_t ssize=0,rsize=0,r_mnodes=0;
     uint64_t blank_num_buckets=0, all_num_buckets=0;
     uint64_t pad_width=1000;

     for (int i=0; i<size; i++)
         {
             os_cnt.str("");
             tmp_buffer.clear();
             assert(kmer_send_buf[i].size() == scounts_kmer[i]);
             for (size_t j=0; j<kmer_send_buf[i].size(); j++) {
                 serialize(os_cnt, kmer_send_buf[i][j]);
             }

             tmp_buffer = os_cnt.str();
             if (tmp_buffer.size())
             {
                 size_t nonPaddedSize = tmp_buffer.size();
                 size_t padded_size = ((nonPaddedSize/pad_width) + (nonPaddedSize%pad_width!=0))*pad_width;

                 tmp_buffer.resize(padded_size, 0);
             }
             scounts[i] = tmp_buffer.length();
             scounts_dd[i] = scounts[i]/pad_width;
             send_buffer.append(tmp_buffer);
         }

    //clear the buffers
    tmp_buffer.clear();
    for (int i=0; i<size; i++) {
        kmer_send_buf[i].clear();
        kmer_send_buf[i].shrink_to_fit();
    }

    for (int t=0; t<size; t++) ssize += scounts[t];
    assert (ssize == send_buffer.length());

    MPI_Barrier(MPI_COMM_WORLD);
    double comm1 = MPI_Wtime ();
    //Sending the serialized data buffer
    MPI_Alltoall(scounts.data(), 1, MPI_UINT64_T, rcounts.data(), 1, MPI_UINT64_T, MPI_COMM_WORLD);
    double comm2 = MPI_Wtime ();
    alltoall_time += (comm2 - comm1);
    
    for (int t=0; t<size; t++) {
        rsize += rcounts[t];
        blank_num_buckets += recv_num_blanks[t];
        all_num_buckets += rcounts_kmer[t];
    }


     for (int t=0; t<size; t++) {
          sdisp_dd[t] = (t>0) ? (scounts_dd[t-1] + sdisp_dd[t-1]) : 0;
          rdisp[t] = (t>0) ? (rcounts[t-1] + rdisp[t-1]) : 0;
	      rcounts_dd[t] = rcounts[t]/pad_width;
	      rdisp_dd[t] = (t>0) ? (rcounts_dd[t-1] + rdisp_dd[t-1]) : 0;
     }

     recv_buffer.resize(rsize, 'F');
     recv_ptr=&recv_buffer[0];
     MPI_Datatype rowtype;

     //create contiguous derived data type
     MPI_Type_contiguous(pad_width, MPI_BYTE, &rowtype);
     MPI_Type_commit(&rowtype);

     MPI_Barrier(MPI_COMM_WORLD);

    double t9 = MPI_Wtime ();
    int result = MPI_Alltoallv(send_buffer.c_str(), scounts_dd.data(), sdisp_dd.data(), rowtype,
                   recv_ptr, rcounts_dd.data(), rdisp_dd.data(), rowtype, 
                   MPI_COMM_WORLD);

    if (result != MPI_SUCCESS) {
         printf("rank: %d, MPI_Alltoallv in Transfer_kmers failed with return value: %d\n", rank, result);
         MPI_Finalize();
         exit(2);
     }
    double t10 = MPI_Wtime ();
    alltoallv_time += (t10 - t9);

    // free datatype
     MPI_Type_free(&rowtype);

    //clear buffers
     os_cnt.str("");
     scounts.clear();
     scounts_dd.clear();
     sdisp_dd.clear();
     send_buffer.clear();
     scounts_kmer.clear();
     send_num_blanks.clear();
     rcounts_dd.clear();
     rdisp_dd.clear();

    //deserialize into temporary buffer
    double t21 = MPI_Wtime ();

    std::vector<KmerBucketPair> kmer_recv_buf;
    std::istringstream is(std::ios::binary | std::ios::out | std::ios::in);
    uint64_t blank_idx=0;
    
    for (int i=0; i<size; i++)
      {
          is.rdbuf()->pubsetbuf(const_cast<char*>(recv_buffer.c_str()+rdisp[i]), rcounts[i]);
          int nnodes = rcounts_kmer[i];
          int blank_nodes = recv_num_blanks[i];
          int idx=0;

          while(nnodes)
          {
              KmerBucket kb;
              deserialize(is, kb);

              if (idx<blank_nodes) {    
                  kmer_recv_buf.push_back(KmerBucketPair{kb.seq, kb.ReadIds, char_to_el('H')});
                  blank_idx++;
              }
              else { 
                  kmer_recv_buf.push_back(KmerBucketPair{kb.seq, kb.ReadIds, (ElType)(extract_pred(kb.seq) & 0x7)});
              }
              idx++;
              nnodes--;

          }
          assert(idx == rcounts_kmer[i]);
      }
      MPI_Barrier(MPI_COMM_WORLD);
      assert(all_num_buckets == kmer_recv_buf.size());
      assert(blank_num_buckets == blank_idx);

      recv_buffer.clear();
      rcounts.clear();
      rdisp.clear();
      recv_num_blanks.clear();
      rcounts_kmer.clear();  

    // sort the recv buffer
    double tm1 = MPI_Wtime ();
    
    std::sort(kmer_recv_buf.begin(), kmer_recv_buf.end());

    double tm2 = MPI_Wtime ();
    unpack_rbuf_sort += (tm2 - tm1);

    // Insert and sort
    double tm3 = MPI_Wtime ();
    process_inplace_merge(kmer_proc_buf, kmer_recv_buf);

    double tm4 = MPI_Wtime ();
    unpack_rbuf_insert += (tm4 - tm3);

    kmer_recv_buf.clear();
    kmer_recv_buf.shrink_to_fit();

    // Combine and aggregate the read_id lists for each k-mer bucket
    double tm5 = MPI_Wtime ();

    SortAndAggregatePairBucket(kmer_proc_buf);
    double tm6 = MPI_Wtime ();
    unpack_rbuf_acc += (tm6 - tm5);

    double t22 = MPI_Wtime ();
    unpack_rbuf_time += (t22-t21);
    num_batch_transfers++;

}


void recalculate_min_lmer_opt (kmer_t kmer_in, lmer_t *m_lmer, lmer_t *m_lmer_freq, int *m_pos)
{
    lmer_t min_lmer=0, tmp_lmer=0;
    lmer_t min_lmer_freq=0, tmp_lmer_freq=0;
    int min_pos=0, k=0, j=0;
    lmer_t lmer_out=0;

    for (k=0; ((KMER_LENGTH-1) - k) >= (LMER_LENGTH-1); k++) {
       if (k==0) { 
          for(int j=k; j<LMER_LENGTH+k; j++) {
              lmer_out = kmer_to_lmer (kmer_in, j, lmer_out);
          }
       }
       else {
           for(; j<LMER_LENGTH+k; j++) {
              lmer_out = kmer_to_lmer (kmer_in, j, lmer_out);
           }
       }

        tmp_lmer = lmer_out;
        tmp_lmer_freq = global_lmer_frequency[tmp_lmer];

        if (k == 0) {
            min_lmer = tmp_lmer;
            min_lmer_freq = tmp_lmer_freq;
            min_pos = 0;
        }
        else {
           if (tmp_lmer_freq < min_lmer_freq) {
               min_lmer = tmp_lmer;
               min_lmer_freq = tmp_lmer_freq;
               min_pos = k;
           }
        }
    }
    assert (k == (KMER_LENGTH-LMER_LENGTH+1));

    *m_lmer = min_lmer;
    *m_lmer_freq = min_lmer_freq;
    *m_pos = min_pos;
}

void recalculate_min_lmer (kmer_t kmer_in, lmer_t *m_lmer, lmer_t *m_lmer_freq, int *m_pos)
{
    lmer_t min_lmer=0, tmp_lmer=0;
    lmer_t min_lmer_freq=0, tmp_lmer_freq=0;
    int min_pos=0, k=0;

    for (k=0; ((KMER_LENGTH-1) - k) >= (LMER_LENGTH-1); k++) {
        lmer_t lmer_out=0;
        for(int j=k; j<LMER_LENGTH+k; j++) {
            lmer_out = kmer_to_lmer (kmer_in, j, lmer_out);
        }

        tmp_lmer = lmer_out;
        tmp_lmer_freq = global_lmer_frequency[tmp_lmer];

        if (k == 0) {
            min_lmer = tmp_lmer;
            min_lmer_freq = tmp_lmer_freq;
            min_pos = 0;
        }
        else {
           if (tmp_lmer_freq < min_lmer_freq) {
               min_lmer = tmp_lmer;
               min_lmer_freq = tmp_lmer_freq;
               min_pos = k;
           }
        }
    }
    assert (k == (KMER_LENGTH-LMER_LENGTH+1));

    *m_lmer = min_lmer;
    *m_lmer_freq = min_lmer_freq;
    *m_pos = min_pos;
}

void Sliding_window_l_variableReadLength (input_read_data & reads) {

#ifdef LMER_DEBUG2
    ElType this_alpha;
    char kmer_out[lmer_len+1];
#endif

    for (size_t idx=0; idx<reads.num_of_reads; idx++) {

        size_t start = reads.read_data_ptr[idx];
        size_t end = reads.read_data_ptr[idx+1];

        if(end - start < LMER_LENGTH) continue; /*too short a read*/

        lmer_t kmer = 0;

        size_t p=start;
        for(size_t i=0; p<end && i<LMER_LENGTH-1; i++) {
            kmer = lmer_shift(kmer, char_to_el(reads.read_data[p++]));
        }

        while(p<end) {
            kmer = lmer_shift(kmer, char_to_el(reads.read_data[p++]));

#ifdef LMER_DEBUG2
            printf("lmer: %lu,");
            for(int j=0; j<LMER_LENGTH; j++) {
                this_alpha = lmerel (kmer, j);
                lmer_out[j] = el_to_char(this_alpha);
            }
            lmer_out[lmer_len] = '\0';
            printf(" lmer_out: %s \n", lmer_out);
#endif

            lmer_frequency[kmer]++;
        }
    }
}


void Sliding_window_l (const char *ptr, size_t length) {

#ifdef LMER_DEBUG2
      ElType this_alpha;
      char kmer_out[lmer_len+1];
#endif

  size_t p=0;

  /*find start of a read*/
  for(; ptr[p]!='>' && p<length; p++) {/*noop*/ }

  lmer_t kmer = 0;

  while(p<length) {
    assert(ptr[p]=='>'); /*this will be true*/

    /*skip till newline*/
    for(; p<length && ptr[p]!='\n'; p++) {/*noop*/ }
    p++; /*skip the newline*/

    if(p+LMER_LENGTH > length) break; /*too short a read*/
    kmer = 0;
    int i;
    for(i=0; ptr[p]!='\n' && i<LMER_LENGTH-1; i++) {
      kmer = lmer_shift(kmer, char_to_el(ptr[p++]));
    }

    while(p<length && ptr[p]!='\n') {
      kmer = lmer_shift(kmer, char_to_el(ptr[p++]));

#ifdef LMER_DEBUG2
      printf("lmer: %lu,");
      for(int j=0; j<LMER_LENGTH; j++) {
          this_alpha = lmerel (kmer, j);
          lmer_out[j] = el_to_char(this_alpha);
      }
      lmer_out[lmer_len] = '\0';
      printf(" lmer_out: %s \n", lmer_out);
#endif

      lmer_frequency[kmer]++;

    }
    p++; /*skip the newline*/
  }

}

void check_kmer_revc(kmer_t fkmer, kmer_t rkmer)
{
    kmer_t tkmer=0;

    for (int k=0; k<KMER_LENGTH; k++) {
        tkmer = kmer_shift_rc(tkmer, rev_comp(el_to_char(kmerel(fkmer, k))));
    }
    if (tkmer != rkmer)
        printf("Error!, k-mers are not rev complements of each other: f:%lu, r:%lu, t:%lu\n",fkmer, rkmer, tkmer);
    assert(tkmer == rkmer);
}

size_t Sliding_window (const char *ptr, size_t length, int *n_kmers, 
                     std::vector<std::vector<KmerPair>> &partial_kmer_counts,
                     std::vector<std::vector<KmerPair>> &partial_blank_kmer_counts,
                     uint64_t reads_disp, size_t read_start_idx, bool isForward)
{
#ifdef LMER_DEBUG
    FILE *fp_d;
    char debug_file_name[25];
    char proc_id[3];

    sprintf(proc_id, "%d", rank); 
    strcpy(debug_file_name,"debug_p");
    strcpy(&debug_file_name[strlen(debug_file_name)],proc_id);
    strcpy(&debug_file_name[strlen(debug_file_name)],".log");
    fp_d = fopen (debug_file_name, "w");
    
    /*check to see if it opened okay */
    if (fp_d == NULL)
    {
		printf ("Error opening proc %d 's dump file \n", rank);
		exit (0);
    }
#endif

  size_t p=0;
  std::vector<int> scounts_kmer (size,0);
  std::vector<int> num_blanks (size, 0);
  std::vector< std::vector<KmerBucket> > kmer_send_buf(size);
  int kpos=0;
  size_t num_reads=0;
  if (rank==0) { printf("rank: %d, reads_size: %lu\n", rank, reads_data.size()); }

  /*find start of a read*/
  for(; ptr[p]!='>' && p<length; p++) {/*noop*/ }

  int num_kmers=*n_kmers;
  read_t read_id=reads_disp;
  int kmers_per_read=0;
  kmer_t kmer = 0;
#ifdef BIDIRECT 
  kmer_t kmer_rc=0;
#endif 
  lmer_t lmer_out = 0, min_lmer=0, tmp_lmer=0;
  uint64_t min_lmer_freq=0, tmp_lmer_freq=0;
  int min_pos=0, tmp_pos=0;

  while(p<length) {
    /* start of a new read */
    assert(ptr[p]=='>'); /*this will be true*/

    /*skip till newline*/
    for(; p<length && ptr[p]!='\n'; p++) {/*noop*/ }
    p++; /*skip the newline*/

    if(p+KMER_LENGTH > length) break; /*too short a read*/

    kmer=0, lmer_out=0;
    min_lmer=0, min_lmer_freq=0;
    tmp_lmer=0, tmp_lmer_freq=0;
    min_pos=0, tmp_pos=0;
    kmers_per_read=0;
    int i=0;
 #ifdef BIDIRECT    
    int r=0;
    size_t rd_off=0, rc_idx=0;
    std::vector<RevKmerPair> kmer_counter;
    size_t read_src=0, rptr=0;
    kmer_rc=0;
 #endif

    std::vector<read_t> tmp_read;
    int read_len=0;

 #ifdef BIDIRECT   
    read_src=p;
 #endif
    int counter = 0;
    for(i=0; ptr[p]!='\n' && i<KMER_LENGTH-1; i++) {
      ElType el;
      if (ptr[p] == 'N') {
          // set a counter that will count down and not output anything until the count down is finished
          // but continue to add new elements to tmp_read

          // reset the counter anytime an N is read
          el = char_to_el('A'); // FIXME: also record where it is to recover it later
          counter = KMER_LENGTH;
      } else {
          if (counter>0) counter--;
          el = char_to_el(ptr[p]);
      }

      kmer = kmer_shift(kmer, el);
      push_back_bp_w(char_to_el(ptr[p]), &read_len, tmp_read);

      if (i<LMER_LENGTH-1) 
          lmer_out = lmer_shift(lmer_out, el);
      else {
           lmer_out = lmer_shift(lmer_out, el);

           tmp_lmer = lmer_out;
           tmp_lmer_freq = global_lmer_frequency[tmp_lmer];
           tmp_pos = i-(LMER_LENGTH-1);

           if (i == LMER_LENGTH-1) {
               min_lmer = tmp_lmer;
               min_lmer_freq = tmp_lmer_freq;
           }
           else {
                if (tmp_lmer_freq < min_lmer_freq) {
                    min_lmer = tmp_lmer;
                    min_lmer_freq = tmp_lmer_freq;
                    min_pos = tmp_pos;
                }
           }
      }
      p++;
#ifdef BIDIRECT 
      rptr++;
#endif
    }

    int num_kmers_added = 0;
    while(p<length && ptr[p]!='\n') {

      ElType el;
      if (ptr[p] == 'N') {
          // set a counter that will count down and not output anything until the count down is finished
          // but continue to add new elements to tmp_read

          // reset the counter anytime an N is read
          el = char_to_el('A'); // FIXME: also record where it is to recover it later
          counter = KMER_LENGTH;
      } else {
          if (counter>0) counter--;
          el = char_to_el(ptr[p]);
      }



      if (!kmers_per_read)
         kmer = kmer_shift(kmer, el);
      else
         kmer = kmerplusone_shift(kmer, el);

      lmer_out = lmer_shift(lmer_out, el);
      uint64_t lmer_out_freq = global_lmer_frequency[lmer_out];
      push_back_bp_w(char_to_el(ptr[p]), &read_len, tmp_read);
      
      if (min_pos < 0) {
          recalculate_min_lmer(extract_kmer(kmer), &min_lmer, &min_lmer_freq, &min_pos);
          num_recalculate_lmer++;
      }

      if (lmer_out_freq < min_lmer_freq) {
          min_lmer = lmer_out;
          min_lmer_freq = lmer_out_freq;
          min_pos = KMER_LENGTH-LMER_LENGTH;
      }
      p++;
      min_pos--;
#ifdef BIDIRECT 
      rptr++;
#endif
      
#ifdef LMER_DEBUG
      fprintf(fp_d, "kmer: %lu, lmer: %lu, lmer_freq: %lu, min_lmer: %lu, min_freq: %lu\n", kmer,lmer_out,lmer_out_freq,min_lmer,min_lmer_freq);
#endif

      double T1 = MPI_Wtime();

      if (counter == 0) {
#ifdef BIDIRECT 
      if (!kmers_per_read) {
           kmer_counter.push_back(RevKmerPair{kmer, min_lmer, true});
      } else {
           kmer_counter.push_back(RevKmerPair{kmer, min_lmer, false});
      }
#else
      if (!kmers_per_read) {
          partial_blank_kmer_counts[retrieve_proc_id(min_lmer)].push_back(KmerPair{kmer,read_id});
      } else {
          partial_kmer_counts[retrieve_proc_id(min_lmer)].push_back(KmerPair{kmer,read_id});
      }
      
#endif
      num_kmers_added ++;
      }

      double T2 = MPI_Wtime();
      vec_insert_time += (T2-T1);

      kmers_per_read++;

      
    } // end of while loop, end of a read
    num_kmers += num_kmers_added;

#ifdef BIDIRECT
    assert(rptr==read_length); 
    if ((p-read_src)!=read_length) 
         printf("rank: %d, read_id: %lu, read_len: %lu, p: %lu, read_src: %lu, rptr: %lu\n", 
                 rank, read_id, (p-read_src), p, read_src, rptr);
    assert((p-read_src)==read_length);
    std::string str{ptr+read_src, ptr+(rptr+read_src)};
    std::string str_rev(str);
    std::reverse(str_rev.begin(), str_rev.end());
    assert(str_rev.size()==read_length);

    //initialize
    kmer_rc=0, lmer_out=0;
    min_lmer=0, min_lmer_freq=0;
    tmp_lmer=0, tmp_lmer_freq=0;
    min_pos=0, tmp_pos=0;
    r=0;
    rd_off=0, rc_idx=0;
    int kmers_in_read=0;

    counter = 0;
    for(int j=0; str_rev[r]!='\n' && j<KMER_LENGTH-1; j++) {
      ElType el;
      if (str_rev[r] == 'N') {
          // set a counter that will count down and not output anything until the count down is finished
          // but continue to add new elements to tmp_read

          // reset the counter anytime an N is read
          el = rev_comp('A'); // FIXME: also record where it is to recover it later
          counter = KMER_LENGTH;
      } else {
          if (counter>0) counter--;
          el = rev_comp(str_rev[r]);
      }


      kmer_rc = kmer_shift(kmer_rc, el);

      if (j<LMER_LENGTH-1) 
          lmer_out = lmer_shift(lmer_out, el);
      else {
           lmer_out = lmer_shift(lmer_out, el);

           tmp_lmer = lmer_out;
           tmp_lmer_freq = global_lmer_frequency[tmp_lmer];
           tmp_pos = j-(LMER_LENGTH-1);

           if (j == LMER_LENGTH-1) {
               min_lmer = tmp_lmer;
               min_lmer_freq = tmp_lmer_freq;
           }
           else {
                if (tmp_lmer_freq < min_lmer_freq) {
                    min_lmer = tmp_lmer;
                    min_lmer_freq = tmp_lmer_freq;
                    min_pos = tmp_pos;
                }
           }
      }
      r++;
    }

    while(r<str_rev.size() && str_rev[r]!='\n') {
      ElType el;
      if (str_rev[r] == 'N') {
          // set a counter that will count down and not output anything until the count down is finished
          // but continue to add new elements to tmp_read

          // reset the counter anytime an N is read
          el = rev_comp('A'); // FIXME: also record where it is to recover it later
          counter = KMER_LENGTH;
      } else {
          if (counter>0) counter--;
          el = rev_comp(str_rev[r]);
      }

      if (!kmers_in_read)
         kmer_rc = kmer_shift(kmer_rc, el);
      else
         kmer_rc = kmerplusone_shift(kmer_rc, el);

      lmer_out = lmer_shift(lmer_out, el);
      uint64_t lmer_out_freq = global_lmer_frequency[lmer_out];
      
      if (min_pos < 0) {
          recalculate_min_lmer(extract_kmer(kmer_rc), &min_lmer, &min_lmer_freq, &min_pos);
          num_recalculate_lmer++;
      }

      if (lmer_out_freq < min_lmer_freq) {
          min_lmer = lmer_out;
          min_lmer_freq = lmer_out_freq;
          min_pos = KMER_LENGTH-LMER_LENGTH;
      }
      r++;
      min_pos--;
      if (counter == 0) {
      rc_idx = num_kmers_added-rd_off-1;

      RevKmerPair this_pair = (!kmers_in_read)? RevKmerPair{kmer_rc, min_lmer, true}: RevKmerPair{kmer_rc, min_lmer, false}; ;
#ifdef CHECK_BIDIRECT
      check_kmer_revc(extract_kmer(kmer_counter[rc_idx].fkmer), extract_kmer(kmer_rc));
#endif      
      RevKmerPair krd_pair = (extract_kmer(kmer_counter[rc_idx].fkmer)<extract_kmer(kmer_rc)) ? kmer_counter[rc_idx] : this_pair;
      if (krd_pair.is_blank) 
           partial_blank_kmer_counts[retrieve_proc_id(krd_pair.min_lmer)].push_back(KmerPair{krd_pair.fkmer, read_id});
      else
           partial_kmer_counts[retrieve_proc_id(krd_pair.min_lmer)].push_back(KmerPair{krd_pair.fkmer, read_id});
      
      rd_off++;
      }
      kmers_in_read++;
    }

#endif

    if (num_kmers > MAX_KMER_COUNT){
          double t25 = MPI_Wtime ();
          // initiate collective communication to pass k-mers and their respective counts to rightful owners
          // calculate global owner of each k-mer and populate the k-mer and count to their respective 'p'th local buffers
          // if global position of a k-mer in my k_map != me, delete k-mer from my k_map
          // iterate through k-mers recieved after collective communication ends, and add k-mers to my k_map 
          // reset num_kmers count to 0.
          
          double T5 = MPI_Wtime();

          for (int t=0; t<size; t++)
          {
            kpos=0;
            sort(partial_blank_kmer_counts[t].begin(), partial_blank_kmer_counts[t].end());
            
            int nblank = SortAndAggregatePair(partial_blank_kmer_counts[t], kmer_send_buf[t]);
            kpos += nblank;
            num_blanks[t] = nblank;
            partial_blank_kmer_counts[t].clear();
            partial_blank_kmer_counts[t].shrink_to_fit();

            sort(partial_kmer_counts[t].begin(), partial_kmer_counts[t].end());
            
            kpos += SortAndAggregatePair(partial_kmer_counts[t], kmer_send_buf[t]);
            partial_kmer_counts[t].clear();
            partial_kmer_counts[t].shrink_to_fit();
            scounts_kmer[t] = kpos;
          }
          

          double T6 = MPI_Wtime();
          tmap_insert_time += (T6-T5);

          transfer_kmers (scounts_kmer, kmer_send_buf, num_blanks);
          num_kmers = 0;

          scounts_kmer.clear();
          num_blanks.clear();
          double t26 = MPI_Wtime ();
          sl_win_time_int += (t26-t25);
       } // end of if condition

    assert(kmers_per_read == (read_length-KMER_LENGTH+1));
    assert(read_len == read_length);
    assert(tmp_read.size() == ceil((double)((double)read_length/nels_per_value_w)));
    assert((read_start_idx*read_size)<reads_data.size());

    std:copy(tmp_read.begin(), tmp_read.end(), reads_data.begin()+(read_start_idx*read_size));

    read_id++;
    read_start_idx+=1+ispaired;
    num_reads++;

#ifdef LMER_DEBUG
    fprintf(fp_d, "----------------------------\n");
#endif
    p++; /*skip the newline*/
  } //end of processing all reads for a given process

#ifdef LMER_DEBUG
          fclose (fp_d);
#endif

  *n_kmers = num_kmers;
  scounts_kmer.shrink_to_fit();

  return num_reads;

}

size_t Sliding_window_variableReadLength (input_read_data &rdata, int *n_kmers, 
                     std::vector<std::vector<KmerPair>> &partial_kmer_counts,
                     std::vector<std::vector<KmerPair>> &partial_blank_kmer_counts,
                     uint64_t reads_disp, bool isForward)
{
    const char *ptr = rdata.read_data;

#ifdef LMER_DEBUG
    FILE *fp_d;
    char debug_file_name[25];
    char proc_id[3];

    sprintf(proc_id, "%d", rank); 
    strcpy(debug_file_name,"debug_p");
    strcpy(&debug_file_name[strlen(debug_file_name)],proc_id);
    strcpy(&debug_file_name[strlen(debug_file_name)],".log");
    fp_d = fopen (debug_file_name, "w");

    /*check to see if it opened okay */
    if (fp_d == NULL)
    {
        printf ("Error opening proc %d 's dump file \n", rank);
        exit (0);
    }
#endif

    size_t p=0;
    std::vector<int> scounts_kmer (size,0);
    std::vector<int> num_blanks (size, 0);
    std::vector< std::vector<KmerBucket> > kmer_send_buf(size);
    int kpos=0;
    size_t num_reads=0;
    if (rank==0) { printf("rank: %d, reads_size: %lu\n", rank, reads_data.size()); }

    int num_kmers=*n_kmers;
    read_t read_id=reads_disp;
    int kmers_per_read=0;
    kmer_t kmer = 0;
#ifdef BIDIRECT 
    kmer_t kmer_rc=0;
#endif 
    lmer_t lmer_out = 0, min_lmer=0, tmp_lmer=0;
    uint64_t min_lmer_freq=0, tmp_lmer_freq=0;
    int min_pos=0, tmp_pos=0;

    for (size_t idx=0; idx<rdata.num_of_reads; idx++) {

        size_t start = rdata.read_data_ptr[idx];
        size_t end = rdata.read_data_ptr[idx+1];
        size_t cur_read_length = end - start;
        size_t p = start;

        kmer=0, lmer_out=0;
        min_lmer=0, min_lmer_freq=0;
        tmp_lmer=0, tmp_lmer_freq=0;
        min_pos=0, tmp_pos=0;
        kmers_per_read=0;
#ifdef BIDIRECT    
        std::vector<RevKmerPair> kmer_counter;
        kmer_rc=0;
#endif
        std::vector<read_t> tmp_read;
        int read_len=0;

        int counter = 0;
        for(size_t i=0; p<end && i<KMER_LENGTH-1; i++) {
            ElType el;
            if (ptr[p] == 'N') {
                // set a counter that will count down and not output anything until the count down is finished
                // but continue to add new elements to tmp_read

                // reset the counter anytime an N is read
                el = char_to_el('A'); // FIXME: also record where it is to recover it later
                counter = KMER_LENGTH;
            } else {
                if (counter>0) counter--;
                el = char_to_el(ptr[p]);
            }

            kmer = kmer_shift(kmer, el);
            push_back_bp_w(char_to_el(ptr[p]), &read_len, tmp_read);

            if (i<LMER_LENGTH-1) 
                lmer_out = lmer_shift(lmer_out, el);
            else {
                lmer_out = lmer_shift(lmer_out, el);

                tmp_lmer = lmer_out;
                tmp_lmer_freq = global_lmer_frequency[tmp_lmer];
                tmp_pos = i-(LMER_LENGTH-1);

                if (i == LMER_LENGTH-1) {
                    min_lmer = tmp_lmer;
                    min_lmer_freq = tmp_lmer_freq;
                }
                else {
                    if (tmp_lmer_freq < min_lmer_freq) {
                        min_lmer = tmp_lmer;
                        min_lmer_freq = tmp_lmer_freq;
                        min_pos = tmp_pos;
                    }
                }
            }
            p++;
        }

        int num_kmers_added = 0;
        while(p<end) {

            ElType el;
            if (ptr[p] == 'N') {
                // set a counter that will count down and not output anything until the count down is finished
                // but continue to add new elements to tmp_read

                // reset the counter anytime an N is read
                el = char_to_el('A'); // FIXME: also record where it is to recover it later
                counter = KMER_LENGTH;
            } else {
                if (counter>0) counter--;
                el = char_to_el(ptr[p]);
            }



            if (!kmers_per_read)
                kmer = kmer_shift(kmer, el);
            else
                kmer = kmerplusone_shift(kmer, el);

            lmer_out = lmer_shift(lmer_out, el);
            uint64_t lmer_out_freq = global_lmer_frequency[lmer_out];
            push_back_bp_w(char_to_el(ptr[p]), &read_len, tmp_read);

            if (min_pos < 0) {
                recalculate_min_lmer(extract_kmer(kmer), &min_lmer, &min_lmer_freq, &min_pos);
                num_recalculate_lmer++;
            }

            if (lmer_out_freq < min_lmer_freq) {
                min_lmer = lmer_out;
                min_lmer_freq = lmer_out_freq;
                min_pos = KMER_LENGTH-LMER_LENGTH;
            }
            p++;
            min_pos--;

#ifdef LMER_DEBUG
            fprintf(fp_d, "kmer: %lu, lmer: %lu, lmer_freq: %lu, min_lmer: %lu, min_freq: %lu\n", kmer,lmer_out,lmer_out_freq,min_lmer,min_lmer_freq);
#endif

            double T1 = MPI_Wtime();

            if (counter == 0) {
#ifdef BIDIRECT 
                if (!kmers_per_read) {
                    kmer_counter.push_back(RevKmerPair{kmer, min_lmer, true});
                } else {
                    kmer_counter.push_back(RevKmerPair{kmer, min_lmer, false});
                }
#else
                if (!kmers_per_read) {
                    partial_blank_kmer_counts[retrieve_proc_id(min_lmer)].push_back(KmerPair{kmer,read_id});
                } else {
                    partial_kmer_counts[retrieve_proc_id(min_lmer)].push_back(KmerPair{kmer,read_id});
                }

#endif
                num_kmers_added ++;
            }

            double T2 = MPI_Wtime();
            vec_insert_time += (T2-T1);

            kmers_per_read++;


        } // end of while loop, end of a read
        num_kmers += num_kmers_added;

#ifdef BIDIRECT
        std::string str{ptr+start, ptr+end};
        std::string str_rev(str);
        std::reverse(str_rev.begin(), str_rev.end());
        assert(str_rev.size()==cur_read_length);

        //initialize
        kmer_rc=0, lmer_out=0;
        min_lmer=0, min_lmer_freq=0;
        tmp_lmer=0, tmp_lmer_freq=0;
        min_pos=0, tmp_pos=0;
        size_t r=0;
        size_t rd_off=0, rc_idx=0;
        int kmers_in_read=0;

        counter = 0;
        for(int j=0; r < cur_read_length && j<KMER_LENGTH-1; j++) {
            ElType el;
            if (str_rev[r] == 'N') {
                // set a counter that will count down and not output anything until the count down is finished
                // but continue to add new elements to tmp_read

                // reset the counter anytime an N is read
                el = rev_comp('A'); // FIXME: also record where it is to recover it later
                counter = KMER_LENGTH;
            } else {
                if (counter>0) counter--;
                el = rev_comp(str_rev[r]);
            }


            kmer_rc = kmer_shift(kmer_rc, el);

            if (j<LMER_LENGTH-1) 
                lmer_out = lmer_shift(lmer_out, el);
            else {
                lmer_out = lmer_shift(lmer_out, el);

                tmp_lmer = lmer_out;
                tmp_lmer_freq = global_lmer_frequency[tmp_lmer];
                tmp_pos = j-(LMER_LENGTH-1);

                if (j == LMER_LENGTH-1) {
                    min_lmer = tmp_lmer;
                    min_lmer_freq = tmp_lmer_freq;
                }
                else {
                    if (tmp_lmer_freq < min_lmer_freq) {
                        min_lmer = tmp_lmer;
                        min_lmer_freq = tmp_lmer_freq;
                        min_pos = tmp_pos;
                    }
                }
            }
            r++;
        }

        while(r<cur_read_length) {
            ElType el;
            if (str_rev[r] == 'N') {
                // set a counter that will count down and not output anything until the count down is finished
                // but continue to add new elements to tmp_read

                // reset the counter anytime an N is read
                el = rev_comp('A'); // FIXME: also record where it is to recover it later
                counter = KMER_LENGTH;
            } else {
                if (counter>0) counter--;
                el = rev_comp(str_rev[r]);
            }

            if (!kmers_in_read)
                kmer_rc = kmer_shift(kmer_rc, el);
            else
                kmer_rc = kmerplusone_shift(kmer_rc, el);

            lmer_out = lmer_shift(lmer_out, el);
            uint64_t lmer_out_freq = global_lmer_frequency[lmer_out];

            if (min_pos < 0) {
                recalculate_min_lmer(extract_kmer(kmer_rc), &min_lmer, &min_lmer_freq, &min_pos);
                num_recalculate_lmer++;
            }

            if (lmer_out_freq < min_lmer_freq) {
                min_lmer = lmer_out;
                min_lmer_freq = lmer_out_freq;
                min_pos = KMER_LENGTH-LMER_LENGTH;
            }
            r++;
            min_pos--;
            if (counter == 0) {
                rc_idx = num_kmers_added-rd_off-1;

                RevKmerPair this_pair = (!kmers_in_read)? RevKmerPair{kmer_rc, min_lmer, true}: RevKmerPair{kmer_rc, min_lmer, false}; ;
#ifdef CHECK_BIDIRECT
                check_kmer_revc(extract_kmer(kmer_counter[rc_idx].fkmer), extract_kmer(kmer_rc));
#endif      
                RevKmerPair krd_pair = (extract_kmer(kmer_counter[rc_idx].fkmer)<extract_kmer(kmer_rc)) ? kmer_counter[rc_idx] : this_pair;
                if (krd_pair.is_blank) 
                    partial_blank_kmer_counts[retrieve_proc_id(krd_pair.min_lmer)].push_back(KmerPair{krd_pair.fkmer, read_id});
                else
                    partial_kmer_counts[retrieve_proc_id(krd_pair.min_lmer)].push_back(KmerPair{krd_pair.fkmer, read_id});

                rd_off++;
            }
            kmers_in_read++;
        }

#endif

        if (num_kmers > MAX_KMER_COUNT){
            double t25 = MPI_Wtime ();
            // initiate collective communication to pass k-mers and their respective counts to rightful owners
            // calculate global owner of each k-mer and populate the k-mer and count to their respective 'p'th local buffers
            // if global position of a k-mer in my k_map != me, delete k-mer from my k_map
            // iterate through k-mers recieved after collective communication ends, and add k-mers to my k_map 
            // reset num_kmers count to 0.

            double T5 = MPI_Wtime();

            for (int t=0; t<size; t++)
            {
                kpos=0;
                sort(partial_blank_kmer_counts[t].begin(), partial_blank_kmer_counts[t].end());

                int nblank = SortAndAggregatePair(partial_blank_kmer_counts[t], kmer_send_buf[t]);
                kpos += nblank;
                num_blanks[t] = nblank;
                partial_blank_kmer_counts[t].clear();
                partial_blank_kmer_counts[t].shrink_to_fit();

                sort(partial_kmer_counts[t].begin(), partial_kmer_counts[t].end());

                kpos += SortAndAggregatePair(partial_kmer_counts[t], kmer_send_buf[t]);
                partial_kmer_counts[t].clear();
                partial_kmer_counts[t].shrink_to_fit();
                scounts_kmer[t] = kpos;
            }


            double T6 = MPI_Wtime();
            tmap_insert_time += (T6-T5);

            transfer_kmers (scounts_kmer, kmer_send_buf, num_blanks);
            num_kmers = 0;

            scounts_kmer.clear();
            num_blanks.clear();
            double t26 = MPI_Wtime ();
            sl_win_time_int += (t26-t25);
        } // end of if condition

        assert(cur_read_length < KMER_LENGTH || kmers_per_read == (cur_read_length-KMER_LENGTH+1));
        assert(read_len == cur_read_length);
        assert(tmp_read.size() == ceil((double)((double)cur_read_length/nels_per_value_w)));

        size_t vertex_data_start_pos = vtxDataPtrs[idx];
        if ( ispaired && !isForward ) {
            size_t fwd_read_length = reads_data[vertex_data_start_pos];
            size_t fwd_vertex_data_size = ceil((double)((double)fwd_read_length/nels_per_value_w));
            vertex_data_start_pos += fwd_vertex_data_size + 1;
        }
        if ( ispaired && isForward )
            assert(vtxDataPtrs[idx+1] > vertex_data_start_pos + tmp_read.size() + 1);
        else
            assert(vtxDataPtrs[idx+1] == vertex_data_start_pos + tmp_read.size() + 1);
        reads_data[vertex_data_start_pos] = cur_read_length;
        std:copy(tmp_read.begin(), tmp_read.end(), reads_data.begin()+vertex_data_start_pos+1);

        read_id++;
        num_reads++;

        #ifdef LMER_DEBUG
        fprintf(fp_d, "----------------------------\n");
        #endif
    } //end of processing all reads for a given process

    #ifdef LMER_DEBUG
    fclose (fp_d);
    #endif

    *n_kmers = num_kmers;
    scounts_kmer.shrink_to_fit();

    return num_reads;

}

void process_remaining_kmers(
                     std::vector<std::vector<KmerPair>> &partial_kmer_counts,
                     std::vector<std::vector<KmerPair>> &partial_blank_kmer_counts) 
{

     std::vector<int> scounts_kmer (size,0);
     std::vector<int> num_blanks (size, 0);
     std::vector< std::vector<KmerBucket> > kmer_send_buf(size);
     int kpos=0;

      double T7 = MPI_Wtime();
          for (int t=0; t<size; t++)
          {
             kpos=0;
             if(partial_blank_kmer_counts[t].size())
             {     
                  sort(partial_blank_kmer_counts[t].begin(), partial_blank_kmer_counts[t].end());
                  int nblank = SortAndAggregatePair(partial_blank_kmer_counts[t], kmer_send_buf[t]);
                  kpos += nblank;
                  num_blanks[t] = nblank;
                  partial_blank_kmer_counts[t].clear();
                  partial_blank_kmer_counts[t].shrink_to_fit();
             }
             if(partial_kmer_counts[t].size())
             {
                 sort(partial_kmer_counts[t].begin(), partial_kmer_counts[t].end());
                kpos += SortAndAggregatePair(partial_kmer_counts[t], kmer_send_buf[t]);
                partial_kmer_counts[t].clear();
                partial_kmer_counts[t].shrink_to_fit();
            }
            scounts_kmer[t] = kpos;
          }


     double T8 = MPI_Wtime();
     tmap_insert_time += (T8-T7);

    //clear the partial counts
    for(int k=0; k<size; k++)
        partial_kmer_counts[k].clear();         
 
    if (rank==0)
        printf("transfer kmers last round\n");

    transfer_kmers (scounts_kmer, kmer_send_buf, num_blanks);
         
    scounts_kmer.clear();
    scounts_kmer.shrink_to_fit();
    num_blanks.clear();
    num_blanks.shrink_to_fit();

}

void print_kmer_count_timers()
{

    MPI_Reduce(&alltoall_time, &global_alltoall_time, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if (rank == 0) printf ("Average time for Alltoall across all procs (secs): %f \n", 
                            (double)global_alltoall_time/(double)size);

    MPI_Reduce(&pack_sbuf_time, &global_pack_sbuf_time, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if (rank == 0) printf ("Average time for pack_sbuf_time across all procs (secs): %f \n", 
                            (double)global_pack_sbuf_time/(double)size);

    MPI_Reduce(&alltoallv_time, &global_alltoallv_time, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if (rank == 0) printf ("Average time for AlltoallV across all procs (secs): %f \n", 
                            (double)global_alltoallv_time/(double)size);

    MPI_Reduce(&unpack_rbuf_time, &global_unpack_rbuf_time, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if (rank == 0) printf ("Average time for unpack_rbuf_time across all procs (secs): %f \n", 
                            (double)global_unpack_rbuf_time/(double)size);

    MPI_Reduce(&unpack_rbuf_sort, &global_unpack_rbuf_sort, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if (rank == 0) printf ("Average time for unpack_rbuf_time:sort across all procs (secs): %f \n", 
                            (double)global_unpack_rbuf_sort/(double)size);

    MPI_Reduce(&unpack_rbuf_insert, &global_unpack_rbuf_insert, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if (rank == 0) printf ("Average time for unpack_rbuf_time:insert across all procs (secs): %f \n", 
                            (double)global_unpack_rbuf_insert/(double)size);

    MPI_Reduce(&unpack_rbuf_acc, &global_unpack_rbuf_acc, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if (rank == 0) printf ("Average time for unpack_rbuf_time:acc across all procs (secs): %f \n", 
                            (double)global_unpack_rbuf_acc/(double)size);

    MPI_Reduce(&sl_win_time_int, &global_sl_win_time_int, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if (rank == 0) printf ("Average time for partial Sliding Window across all procs (secs): %f \n", 
                            (double)global_sl_win_time_int/(double)size);

    MPI_Reduce(&sl_win_time, &global_sl_win_time, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if (rank == 0) printf ("Average time for Sliding Window across all procs (secs): %f \n", 
                            (double)global_sl_win_time/(double)size);

    MPI_Reduce(&sl_lmer_freq, &global_sl_lmer_freq, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if (rank == 0) printf ("Average time for lmer Sliding Window across all procs (secs): %f \n", 
                            (double)global_sl_lmer_freq/(double)size);

    MPI_Reduce(&vec_insert_time, &global_vec_insert_time, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if (rank == 0) printf ("Average time for vector insert across all procs (secs): %f \n", 
                            (double)global_vec_insert_time/(double)size);

    MPI_Reduce(&tmap_insert_time, &global_tmap_insert_time, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if (rank == 0) printf ("Average time for temp_map insert across all procs (secs): %f \n", 
                            (double)global_tmap_insert_time/(double)size);

    MPI_Reduce(&num_recalculate_lmer, &global_num_recalculate_lmer, 1, MPI_UINT64_T, MPI_SUM, 0, MPI_COMM_WORLD);
    if (rank == 0) printf ("Average number of re-calculations of min l-mer across all procs: %f \n",
                            (double)global_num_recalculate_lmer/(double)size);

    double global_kmer_count_time_max;
    MPI_Reduce(&kmer_count_time, &global_kmer_count_time, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&kmer_count_time, &global_kmer_count_time_max, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    if (rank == 0) {
        printf ("Average time for performing k-mer bucketing across all procs (secs): %f \n",
                            (double)global_kmer_count_time/(double)size);
        printf ("Maximum time for performing k-mer bucketing across all procs (secs): %f \n",
                            (double)global_kmer_count_time_max);

        fprintf (stderr, "Average time for performing k-mer bucketing across all procs (secs): %f \n",
                            (double)global_kmer_count_time/(double)size);
        fprintf (stderr, "Maximum time for performing k-mer bucketing across all procs (secs): %f \n",
                            (double)global_kmer_count_time_max);
    }

}

void print_prof_size()
{
    uint64_t all_proc_kmer_count = 0;
    uint64_t tmp_kmer_count = kmer_proc_buf.size();

    MPI_Reduce(&tmp_kmer_count, &all_proc_kmer_count, 1, MPI_UINT64_T, MPI_SUM, 0, MPI_COMM_WORLD);

    if (rank == 0) { 
        printf("Total distinct (k+1)-mer entries across all proc's: %lu\n", all_proc_kmer_count); 
        printf("Number of batch iterations: %d \n", num_batch_transfers);
    }
}

void free_proc_buffers()
{
    kmer_proc_buf.clear();
    kmer_proc_buf.shrink_to_fit();
}

void free_kmer_count_buffers()
{
   
    col_ind_readid.clear();
    col_ind_readid.shrink_to_fit();

    bucket_list.clear();
    bucket_list.shrink_to_fit();

    row_bucket_ptr.clear();
    row_bucket_ptr.shrink_to_fit();

}

void print_kmer_buckets (int klen)
{

    std::string output_filename("kmers_p" + std::to_string(rank) + "_k" + std::to_string(klen) + ".csv");

    FILE *f = fopen(output_filename.c_str(), "w");
    if (f == NULL)
    {
        printf("Error opening kmers file!\n");
        exit(1);
    }

    for (size_t it=0; it<kmer_proc_buf.size(); it++)
    {
      
        kmer_t full_kmer = kmer_proc_buf[it].seq;
        kmer_t pred_k=extract_pred(full_kmer);
        kmer_t succ_k=extract_kmer(full_kmer);
        int pbp=(int)kmer_proc_buf[it].pred_bp;

        char kmer_out[klen+1];
        for (size_t l=0; l<KMER_LENGTH; l++) {
              kmer_out[l] = el_to_char(kmerel (succ_k, l));
         }
         kmer_out[KMER_LENGTH] = '\0';

        fprintf(f, "%lu, %lu, %lu, %lu, %s, %lu, (", full_kmer, pred_k, pbp, succ_k, kmer_out, kmer_proc_buf[it].ReadIds.size());
        for (size_t j=0; j<kmer_proc_buf[it].ReadIds.size(); j++)
             fprintf(f, "%lu ", kmer_proc_buf[it].ReadIds[j]);
        fprintf(f, ")\n");

    
    }

    fclose(f);
}

#ifdef DEBUG_PRINT_BUCKET
//Printing buffer contents for debugging
void print_kmer_buckets_pair (int klen,
                              std::vector < std::vector<KmerBucket> > &kmer_send_buf)
{

    std::string output_filename("kmerpair_p" + std::to_string(rank) + "_k" + std::to_string(klen) + ".csv");

    FILE *f = fopen(output_filename.c_str(), "w");
    if (f == NULL)
    {
        printf("Error opening kmers file!\n");
        exit(1);
    }

    for (int k=0; k<size; k++) 
    {
    for (size_t it=0; it<kmer_send_buf[k].size(); it++)
     {
      
        kmer_t full_kmer = kmer_send_buf[k][it].seq;
        kmer_t pred_k=extract_pred(full_kmer);
        kmer_t succ_k=extract_kmer(full_kmer);

        char kmer_out[klen+1];
        for (size_t l=0; l<KMER_LENGTH; l++) {
              kmer_out[l] = el_to_char(kmerel (succ_k, l));
         }
         kmer_out[KMER_LENGTH] = '\0';

        fprintf(f, "%lu, %lu, %lu, %s, (", full_kmer, pred_k, succ_k, kmer_out);
        for (size_t j=0; j<kmer_send_buf[k][it].ReadIds.size(); j++)
             fprintf(f, "%lu ", kmer_send_buf[k][it].ReadIds[j]);
        fprintf(f, ")\n");

    
     }
    }

    fclose(f);
}

void print_kmer_buckets_rank (int klen,
                              std::vector <KmerPair> &kmer_send_buf)
{

    std::string output_filename("kmerpair1d_p" + std::to_string(rank) + "_k" + std::to_string(klen) + ".csv");

    FILE *f = fopen(output_filename.c_str(), "w");
    if (f == NULL)
    {
        printf("Error opening kmers file!\n");
        exit(1);
    }

    for (size_t it=0; it<kmer_send_buf.size(); it++)
     {
      
        kmer_t full_kmer = kmer_send_buf[it].bkmer;
        kmer_t pred_k=extract_pred(full_kmer);
        kmer_t succ_k=extract_kmer(full_kmer);

        char kmer_out[klen+1];
        for (size_t l=0; l<KMER_LENGTH; l++) {
              kmer_out[l] = el_to_char(kmerel (succ_k, l));
         }
         kmer_out[KMER_LENGTH] = '\0';

        fprintf(f, "%lu, %lu, %lu, %s, %lu\n", full_kmer, pred_k, succ_k, kmer_out, kmer_send_buf[it].read_id);
    
     }

    fclose(f);
}

void print_kmer_buckets_list (int klen)
{

    std::string output_filename("bucketlist_p" + std::to_string(rank) + "_k" + std::to_string(klen) + ".csv");

    FILE *f = fopen(output_filename.c_str(), "w");
    if (f == NULL)
    {
        printf("Error opening kmers file!\n");
        exit(1);
    }

    for (size_t it=0; it<ReadBucketList.size(); it++)
     {
      
        kmer_t key = ReadBucketList[it].seq;

        char kmer_out[klen+1];
        for (size_t l=0; l<KMER_LENGTH; l++) {
              kmer_out[l] = el_to_char(kmerel (key, l));
         }
         kmer_out[KMER_LENGTH] = '\0';

        fprintf(f, "%lu, %s, (", key, kmer_out);
        for (size_t j=0; j<ReadBucketList[it].ReadIds.size(); j++)
             fprintf(f, "%lu ", ReadBucketList[it].ReadIds[j]);
        fprintf(f, ")\n");

    
     }

    fclose(f);
}

void print_csr_bucket_data ()
{
   std::string output_filename("bucketlist_p" + std::to_string(rank) + ".csv");

    FILE *f = fopen(output_filename.c_str(), "w");
    if (f == NULL)
    {
        printf("Error opening kmers file!\n");
        exit(1);
    }

    int offset_start=0, offset_end=0;
    for (size_t it=0; it<bucket_list.size(); it++)
    {
        kmer_t key = bucket_list[it];
        
        char kmer_out[KMER_LENGTH+1];
        for (size_t l=0; l<KMER_LENGTH; l++) {
              kmer_out[l] = el_to_char(kmerel (key, l));
         }
         kmer_out[KMER_LENGTH] = '\0';

        fprintf(f, "%lu, %s, ", key, kmer_out);
        offset_start = row_bucket_ptr[it];
        offset_end = row_bucket_ptr[it+1];
        fprintf(f, "%lu, (", (offset_end-offset_start));
        for (size_t j=offset_start; j<offset_end; j++)
            fprintf(f, "%lu ", col_ind_readid[j]);
        fprintf(f, ")\n");

    } 

    fclose(f);

}

#endif

void print_hgraph_bucket_data (HGRAPH_DATA &hg)
{
   std::string output_filename(outputDirName + "/bucketlist_p" + std::to_string(rank) + ".csv");

    FILE *f = fopen(output_filename.c_str(), "w");
    if (f == NULL)
    {
        printf("Error opening kmers file!\n");
        exit(1);
    }


    int offset_start=0, offset_end=0;
    for (int it=0; it<hg.getNumMyHEdges(); it++)
    {
        ZOLTAN_ID_TYPE key = hg.edgeGID[it];
        
        char kmer_out[KMER_LENGTH+1];
        for (size_t l=0; l<KMER_LENGTH; l++) {
              kmer_out[l] = el_to_char(kmerel (key, l));
         }
         kmer_out[KMER_LENGTH] = '\0';

        fprintf(f, "%lu, %s, ", key, kmer_out);
        offset_start = hg.nborIndex[it];
        offset_end = hg.nborIndex[it+1];
        fprintf(f, "%lu, (", (offset_end-offset_start));
        for (size_t j=offset_start; j<offset_end; j++)
            fprintf(f, "%lu ", hg.nborGID[j]);
        fprintf(f, ")\n");

    } 

    fclose(f);

}

int nonce=0;
void filesystem_barrier()
{
    nonce++;
    system("sync");
    std::filesystem::create_directory(outputDirName + "/" + std::to_string(nonce) + "-" + std::to_string(rank));
    system("sync");
    for (int i = 0; i < size; ++i) {
        while (!std::filesystem::is_directory(outputDirName + "/" + std::to_string(nonce) + "-" + std::to_string(i))) {
            std::this_thread::sleep_for(std::chrono::milliseconds(10));
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    std::filesystem::remove(outputDirName + "/" + std::to_string(nonce) + "-" + std::to_string(rank));
}

std::string read_tVecToStr(std::vector<read_t> & vtxData, size_t start, size_t end)
{
    assert(start < end);
    assert(end <= vtxData.size());
    size_t this_read_length = vtxData[start]; // in the begining of each read in read_t, we should encode the size in bp
    assert(this_read_length <= (end - start - 1) * nels_per_value_w);
    std::string out_read;
    for(int j=0; j<this_read_length; j++) {
        out_read += el_to_char(get_bp(j, vtxData, start+1, this_read_length));
    }
    return out_read;
}

std::vector<std::string> read_tVecToStr(std::vector<read_t> & vtxData, size_t start, size_t end, bool ispaired)
{
    std::vector<std::string> read;
    read.push_back(read_tVecToStr(vtxData, start, end));
    if (ispaired) {
        size_t fwd_read_length = vtxData[start];
        assert(read[0].size() == fwd_read_length);
        start = start + ceil((double)((double)fwd_read_length/nels_per_value_w)) + 1;
        read.push_back(read_tVecToStr(vtxData, start, end));
    }
    return read;
}

void print_read_data(size_t numMyVertices, std::vector<size_t> & partsPtr, std::vector<int> & myParts, std::vector<size_t> & vtxDataPtrs, std::vector<read_t> & vtxData)
{
    size_t numMyParts = partsPtr.size() - 1;
    uint64_t local_read_count=numMyVertices, global_read_count=0;
    uint64_t local_parts_count=0, global_parts_count=0;
    MPI_Reduce(&local_read_count, &global_read_count, 1, MPI_UINT64_T, MPI_SUM, 0, MPI_COMM_WORLD);
    if (rank == 0) {
        printf("Total number of reads partitioned across %d procs: %lu\n", size, global_read_count);
        assert(global_read_count == global_nvertices);
    }

    assert(numMyParts == myParts.size());
    for (size_t partID=0; partID<myParts.size(); partID++) {
        if (myParts[partID] > local_parts_count) local_parts_count = myParts[partID];
    }
    MPI_Allreduce(&local_parts_count, &global_parts_count, 1, MPI_UINT64_T, MPI_MAX, MPI_COMM_WORLD);
    global_parts_count += 1;
    if (rank == 0) {
        printf("Total number of parts across %d procs: %lu\n", size, global_parts_count);
    }

    for (size_t partID = 0; partID < numMyParts; partID++) {
        std::string output_filename(outputDirName + "/reads_r" + std::to_string(rank) + "_p" + std::to_string(myParts[partID]) + ".fa");
        if (size > global_parts_count) output_filename += ".temp";

        FILE *f = fopen(output_filename.c_str(), "w");

        size_t partStart, partEnd;
        partStart = partsPtr[partID];
        partEnd = partsPtr[partID+1];
        for (int i=partStart; i<partEnd; i++)
        {
            std::vector<std::string> out_reads = read_tVecToStr(vtxData, vtxDataPtrs[i], vtxDataPtrs[i+1], ispaired);
            if (!ispaired) {
                assert(out_reads.size() == 1);
                std::string output_read_name(">read_" + std::to_string(i-partStart) + "_p_" + std::to_string(rank) + ".fa");
                fprintf(f, "%s\n", output_read_name.c_str());
                fprintf(f, "%s\n", out_reads[0].c_str());
            } else {
                assert(out_reads.size() == 2);
                std::string output_read_name1(">read_" + std::to_string(i-partStart) + "_one" "_p_" + std::to_string(rank) + ".fa");
                std::string output_read_name2(">read_" + std::to_string(i-partStart) + "_two" "_p_" + std::to_string(rank) + ".fa");
                fprintf(f, "%s\n", output_read_name1.c_str());
                fprintf(f, "%s\n", out_reads[0].c_str());
                fprintf(f, "%s\n", output_read_name2.c_str());
                fprintf(f, "%s\n", out_reads[1].c_str());
            }
        }

        fclose(f);
    }

    MPI_Barrier(MPI_COMM_WORLD);

    // When the number of MPI rank is larger than the number of
    // partitions, parts will be splitted among MPI ranks.
    // So we will write to a temporary location (*.temp) for parallel writing
    // of MPI ranks first, then merge the *temp files that belong to the
    // same part.
    if (size > global_parts_count) {
        filesystem_barrier();

        assert(numMyParts <= 1);

        int64_t partIDs[size], all_partIDs[size];
        for (size_t i=0; i<size; i++) partIDs[i] =-1;
        if (numMyParts > 0)
            partIDs[rank] = myParts[0];
        MPI_Allreduce(&partIDs, &all_partIDs, size, MPI_INT64_T, MPI_MAX, MPI_COMM_WORLD);
        if (rank < global_parts_count) {
            // merge if necessary
            int partID = rank;
            std::vector<int> ranksToMerge;
            for (int i=0; i<size; i++) if (all_partIDs[i] == partID) ranksToMerge.push_back(i);

            assert(ranksToMerge.size() >= 1);
            std::string outputfilename(outputDirName + "/reads_r" + std::to_string(partID) + "_p" + std::to_string(partID) + ".fa");

            try {
                std::filesystem::rename(outputDirName + "/reads_r" + std::to_string(ranksToMerge[0]) + "_p" + std::to_string(partID) + ".fa.temp", outputfilename);
            } catch (std::filesystem::filesystem_error& e) {
                std::cout << e.what() << '\n';
            }

            for (size_t i=1; i<ranksToMerge.size(); i++) {
                std::string catstr("cat " + outputDirName + "/reads_r" + std::to_string(ranksToMerge[i]) + "_p" + std::to_string(partID) + ".fa.temp >> " + outputfilename);
                system(catstr.c_str());
                std::filesystem::remove(outputDirName + "/reads_r" + std::to_string(ranksToMerge[i]) + "_p" + std::to_string(partID) + ".fa.temp" );
            }
        }
    }
}

void print_hgraph_read_data (HGRAPH_DATA &hg)
{
    print_read_data(hg.numMyVertices, hg.partPtr, hg.myParts, hg.vtxDataPtrs, hg.vtxData);
}

void print_graph_read_data (GRAPH_DATA &g)
{
    print_read_data(g.numMyVertices, g.partPtr, g.myParts, g.vtxDataPtrs, g.vtxData);
}



std::pair<uint64_t,uint64_t> generate_input_to_Zoltan()
{
    //Generate into HyperGraph: List of buckets and corresponding read id's
    std::vector<int> hist_bucket(max_reads_per_bucket,0);
    std::vector<int> global_hist_bucket(max_reads_per_bucket,0);
    uint64_t num_outliers=0, num_buckets;
    GenerateSortedCSRListPerBucket(kmer_proc_buf, col_ind_readid, bucket_list, row_bucket_ptr, hist_bucket, &num_outliers, &num_buckets);


#ifdef PRINT_BUCKET_HISTOGRAM

    MPI_Reduce (hist_bucket.data(), global_hist_bucket.data(), max_reads_per_bucket, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);

    if (rank == 0) {
        std::string histogram=outputDirName + "/reads_histogram.txt";
        FILE *fr = fopen(histogram.c_str(), "w");
        if (fr == NULL)
        {
            printf("Error opening kmers file!\n");
            exit(1);
        }
        
        for (int i=0; i<(int)global_hist_bucket.size(); i++)
             fprintf(fr,"i: %d, %d\n", i, global_hist_bucket[i]);
        fclose(fr);
    }

    hist_bucket.clear();
    global_hist_bucket.clear();
#endif

    print_prof_size();
    uint64_t global_distinct_kmer_count = 0;
    uint64_t local_kmer_count = bucket_list.size();

    MPI_Reduce(&local_kmer_count, &global_distinct_kmer_count, 1, MPI_UINT64_T, MPI_SUM, 0, MPI_COMM_WORLD);

    if (rank == 0) { 
        printf("Total distinct k-mer bucket entries across all proc's: %lu\n", global_distinct_kmer_count);
        fprintf (stderr, "Total distinct k-mer bucket entries across all proc's: %lu\n", global_distinct_kmer_count); 
    }

    // free the local proc buffers storing the k+1-mer's
    free_proc_buffers();

    return std::make_pair(num_outliers,num_buckets);
    
}

std::pair<uint64_t,uint64_t> generate_input_to_ParMetis(uint64_t reads_disp, uint64_t num_reads_per_p)
{

    //Generate input to ParMetis
    std::vector < std::vector< std::pair<read_t, read_t>> > local_edge_list;
    uint64_t num_outliers=0, num_buckets=0;

    local_edge_list = GenerateEdgeListPerProcess(kmer_proc_buf, &num_outliers, &num_buckets);
    print_prof_size();

    // free the local proc buffers storing the k+1-mer's
    free_proc_buffers();

    transfer_edge_list(local_edge_list, reads_disp, num_reads_per_p);
    
    return std::make_pair(num_outliers,num_buckets);

}

size_t process_sliding_window(input_read_data &rdata, uint64_t reads_disp, size_t read_start_idx, bool isForward)
{
    std::vector< std::vector<KmerPair> > partial_kmer_counts(size);
    std::vector< std::vector<KmerPair> > partial_blank_kmer_counts(size);
    int num_kmers = 0;

    if (rank==0)
        fprintf(stderr, "Starting sliding window\n");
        

    size_t num_treads = Sliding_window_variableReadLength (rdata, &num_kmers, partial_kmer_counts, partial_blank_kmer_counts, reads_disp, isForward);

    #ifdef DEUBG
    fprintf(stderr, "Completed going through the reads data\n");
    #endif

    // initiate communication for the residual kmers
    //if (num_kmers) { // this should be avoided as there could be a rank
                       // having zero num_kmers and inside of the process_remaining_kmers,
                       // there are communications among all ranks
                       // remains are not necessarily the same among all ranks

          process_remaining_kmers (partial_kmer_counts, partial_blank_kmer_counts);
    //}

    if (rank==0)
        fprintf(stderr, "Completed sliding window\n");

    for (int i=0; i<size; i++)
         partial_kmer_counts[i].clear();
    partial_kmer_counts.shrink_to_fit();

    return num_treads;

}

size_t getCompressedSize(input_read_data &rdata_fwd, input_read_data &rdata_rev) {
    // assign vtxDataPtrs as well
    if (ispaired) {
        assert(rdata_fwd.num_of_reads == rdata_rev.num_of_reads);
    } else {
        assert(rdata_rev.num_of_reads == 0);
    }
    size_t num_of_vertices = rdata_fwd.num_of_reads;

    vtxDataPtrs.resize(num_of_vertices+1);
    vtxDataPtrs[0] = 0;
    for (size_t i=0; i<num_of_vertices; i++) {
        vtxDataPtrs[i+1]=vtxDataPtrs[i];

        size_t fwd_read_length = rdata_fwd.read_data_ptr[i+1]-rdata_fwd.read_data_ptr[i];
        vtxDataPtrs[i+1]+=ceil((double)((double)fwd_read_length/nels_per_value_w)) + 1;

        if (ispaired) {
            size_t rev_read_length = rdata_rev.read_data_ptr[i+1]-rdata_rev.read_data_ptr[i];
            vtxDataPtrs[i+1]+=ceil((double)((double)rev_read_length/nels_per_value_w)) + 1;
        }
    }
    return vtxDataPtrs[num_of_vertices];
}

void perform_kmer_counting(input_read_data &rdata_fwd, input_read_data &rdata_rev, uint64_t reads_disp, uint64_t local_nlines)
{

    double start_t = MPI_Wtime ();

    //pre-allocate buffer to hold all the reads from both the frwd and rev read files
    size_t compressed_read_size = getCompressedSize(rdata_fwd, rdata_rev);
    reads_data.resize(compressed_read_size, 0);
    fprintf(stderr,"rank:%d, reads_data_size: %lu, local_nlines: %lu, reads_disp: %lu\n", rank, reads_data.size(), local_nlines, reads_disp);

    #ifdef DEBUG
    assert(vtxDataPtrs[rdata_fwd.num_of_reads] == reads_data.size());
    for (size_t idx=0; idx<vtxDataPtrs.size(); idx++)
        fprintf(stderr,"The %luth vtxDataPtrs : %lu\n", idx, vtxDataPtrs[idx]);
    #endif

    size_t read_start_idx=0;

    double t23 = MPI_Wtime ();

    //processing the fwd reads file
    size_t num_freads=process_sliding_window(rdata_fwd, reads_disp, read_start_idx, true);
    #ifdef DEBUG
    fprintf(stderr,"Finish processing sliding window forward\n");
    MPI_Barrier(MPI_COMM_WORLD);
    #endif

    //processing the rev reads file
    size_t num_rreads=0;
    if (ispaired) num_rreads = process_sliding_window(rdata_rev, reads_disp, read_start_idx+1, false);

    #ifdef DEBUG
    MPI_Barrier(MPI_COMM_WORLD);
    for (size_t idx=1; idx<vtxDataPtrs.size(); idx++) {
        fprintf(stderr,"The %luth (paired) read: \n", idx-1);
        std::vector<std::string> r = read_tVecToStr(reads_data, vtxDataPtrs[idx-1], vtxDataPtrs[idx], ispaired);
        for (size_t j=0; j<r.size(); j++)
            fprintf(stderr,"%s\n", r[j].c_str());
    }
    MPI_Barrier(MPI_COMM_WORLD);
    #endif

    double t24 = MPI_Wtime ();
    sl_win_time = t24 - t23;

    if (ispaired)
        assert((num_freads+num_rreads)==local_nlines);
    else
        assert(num_freads==local_nlines);

    lmer_frequency.clear();
    global_lmer_frequency.clear();

    //free the reads
    free(rdata_fwd.read_data);
    if (ispaired)
        free(rdata_rev.read_data);

    std::pair<uint64_t,uint64_t> binfo(0,0);
    switch (ptype) {
    case Zoltan:
         if (rank == 0) printf("Generating input for Zoltan\n");
         binfo=generate_input_to_Zoltan();
         break;
    case ParMetis:
         if (rank == 0) printf("Generating input for ParMetis\n");
         binfo=generate_input_to_ParMetis(reads_disp, num_freads);
         break;
    default:
         assert(0 && "Should not reach here!! Please provide a correct option for the type of partitioner.");
         break;
    }

    //calculate number of buckets discarded
    uint64_t num_outliers=binfo.first, global_num_outliers=0;
    MPI_Reduce (&num_outliers, &global_num_outliers, 1, MPI_UINT64_T, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    if (rank == 0) fprintf(stderr,"Number of outliers (buckets discarded): %lu\n", global_num_outliers);
    if (rank == 0) printf("Number of outliers (buckets discarded): %lu\n", global_num_outliers);

    //calculate total number of buckets
    uint64_t num_buckets=binfo.second, global_num_buckets=0;
    MPI_Reduce (&num_buckets, &global_num_buckets, 1, MPI_UINT64_T, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    if (rank == 0) fprintf(stderr,"Number of buckets: %lu\n", global_num_buckets);
    if (rank == 0) printf("Number of buckets: %lu\n", global_num_buckets);
    if (rank == 0) fprintf(stderr,"percentage of buckets discarded: %f\n", 
                                   (double)((double)global_num_outliers/(double)global_num_buckets)*100);
    if (rank == 0) printf("percentage of buckets discarded: %f\n", 
                                   (double)((double)global_num_outliers/(double)global_num_buckets)*100);

    double end_t = MPI_Wtime ();

    kmer_count_time = end_t - start_t;

    MPI_Barrier(MPI_COMM_WORLD);

    if (rank==0)
        fprintf(stderr, "Completed kmer counting\n");

    print_kmer_count_timers();
}
