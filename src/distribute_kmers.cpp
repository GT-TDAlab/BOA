#include <mpi.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <assert.h>
#include <unistd.h>
#include <math.h>
#include <inttypes.h>
#include <vector>
#include <algorithm>
#include <unordered_map>
#include <parallel/algorithm>
#include <numeric>
#include <omp.h>
#include "distribute_kmers.h"
#include "timers.h"
#include "zoltan_partition.h"

#include <fstream>
#include <filesystem>

long int MAX_KMER_COUNT=0;
int rank, size;
int coverage=0;
std::string inputFileName;
std::string inputFileNameMate;
bool ispaired=true;
int read_length=0;
int num_buckets=0;
int num_threads=0;
int num_partitions=0;
unsigned int partitioning_seed=0;
std::string outputDirName;
PartitionerType ptype;

int num_batch_transfers=0;
uint64_t global_nlines=0;
uint64_t global_nvertices=0;

std::vector<lmer_t> lmer_frequency(LMER_SIZE,0);
std::vector<lmer_t> global_lmer_frequency(LMER_SIZE,0);

// Buffer holding the reads displacement for each process
std::vector<uint64_t> global_reads_disp;

// Buffer containing list of read id's for each k-mer bucket
std::vector<read_t> col_ind_readid;
std::vector<kmer_t> bucket_list;
std::vector<int> row_bucket_ptr;

// Buffer holding the reads data
std::vector<size_t> vtxDataPtrs;
std::vector<read_t> reads_data_t;
std::vector<read_t> reads_data;
int read_size=0;

uint64_t hasha = 68111;
uint64_t hashb = 105929;

void parseCommandLine(const int argc, char * const argv[]);

int retrieve_read_proc_owner(read_t read_id)
{
    int proc=-1;
    std::vector<read_t>::iterator low;

    low=std::lower_bound (global_reads_disp.begin(), global_reads_disp.end(), read_id);
    int dist=low-global_reads_disp.begin();
    proc = global_reads_disp[dist] > read_id ? (dist-1) : dist;
    
    if (!(proc<size))
         printf("read_id: %lu, size: %d, proc: %d, rank: %d\n", read_id, size, proc, rank);
    assert(proc < size);
    return proc;
}

// determine which proc holds the l-mer bucket based on a block-cyclic
// distribution of buckets across the processes
int retrieve_proc_id (lmer_t min_lmer)
{
    int proc_num = (int)((int)(hash31(hasha, hashb, min_lmer)) % (int)size); // randomized

    return proc_num;
}
#ifdef EXTEND_KMER
int retrieve_proc_id (kmer_t min_lmer)
{
    int proc_num = (int)((int)(uhash31(hasha, hashb, min_lmer)) % (int)size); // randomized

    return proc_num;
}
#endif

void set_num_threads()
{

#pragma omp parallel
{
    num_threads = omp_get_num_threads();
}

}

void intiate_zoltan_processing_hg(int argc, char **argv,
                               uint64_t reads_disp, uint64_t num_reads_per_p)
{
    // construct the HGRAPH_DATA data-structure
    HGRAPH_DATA hg = populate_hgraph_data(reads_disp, num_reads_per_p);
    perform_zoltan_partition_hg (argc, argv, hg, num_partitions);

    // output the read files to disk
    print_hgraph_read_data(hg);
}

void intiate_zoltan_processing_g(int argc, char **argv,
                               uint64_t reads_disp, uint64_t num_reads_per_p)
{
    //construct the GRAPH_DATA data-structure
    GRAPH_DATA g = populate_graph_data(reads_disp, num_reads_per_p);
    perform_zoltan_partition_g (argc, argv, g, num_partitions);

    // output the read files to disk
    print_graph_read_data(g);
}

int main(int argc, char **argv) {

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    parseCommandLine(argc, argv);
    set_num_threads();


    if (rank == 0) {
        std::string pathname(outputDirName);
        std::filesystem::create_directories(pathname);
    }
    MPI_Barrier(MPI_COMM_WORLD);


    uint64_t reads_disp_fwd=0, reads_disp_rev=0; 
    uint64_t num_reads_per_p_fwd=0, num_reads_per_p_rev=0; // number of reads within the rank
    uint64_t fwd_nlines=0, rev_nlines=0; // total number of reads across all ranks
        
    input_read_data rdata_fwd = perform_input_reading(rank, size, inputFileName, &num_reads_per_p_fwd, &fwd_nlines);
    input_read_data rdata_rev = perform_input_reading(rank, size, inputFileNameMate, &num_reads_per_p_rev, &rev_nlines);

    if (ispaired) {
        // Since we are reading read-pairs, the number of reads assigned to each process and the
        // total number of reads across all procs should be the same.
        assert(fwd_nlines==rev_nlines);
        // assume the number of reads within each pair is the same and the comments are the same as well, thus the parallel input of the two file is the split the same among the rank
        assert(num_reads_per_p_fwd==num_reads_per_p_rev);
    }

    // if it is not paired, rev_nlines should be equal to zero
    global_nlines=fwd_nlines+rev_nlines;
    // if it is not paired, num_reads_per_p_rev should be equal to zero
    uint64_t local_nlines=num_reads_per_p_fwd+num_reads_per_p_rev;
    global_nvertices=fwd_nlines;

    uint64_t scan_nlines_pair=0;
    // calculate the read displacement for each read-pair per process
    MPI_Scan (&num_reads_per_p_fwd, &scan_nlines_pair, 1, MPI_UINT64_T, MPI_SUM, MPI_COMM_WORLD);
    uint64_t reads_disp_pair=scan_nlines_pair-num_reads_per_p_fwd; // the starting read (pair) id

    // gather all the read displacements across all procs
    global_reads_disp.resize(size,0);
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Allgather(&reads_disp_pair, 1, MPI_UINT64_T, global_reads_disp.data(), 1, MPI_UINT64_T, MPI_COMM_WORLD);
    assert(global_reads_disp[rank]==reads_disp_pair);
    global_reads_disp.push_back(fwd_nlines);
    assert(global_reads_disp.size()==size+1);

    // calculate frequencies of all l-mers in the Read dataset
    double time_l1 = MPI_Wtime ();
    //Sliding_window_l(rdata_fwd.read_data, rdata_fwd.read_data_size);
    //if (ispaired) Sliding_window_l(rdata_rev.read_data, rdata_rev.read_data_size);
    Sliding_window_l_variableReadLength(rdata_fwd);
    if (ispaired) Sliding_window_l_variableReadLength(rdata_rev);

    double time_l2 = MPI_Wtime ();
    sl_lmer_freq = time_l2 - time_l1;

    // Perform Allreduce to obtain global lmer counts across all proc's
    int num_lmers = pow(4, LMER_LENGTH);

    MPI_Allreduce (lmer_frequency.data(), global_lmer_frequency.data(), num_lmers, MPI_UINT64_T, MPI_SUM, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);

    // perform bucketing for forward reads and reverse reads, wherein the k-mers represent buckets
    perform_kmer_counting (rdata_fwd, rdata_rev, reads_disp_pair, local_nlines);

    switch (ptype) {
    case Zoltan:
         if (rank==0)
             fprintf(stderr, "Begin Zoltan partitioning\n");
         intiate_zoltan_processing_hg(argc, argv, reads_disp_pair, num_reads_per_p_fwd);
         break;
    case ParMetis:
         if (rank==0)
             fprintf(stderr, "Begin ParMetis partitioning\n");
         intiate_zoltan_processing_g(argc, argv, reads_disp_pair, num_reads_per_p_fwd);
         break;
    default:
         assert(0 && "Should not reach here!! Please provide a correct option for the type of partitioner.");
         break;
    }

    free_kmer_count_buffers();

    MPI_Barrier(MPI_COMM_WORLD);

    MPI_Finalize();
    return 0;
}

void parseCommandLine(const int argc, char * const argv[])
{
  int ret;
  
  while ((ret = getopt(argc, argv, "f:m:b:r:p:o:t:d:")) != -1) {
    switch (ret) {
    case 'f':
       inputFileName.assign(optarg);
       break;
    case 'm':
       inputFileNameMate.assign(optarg);
       break;
    case 'b':
       MAX_KMER_COUNT = atol(optarg);
       break;
    case 'r':
       read_length = atol(optarg);
       break;
    case 'p':
       num_partitions = atoi(optarg);
       break;
    case 'o':
       outputDirName.assign(optarg);
       break;
    case 't':
       if (std::string(optarg) == "Zoltan")
           ptype = Zoltan;
       else if (std::string(optarg) == "ParMetis")
           ptype = ParMetis;
       else
           assert(0 && "Partition Type input is invalid.");
       break;
    case 'd':
       partitioning_seed = (unsigned int)std::stoul(optarg);
       break;
    default:
       assert(0 && "Should not reach here!!");
       break;
    }
  }

  if (!num_partitions) {
      num_partitions = size;
      if (rank == 0) std::cerr << "Number of partitions not specified with -p, set to default value=number of prcos, which is: " << size << std::endl;
  }

  if (rank==0 && inputFileName.empty()) {
      std::cerr << "Must specify an input FASTA file name with -f" << std::endl;
      MPI_Abort(MPI_COMM_WORLD, -99);
  }

  if (inputFileNameMate.empty()) {
      ispaired=false;
      if (rank==0) std::cout << "Paired-end reads have not been provided with option -m. StringBOA will proceed with the default single-end read processing." << std::endl;
  }

  if (outputDirName.empty()) {
      outputDirName.assign("output");
  }

  if (!MAX_KMER_COUNT) {
      MAX_KMER_COUNT=100000000;
      if (rank==0) std::cout << "batch size not specified with -b, set to default value MAX_KMER_COUNT=100000000" << std::endl;
  }

  if (!ptype) {
      ptype = Zoltan;
      if (rank==0 ) std::cerr << "Proceeding with partitioner Zoltan" << std::endl;
  }

  if (rank == 0) {
           printf("K-mer size: %d, L-mer size: %d, Number of Processes: %d, MAX_KMER_COUNT: %ld, Number of Partitions: %d,",
           WINDW_SIZE+1, LMER_LENGTH, size, MAX_KMER_COUNT, num_partitions);
  }

} // parseCommandLine
