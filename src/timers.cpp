#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <assert.h>

double kmer_count_time=0.0, global_kmer_count_time=0.0;
double alltoall_time=0.0, global_alltoall_time=0.0;
double alltoallv_time=0.0, global_alltoallv_time=0.0;
double partial_contig_time=0.0, global_partial_contig_time=0.0;
double preproc_time=0.0, global_preproc_time=0.0;
double gather_time=0.0, global_gather_time=0.0;
double barrier_time=0.0, global_barrier_time=0.0;
double sort_time=0.0, global_sort_time=0.0;

double contigs_time=0.0, global_contigs_time=0.0;
double io_time=0.0, cleanup_time=0.0, vector_time=0.0;
double count_time=0.0, global_count_time=0.0;
double read_input_time=0.0, global_read_input_time=0.0;
double pack_sbuf_time=0.0, global_pack_sbuf_time=0.0;
double unpack_rbuf_time=0.0, global_unpack_rbuf_time=0.0;
double unpack_rbuf_sort=0.0, global_unpack_rbuf_sort=0.0;
double unpack_rbuf_insert=0.0, global_unpack_rbuf_insert=0.0;
double unpack_rbuf_acc=0.0, global_unpack_rbuf_acc=0.0;
double sl_win_time=0.0, global_sl_win_time=0.0;
double sl_win_time_int=0.0, global_sl_win_time_int=0.0;
double sl_lmer_freq=0.0, global_sl_lmer_freq=0.0;

double vec_insert_time=0.0, global_vec_insert_time=0.0;
double tmap_insert_time=0.0, global_tmap_insert_time=0.0;


//P2 Timers
double p2alltoall_time=0.0, p2global_alltoall_time=0.0;
double p2alltoallv_time=0.0, p2global_alltoallv_time=0.0;
