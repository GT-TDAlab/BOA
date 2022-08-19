#include <stdio.h>
#include <mpi.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <inttypes.h>
#include <vector>
#include <algorithm>
#include <unordered_map>
#include <parallel/algorithm>
#include <numeric>
#include <omp.h>
#include "distribute_kmers.h"


double file_open_time=0.0, global_file_open_time=0.0;
double file_read_time=0.0, global_file_read_time=0.0;
double file_copy_time=0.0, global_file_copy_time=0.0;
double file_close_time=0.0, global_file_close_time=0.0;
double input_process_time=0.0, global_input_process_time=0.0;

void VarLenDivideReads(input_read_data &out, MPI_File &in, const int rank, const int size, uint64_t &nlines) {
    // Open the file and move to the offset from direct division
    MPI_Offset filesize;

    double t1 = MPI_Wtime();
    MPI_File_get_size(in, &filesize);
    MPI_Offset byte_localsize = filesize/size;
    MPI_Offset byte_start = rank * byte_localsize;
    MPI_Offset byte_end = (rank + 1) * byte_localsize;
    if (rank == size-1) byte_end = filesize;

#ifdef DEBUG_PARLLEL_IN
    fprintf(stderr, "filesize=%lu, byte_start=%lu, byte_end=%lu\n", filesize, byte_start, byte_end);
#endif

    // Newlines divide the reads
    // Real starting point:
    //  - if rank 0, beginning of the file
    //  - otherwise, first newline found in byte_local region
    //  -       if no newlines, then no assigned area (empty rank)
    // Real ending point:
    //  - if last rank, end of the file
    //  - otherwise, first newline after current rank's file ending

    // Then, share pos of first newline with neighbor backwards

    // This will consider all "comments" to be reads, which will later be removed
    MPI_Offset pos_of_first_newline = filesize;

    MPI_File_seek(in, byte_start, MPI_SEEK_SET);
    // Find the first newline
    const MPI_Offset buf_size = 1<<14;
    char* buffer = new char[buf_size];
    
    for (MPI_Offset pos = byte_start; pos < byte_end; ) {
        MPI_Status status;
        memset(&status, 0xff, sizeof(MPI_Status));

        MPI_File_read(in, buffer, buf_size, MPI_CHAR, &status);

        int chars_actually_read;
        MPI_Get_count(&status, MPI_CHAR, &chars_actually_read);
        // Go through and find the first newline
        for (size_t buf_pos = 0; buf_pos < chars_actually_read; ++pos, ++buf_pos) {
            if (buffer[buf_pos] == '\n') {
                // We have a newline at "pos"
                pos_of_first_newline = pos;
                break;
            }
        }
        if (pos_of_first_newline != filesize) break;
    }

    MPI_Offset start = (rank == 0) ? byte_start : pos_of_first_newline;
    MPI_Offset end = (rank == size-1) ? byte_end : 0;

    uint64_t pos_of_first_newline_int = pos_of_first_newline;
    uint64_t end_int = end;

    // Now share the pos of first newline
    if (rank % 2 == 0) {
        if (rank != size-1)
            MPI_Recv(&end_int, 1, MPI_UINT64_T, rank+1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        if (rank != 0)
            MPI_Send(&pos_of_first_newline_int, 1, MPI_UINT64_T, rank-1, 1, MPI_COMM_WORLD);
    } else {
        MPI_Send(&pos_of_first_newline_int, 1, MPI_UINT64_T, rank-1, 0, MPI_COMM_WORLD);
        if (rank != size-1)
            MPI_Recv(&end_int, 1, MPI_UINT64_T, rank+1, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

    end = end_int;

#ifdef DEBUG_PARLLEL_IN
    fprintf(stderr, "filesize=%lu, start=%lu, end=%lu\n", filesize, start, end);
#endif

    // Now, proceed with the actual reading
    // First, find the size to allocate
    // Take a pass over the input, and find offsets and sizes
    std::vector<MPI_Offset> read_offsets;
    read_offsets.reserve(1<<25);
    std::vector<size_t> read_sizes;
    read_sizes.reserve(1<<25);

    MPI_File_seek(in, start, MPI_SEEK_SET);
    bool check_for_comment = true;
    bool inside_of_comment = false;
    for (MPI_Offset pos = start; pos < end; ) {
        MPI_Status status = {0xff,};

        MPI_File_read(in, buffer, buf_size, MPI_CHAR, &status);

        int chars_actually_read;
        MPI_Get_count(&status, MPI_CHAR, &chars_actually_read);
        // Check if the starting is correct
        if (rank != 0 && pos == start) {
            assert(buffer[0] == '\n');
        }
        // Go through and find the first newline
        for (size_t buf_pos = 0; buf_pos < chars_actually_read && pos < end; ++pos, ++buf_pos) {

            if (buffer[buf_pos] == '\n') {
                // End any current element, setting the size
                if (inside_of_comment) {
                    inside_of_comment = false;
                } else {
                    if (read_offsets.size() > 0) read_sizes.emplace_back(pos-read_offsets[read_offsets.size()-1]);
                }
                check_for_comment = true;
            } else if (check_for_comment) {
                if (buffer[buf_pos] == '>') {
                    inside_of_comment = true;
                } else {
                    read_offsets.emplace_back(pos);
                }
                check_for_comment = false;
            }
        }
    }
    if (!inside_of_comment && !check_for_comment && read_offsets.size() > 0) read_sizes.emplace_back(end-read_offsets[read_offsets.size()-1]);

    if (read_sizes.size() != read_offsets.size()) {
        fprintf(stderr, "rs=%lu ro=%lu\n", read_sizes.size(), read_offsets.size());
    }
    if (read_sizes.size() != read_offsets.size()) throw -2;

    // Allocate the space to hold the reads
    size_t total_bytes = 0;
    size_t max_read_len = 0;
    for (size_t idx = 0; idx < read_sizes.size(); ++idx) {
        total_bytes += read_sizes[idx];
        if (read_sizes[idx] > max_read_len) max_read_len = read_sizes[idx];
#ifdef DEBUG_PARLLEL_IN
        fprintf(stderr, "read size %lu\n", read_sizes[idx]);
#endif
    }
#ifdef DEBUG_PARLLEL_IN
    fprintf(stderr, "total bytes= %lu\n", total_bytes);
#endif
    out.read_data = (char*)malloc(sizeof(char)*total_bytes);
    if (out.read_data == nullptr) throw std::bad_alloc();

    // Take a second pass over the input, and copy out the data
    char* read_vals = out.read_data;
    check_for_comment = true;
    inside_of_comment = false;
    MPI_File_seek(in, start, MPI_SEEK_SET);
    size_t read_pos = 0;
    for (MPI_Offset pos = start; pos < end; ) {
        MPI_Status status = {0xff,};

        MPI_File_read(in, buffer, buf_size, MPI_CHAR, &status);

        int chars_actually_read;
        MPI_Get_count(&status, MPI_CHAR, &chars_actually_read);
        // Go through and find the first newline
        for (size_t buf_pos = 0; buf_pos < chars_actually_read && pos < end; ++pos, ++buf_pos) {
            if (buffer[buf_pos] == '\n') {
                // End any current element, setting the size
                if (inside_of_comment) {
                    inside_of_comment = false;
                }
                check_for_comment = true;
            } else if (check_for_comment) {
                if (buffer[buf_pos] == '>')
                    inside_of_comment = true;
                check_for_comment = false;
            }

            if (!check_for_comment && !inside_of_comment) {
                assert(read_pos < total_bytes);
                read_vals[read_pos++] = buffer[buf_pos];
            }
        }
    }

    delete [] buffer;

    // Full out the rest of the data in out
    out.read_data_size = total_bytes;
    out.num_of_reads = read_sizes.size();
    //   read_data_ptr
    out.read_data_ptr = (size_t *)malloc(sizeof(size_t)*(out.num_of_reads+1));
    if (out.read_data_ptr == nullptr) throw std::bad_alloc();
    size_t * read_data_ptr = out.read_data_ptr; // allocated with out.alloct
    read_data_ptr[0] = 0;
    for (size_t idx = 0; idx < read_sizes.size(); ++idx)
        read_data_ptr[idx+1] = read_data_ptr[idx] + read_sizes[idx];

    double t2 = MPI_Wtime();
    file_read_time = t2-t1;
    MPI_Reduce(&file_read_time, &global_file_read_time, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if (rank == 0) printf ("Average time for MPI parallel file reading (secs): %f \n",
                            (double)global_file_read_time/(double)size);
    nlines = out.num_of_reads;

#ifdef DEBUG_PARLLEL_IN
    read_pos =0;
    read_vals = out.read_data;
    for (size_t idx = 0; idx < read_sizes.size(); ++idx) {
        char str[10000] = { '\n' };
        strncpy(str, read_vals+read_pos, read_sizes[idx]);
        read_pos += read_sizes[idx];
        str[read_sizes[idx]] = 0;
        fprintf(stderr, "the %luth read: %s\n", idx, str);
    }

    MPI_Barrier(MPI_COMM_WORLD);
    fprintf(stderr, "number of reads: %lu, number of char: %lu\n", out.num_of_reads, out.read_data_size);
    for (size_t idx = 0; idx < read_sizes.size()+1; ++idx) {
        fprintf(stderr, "ptr at %lu: %lu\n", idx, out.read_data_ptr[idx]);
    }

    MPI_Barrier(MPI_COMM_WORLD);
    fprintf(stderr, "finished\n");
#endif
}

char* DivideReads(MPI_File *in, const int rank, const int size, const int overlap,
                 uint64_t *nlines, size_t *data_size) {
    MPI_Offset filesize;
    MPI_Offset localsize, ori_localsize;
    MPI_Offset start;
    MPI_Offset end;
    char *chunk;
    uint64_t pad_width=1000000000; // 1GB

#ifdef DEBUG_INPUT
    char proc_id[3];
    char output_file_name[25];

    sprintf(proc_id, "%d", rank); 
    strcpy(output_file_name,"./output/out_");
    strcpy(&output_file_name[strlen(output_file_name)],proc_id);
    strcpy(&output_file_name[strlen(output_file_name)],".log");

    char debug_file_name[25];
    strcpy(debug_file_name,"./output/debug_p");
    strcpy(&debug_file_name[strlen(debug_file_name)],proc_id);
    strcpy(&debug_file_name[strlen(debug_file_name)],".log");

    FILE *f = fopen(output_file_name, "w");
    if (f == NULL)
    {
        printf("Error opening file!\n");
        exit(1);
    }

    FILE *f0 = fopen(debug_file_name, "w");
    if (f0 == NULL)
    {
        printf("Error opening p0 file!\n");
        exit(1);
    }
#endif

    /* figure out who reads what */

    double t3 = MPI_Wtime ();
    MPI_File_get_size(*in, &filesize);
    localsize = filesize/size;
    start = rank * localsize;
    end   = start + localsize - 1;
 
#ifdef DEBUG_IO
    printf("Proc:%d size of file is: %lld, start: %lld, end: %lld \n", rank, filesize, start, end);
#endif
    /* add overlap to the end of everyone's chunk... */
    end += overlap;

    /* except the last processor, of course */
    if (rank == size-1) end = (filesize-1);

    localsize =  end - start + 1;

    if (rank == 0) fprintf (stderr, "filesize: %llu, localsize: %llu\n", filesize, localsize);

#ifdef DEBUG_IO
    printf("Proc:%d, size of file is: %lld, start: %lld, end: %lld, localsize: %lld \n", rank, filesize, start, end, localsize);
#endif
    /* allocate memory */
    chunk = (char*) malloc( (localsize + 1)*sizeof(char));

    //create contiguous derived data type
    MPI_Datatype rowtype;
    MPI_Type_contiguous(pad_width, MPI_BYTE, &rowtype);
    MPI_Type_commit(&rowtype);

    /* Check the size of localsize */
    ori_localsize=localsize;
    size_t dest_offset=0;
    if (localsize > pad_width)
    {
      int count_size = floor(localsize/pad_width);

      if (rank == 0) fprintf (stderr, "count_size: %d\n", count_size);

      /* everyone reads in their part */
      MPI_File_read_at_all(*in, start, chunk, count_size, rowtype, MPI_STATUS_IGNORE);
      localsize -= (count_size*pad_width);
      
      assert (ori_localsize == ((count_size*pad_width) + localsize));
      start += count_size*pad_width;
      dest_offset = count_size*pad_width;
    }

    MPI_File_read_at_all(*in, start, &chunk[dest_offset], localsize, MPI_BYTE, MPI_STATUS_IGNORE);
    chunk[ori_localsize] = '\0';

    if (rank == 0) fprintf (stderr, "dest_offset: %lu, localsize: %llu\n", dest_offset, localsize);

    // free datatype
    MPI_Type_free(&rowtype);

    localsize = ori_localsize;
    double t4 = MPI_Wtime ();
    file_read_time = t4 - t3;
    MPI_Reduce(&file_read_time, &global_file_read_time, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if (rank == 0) printf ("Average time for MPI_FIle_Read (secs): %f \n", 
                            (double)global_file_read_time/(double)size);

    /*
     * everyone calculate what their actual start and end positions on a read
     * from the first '>' after start to the first '>' after the
     * overlap region starts (eg, after end - overlap + 1)
     */

    uint64_t locstart=0, locend=localsize;

#ifdef DEBUG_INPUT
    /* Print data */
    for (int i=0; i<localsize; i++)
         fprintf(f0,"%c", chunk[i]);
    printf("\n");
        
    fclose(f0);

    printf("Before: Proc:%d locstart: %lld, locend: %lld, localsize: %lld, overlap: %d \n", rank, locstart, locend, localsize, overlap);
#endif

    double t11 = MPI_Wtime();

    if (rank == 0) fprintf (stderr, "initial locstart: %lu, locend: %lu\n", locstart, locend);

    if (rank != 0) {
        while(chunk[locstart] != '>') locstart++;
    }
    if (rank != size-1) {
        locend-=overlap;
        while(chunk[locend] != '>') locend++;
    }
    localsize = locend-locstart;
    
    if (rank == 0) fprintf (stderr, "new localsize: %llu, locstart: %lu, locend: %lu\n", localsize, locstart, locend);

#ifdef DEBUG_INPUT
    printf("After: Proc:%d locstart: %lld, locend: %lld, localsize: %lld \n", rank, locstart, locend, localsize);
#endif

    /* Now let's copy our actual data over into a new array, with no overlaps */
    char *data = (char *)malloc((localsize+1)*sizeof(char));
    memcpy(data, &(chunk[locstart]), localsize);
    data[localsize] = '\0';
    free(chunk);

    /* Now we'll count the number of lines */
    *nlines = 0;
    for (int i=0; i<localsize; i++)
        if (data[i] == '>') (*nlines)++;

    *data_size = localsize;
#ifdef DEBUG_INPUT
    /* Print data to corresponding output files */
    for (int i=0; i<localsize; i++)
         fprintf(f,"%c", data[i]);

    fclose(f);
#endif

    double t12 = MPI_Wtime();
    file_copy_time = t12-t11;
    MPI_Reduce(&file_copy_time, &global_file_copy_time, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if (rank == 0) printf ("Average time for MPI_FIle_copy (secs): %f \n", 
                            (double)global_file_copy_time/(double)size);

    return data;
}


input_read_data perform_input_reading (const int rank, const int size,
                                       std::string &fileName,
                                       uint64_t *num_reads_per_p, uint64_t *paired_nlines)
{

    MPI_File in;
    int ierr;
   
    //const int overlap=read_length+100;
    uint64_t nlines=0, scan_nlines=0;
    size_t read_data_size=0;
    input_read_data input_rdata;

    if (fileName.empty()) return input_rdata;

    double start_t = MPI_Wtime ();

    double t1 = MPI_Wtime ();
    ierr = MPI_File_open(MPI_COMM_WORLD, fileName.c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &in);
    if (ierr) {
        if (rank == 0) fprintf(stderr, "Couldn't open the input FASTA file\n");
        MPI_Finalize();
        exit(2);
    }

    double t2 = MPI_Wtime ();
    file_open_time = t2 - t1;

    MPI_Reduce(&file_open_time, &global_file_open_time, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if (rank == 0) printf ("Average time for MPI_FIle_Open (secs): %f \n", 
                            (double)global_file_open_time/(double)size);

    if (rank==0)
        fprintf(stderr, "Start Reading the input dataset\n");
 
    VarLenDivideReads(input_rdata, in, rank, size, nlines);
    //input_rdata.read_data = DivideReads(&in, rank, size, overlap, &nlines, &read_data_size);
    //input_rdata.read_data_size = read_data_size;

#ifdef DEBUG_OUT
    printf("Rank %d has %d lines\n", rank, nlines);
#endif

    double t5 = MPI_Wtime ();
    MPI_File_close(&in);

    MPI_Barrier(MPI_COMM_WORLD);
    double t6 = MPI_Wtime ();
    file_close_time = t6 - t5;
    MPI_Reduce(&file_close_time, &global_file_close_time, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if (rank == 0) printf ("Average time for MPI_FIle_Close (secs): %f \n", 
                            (double)global_file_close_time/(double)size);

    double end_t = MPI_Wtime ();
    input_process_time = end_t - start_t;

    MPI_Barrier(MPI_COMM_WORLD);

    MPI_Reduce(&input_process_time, &global_input_process_time, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if (rank == 0) printf ("Average time for reading and storing the Reads in each proc's memory (secs): %f \n", 
                            (double)global_input_process_time/(double)size);

    uint64_t global_reads=0;
    MPI_Allreduce (&nlines, &global_reads, 1, MPI_UINT64_T, MPI_SUM, MPI_COMM_WORLD);
    if (rank == 0) printf ("Total number of reads: %lu\n", global_reads);
    if (rank == 0) fprintf (stderr, "Total number of reads: %lu\n", global_reads);
    
    *num_reads_per_p=nlines;
    *paired_nlines=global_reads;

    if (rank==0)
        fprintf(stderr, "Completed Reading the input dataset\n");

    return input_rdata;
}
