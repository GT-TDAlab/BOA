#ifndef DEFINES_H
#define DEFINES_H

#include <stdint.h>
typedef uint64_t read_t;

#define TIME_PROFILE

//hash function based on read location
extern std::vector<uint64_t> global_reads_disp;
int retrieve_read_proc_owner(read_t read_id);

#endif
