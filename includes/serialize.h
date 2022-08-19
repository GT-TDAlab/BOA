#ifndef SERIALIZE_H
#define SERIALIZE_H

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

template<typename POD>
std::ostream& serialize(std::ostream& os, const POD& v);


template<typename POD1, typename POD2>
std::ostream& serialize(std::ostream& os, const std::pair<POD1, POD2>& p);

template<typename POD>
std::ostream& serialize(std::ostream& os, std::vector<POD> const& vec);

std::ostream& serialize(std::ostream& os, std::vector<std::pair<int, int>> const& pr);

std::ostream& serialize(std::ostream& os, std::vector<bool> const& tvec);

std::ostream& serialize(std::ostream& os, const KmerBucket &bi);



template<typename POD>
std::istream& deserialize(std::istream& is, POD& v);

template<typename POD1, typename POD2>
std::istream& deserialize(std::istream& is, std::pair<POD1, POD2>& p);

template<typename POD>
std::istream& deserialize(std::istream& is, std::vector<POD>& vec);


std::istream& deserialize(std::istream& is, std::vector<std::pair<int, int>>& pr);

std::istream& deserialize(std::istream& is, std::vector<bool>& tvec);

std::istream& deserialize(std::istream& is, KmerBucket &bi);



#endif
