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
#include <stdlib.h>
#include <assert.h>
#include <errno.h>
#include <stdbool.h>
#include <math.h>
#include <time.h>
#include <omp.h>
#include <mpi.h>
#include <iostream>
#include <sstream>
#include <random>
#include <utility>
#include "distribute_kmers.h"
#include "serialize.h"
#include "bucketing.h"

extern int rank, size;

template<typename POD>
std::ostream& serialize(std::ostream& os, const POD& v)
{
    static_assert(std::is_trivial<POD>::value && std::is_standard_layout<POD>::value,
        "Can only serialize POD types with this function");

  os.write(reinterpret_cast<const char*>(&v), sizeof(POD));
  return os;
}

template<typename POD1, typename POD2>
std::ostream& serialize(std::ostream& os, const std::pair<POD1, POD2>& p)
{
    static_assert(std::is_trivial<POD1>::value && std::is_standard_layout<POD1>::value,
        "Can only serialize POD types with this function");

    static_assert(std::is_trivial<POD2>::value && std::is_standard_layout<POD2>::value,
        "Can only serialize POD types with this function");

   serialize(os, std::get<0>(p));
   serialize(os, std::get<1>(p));
   return os;

}

template<typename POD>
std::ostream& serialize(std::ostream& os, std::vector<POD> const& vec)
{
    // this only works on built in data types (PODs)
    static_assert(std::is_trivial<POD>::value && std::is_standard_layout<POD>::value,
        "Can only serialize POD types with this function");

    auto size = vec.size();
    os.write(reinterpret_cast<char const*>(&size), sizeof(size));
    for (auto& v:  vec) {
         serialize(os, v);
    }
    return os;
}


std::ostream& serialize(std::ostream& os, std::vector<std::pair<int, int>> const& pr)
{
  auto size = pr.size();
  os.write(reinterpret_cast<char const*>(&size), sizeof(size));
  for (auto& v:  pr) {
       serialize(os, v);
  }
  return os;

}

std::ostream& serialize(std::ostream& os, std::vector<bool> const& tvec)
{
  auto size = tvec.size();
  os.write(reinterpret_cast<char const*>(&size), sizeof(size));
  for(std::vector<bool>::size_type i = 0; i < size;)
  {
        unsigned char aggr = 0;
        for(unsigned char mask = 1; mask > 0 && i < size; ++i, mask <<= 1)
            if(tvec.at(i))
                aggr |= mask;
        os.write(reinterpret_cast<char const*>(&aggr), sizeof(unsigned char));
  }
  return os;

}



std::ostream& serialize(std::ostream& os, const KmerBucket &bi)
{
  serialize(os, bi.seq);
  serialize(os, bi.ReadIds);
  return os;

}



template<typename POD>
std::istream& deserialize(std::istream& is, POD& v)
{
  static_assert(std::is_trivial<POD>::value && std::is_standard_layout<POD>::value,
        "Can only deserialize POD types with this function");

  is.read(reinterpret_cast<char*>(&v), sizeof(POD));
  return is;
}

template<typename POD1, typename POD2>
std::istream& deserialize(std::istream& is, std::pair<POD1, POD2>& p)
{
  static_assert(std::is_trivial<POD1>::value && std::is_standard_layout<POD1>::value,
        "Can only deserialize POD types with this function");

  static_assert(std::is_trivial<POD2>::value && std::is_standard_layout<POD2>::value,
        "Can only deserialize POD types with this function");

  deserialize(is, std::get<0>(p));
  deserialize(is, std::get<1>(p));
  return is;

}

template<typename POD>
std::istream& deserialize(std::istream& is, std::vector<POD>& vec)
{
    static_assert(std::is_trivial<POD>::value && std::is_standard_layout<POD>::value,
        "Can only deserialize POD types with this function");

    decltype(vec.size()) size;
    is.read(reinterpret_cast<char*>(&size), sizeof(size));
    vec.resize(size);
    for (auto& v: vec) {
         deserialize(is, v);
    }
    return is;
}



std::istream& deserialize(std::istream& is, std::vector<std::pair<int, int>>& pr)
{
  decltype(pr.size()) size;
  is.read(reinterpret_cast<char*>(&size), sizeof(size));
  pr.resize(size);
  for (auto& v: pr) {
       deserialize(is, v);
  }
  return is;

}

std::istream& deserialize(std::istream& is, std::vector<bool>& tvec)
{
  decltype(tvec.size()) size;
  is.read(reinterpret_cast<char*>(&size), sizeof(size));
  tvec.resize(size);
  for(std::vector<bool>::size_type i = 0; i < size;)
  {
        unsigned char aggr;
        is.read(reinterpret_cast<char*>(&aggr), sizeof(unsigned char));
        for(unsigned char mask = 1; mask > 0 && i < size; ++i, mask <<= 1)
            tvec.at(i) = aggr & mask;
  }
  return is;

}


std::istream& deserialize(std::istream& is, KmerBucket &bi)
{
  deserialize(is, bi.seq);
  deserialize(is, bi.ReadIds);
  return is;

}
