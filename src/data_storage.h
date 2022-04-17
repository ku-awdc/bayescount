#ifndef DATA_STORAGE_H
#define DATA_STORAGE_H

// This is required by distribution.h and precludes Rcpp.h:
#include <RcppDist.h>
#include <array>

#include "enums.h"
#include "fecrt.h"

// The data_storage class takes data as input and stores it
// The data_summary class calculates mean/var/k as needed

template<bool t_paired, typename t_cont_type, containers t_container>
class data_storage;

template<bool t_paired, ktypes t_ktype>
class data_summary;

// Paired data storage:
template<>
class data_storage<true>
{
private:
  // Note: these are not actually needed if !s_store_count
  t_cont_type m_pre;
  t_cont_type m_post;
  
public:
  const size_t m_maxN;
};



#endif // DATA_STORAGE_H