///
/// @file solomon/acc.hpp
/// @author Yohei MIKI (The University of Tokyo)
/// @brief pragmas to use GPU by OpenACC directives
///
/// @copyright Copyright (c) 2024 Yohei MIKI
///
/// The MIT License is applied to this software, see LICENSE.txt
///
#if !defined(SOLOMON_ACC_HPP)
#define SOLOMON_ACC_HPP

#if defined(OFFLOAD_BY_OPENACC) && defined(_OPENACC)
// set OpenACC environments
#include "acc/clause.hpp"
#include "acc/directive.hpp"
#include "acc/runtime.hpp"

// set utilities for OpenACC
#include "acc/utils.hpp"

// load converter from OpenMP target to OpenACC
#include "convert/omp2acc/clause.hpp"
#include "convert/omp2acc/directive.hpp"

// set utilities for OpenMP target
#include "omp/target/utils.hpp"
#endif  // defined(OFFLOAD_BY_OPENACC) && defined(_OPENACC)

#endif  // !defined(SOLOMON_ACC_HPP)
