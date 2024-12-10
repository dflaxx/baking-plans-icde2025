#pragma once

/* 
 * A physical algebra to demonstrate and measure hash joins.
 * This is version 2, a reimplementation of the old physical algebra (version 1).
 *
 * Properties:
 * - push-based (producer and comsumer operator)
 * - materialization granularity: tuple-wise
 * - templated
 * - storage and processing format: row store
 *
 * Operators:
 * - Top
 * - Table Scan
 * - Chaining Hash Join with chained hash table (build + probe)
 * - 3D Hash Join with nested/3D hash table (build + probe/unnest)
 * - Maybe TODO: 3D Hash Join where probe and unnest are separate operators.
 * - Maybe TODO: Materialization
 *
 * All operators inherit from a common base class AlgOpBase2 which provides some
 * debugging functionalities. Otherwise, there is no inheritance and also no
 * overriding of virtual functions.
 * Instead, the algebra is templated, i.e., all operators are class templates.
 * To ensure that only compatible types are combined, C++20 concepts are used.
 *
 * Operator interface specification
 *
 * Consumers
 * - mem_alloc
 * - mem_init
 * - init
 * - step
 * - fin
 * - clear
 * - free
 *
 * Producers
 * - run()
 *
 * Storage and processing format:
 * See RelationRS.hh.
 * - Row store relation RelationRS2, using a FlexVectorNC.
 * - structs for tuples, hash functions and predicates.
 * XXX Assumption: Can only process tuples (binary units, BUNs) of the form [T k, T a]
 *                 where T is either uint32_t or uint64_t.
 *                 k is the key attribute, a is the foreign key attribute.
 *
 */

#include "algop_v2_base.hh"
#include "algop_v2_top.hh"
#include "algop_v2_scan.hh"
#include "algop_v2_join_hash_ch_orig.hh"
#include "algop_v2_join_hash_3d_orig.hh"
#include "algop_v2_join_hash_ch_amac.hh"
#include "algop_v2_join_hash_3d_amac.hh"
#include "algop_v2_join_hash_ch_rp.hh"
#include "algop_v2_join_hash_3d_rp.hh"

#include "RelationRS.hh"
