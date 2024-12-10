#include "line.hh"
#include "dfinfra/exceptions.hh"
#include <format>

std::ostream&
rt_stat_t::print(std::ostream& os, const char aSep) const {
  os << _cc_min << aSep
     << _cc_max << aSep
     << _cc_med << aSep
     << _cc_avg << aSep
     << _cc_mx2 << aSep
     << _cc_fst;
  return os;
}

std::ostream&
line_t::print(std::ostream& os, const char aSep) const {
  os << _id               << aSep // (1)
     << _no_rep_key       << aSep // (2)
     << _no_rep_fk        << aSep // (3)
     << _attr_size        << aSep // (4)
     << _log_card_fk      << aSep // (5)
     << _card_fk          << aSep // (6)
     << _log_dom_red_fk   << aSep // (7)
     << _log_dom_fk       << aSep // (8)
     << _dom_fk           << aSep // (9)
     << _skew_fk          << aSep // (10)
     << _log_card_key     << aSep // (11)
     << _card_key         << aSep // (12)
     << _nodv_fk_hll      << aSep // (13)
     << _hashdir_size_fk  << aSep // (14)
     << _hashdir_size_key << aSep // (15)
     << _card_sj_sr       << aSep // (16)
     << _card_result      << aSep // (17)
     << _time;                    // (18)
  // cc_...
  for(uint i = 0; i < kNoPlan; ++i) {
    os << aSep << _rt_stat_build[i];
  }
  for(uint i = 0; i < kNoPlan; ++i) {
    os << aSep << _rt_stat_probe[i];
  }
  // mem_...
  for(uint i = 0; i < kNoPlan; ++i) {
    os << aSep << _mem_stat_ht.at(i);
  }
  // nomn...  // number of main nodes
  for(uint i = 0; i < kNoPlan; ++i) {
    os << aSep << _ht_num_nodes_main.at(i);
  }
  // nosn...  // number of sub nodes (3d HT only)
  for(uint i = 0; i < kNoPlan; ++i) {
    os << aSep << _ht_num_nodes_sub.at(i);
  }
  // dir_fill_...  // fraction of non-empty hash table directory entries
  for(uint i = 0; i < kNoPlan; ++i) {
    os << aSep << _ht_dir_fill.at(i);
  }
  // cnt_bld_...
  for(uint i = 0; i < kNoPlan; ++i) {
    os << aSep << _cnt_build.at(i);
  }
  // cnt_prb_...
  for(uint i = 0; i < kNoPlan; ++i) {
    os << aSep << _cnt_top.at(i);
  }
  return os;
}

std::ostream&
line_t::print_build_rs(std::ostream& os, const char* aSep) const {
  // TODO
  throw df::infra::NotImplementedException();
  return os;
}

std::ostream&
line_t::print_build_sr(std::ostream& os, const char* aSep) const {
  // TODO
  throw df::infra::NotImplementedException();
  return os;
}

std::ostream&
line_t::print_probe_rs(std::ostream& os, const char* aSep) const {
  // TODO
  throw df::infra::NotImplementedException();
  return os;
}

std::ostream&
line_t::print_probe_sr(std::ostream& os, const char* aSep) const {
  // TODO
  throw df::infra::NotImplementedException();
  return os;
}


std::ostream&
line_t::print_header_line(std::ostream& os) {
  static const char* gNst[2] = { "ch", "3d" };
  //static const char* gImc[2][2] = { { "orig", "impr" }, { "amac", "imac" } };
  /*
   * orig = upk nopf
   * impr = pkd nopf
   * uprp = upk rp
   * pkrp = pkd rp
   * amac = upk amac
   * imac = pkd amac
   */
  static const char* gImc[3][2] = { { "orig", "impr" }, { "uprp", "pkrp" }, { "amac", "imac" } };
  static const char* gJor[2] = { "rs", "sr" };
  static const char* gBld[2] = { "bld", "prb" };
  static const char* gAgg[6] = { "min", "max", "med", "avg", "fst", "mx2" };
  os << "# table R { "
     << "int id, "            // (1)
     << "int n_rep_key, "     // (2)
     << "int n_rep_fk, "      // (3)
     << "int attr_size, "     // (4)  attribute byte size
     << "double l_fk, "       // (5)
     << "int c_fk, "          // (6)
     << "int red_dom, "       // (7)
     << "double l_dom, "      // (8)
     << "int c_dom, "         // (9)
     << "int skew, "          // (10)
     << "double l_key, "      // (11)
     << "int c_key, "         // (12)
     << "int ndv_fk_hll, "    // (13)
     << "int hts_fk, "        // (14)
     << "int hts_key, "       // (15)
     << "int c_mid, "         // (16)
     << "int c_res, "         // (17)
     << "double time";        // (18)
  {
    // debug
    std::cout << "# [debug] ouput header vs. plan ID" << std::endl;
    uint i = 0;
    std::cout << std::format("#         {:20s} {:20s}", "header", "PlanID") << std::endl;
    for(uint lNested = 0; lNested <= 1; ++lNested) {
      for(uint lAmac = 0; lAmac <= 2; ++lAmac) {
        for(uint lImpr = 0; lImpr <= 1; ++lImpr) {
          for(uint lJoinOrder = 0; lJoinOrder < kNoJoinOrder; ++lJoinOrder) {
            std::stringstream s;
            s << "XX" << '_' << gNst[lNested] << '_' << gImc[lAmac][lImpr] << '_' << gJor[lJoinOrder];
            std::cout << std::format("#         {:20s} {:20s}", s.str(), to_string(static_cast<plan_et>(i))) << std::endl;
            ++i;
          }
        }
      }
    }
  }
  // cc_...
  for(uint lBldPrb = 0; lBldPrb < 2; ++lBldPrb) {
    for(uint lNested = 0; lNested <= 1; ++lNested) {
      for(uint lAmac = 0; lAmac <= 2; ++lAmac) {
        for(uint lImpr = 0; lImpr <= 1; ++lImpr) {
          for(uint lJoinOrder = 0; lJoinOrder < kNoJoinOrder; ++lJoinOrder) {
            for(uint lAgg = 0; lAgg < 6; ++lAgg) {
              os << ", ";
              if(2 > lAgg) {
                os << "int";
              } else {
                os << "double";
              }
              os << " cc"
                 << '_' << gAgg[lAgg]
                 << '_' << gNst[lNested]
                 << '_' << gImc[lAmac][lImpr]
                 << '_' << gJor[lJoinOrder]
                 << '_' << gBld[lBldPrb];
            }
          }
        }
      }
    }
  }
  // mem_...
  for(uint lNested = 0; lNested <= 1; ++lNested) {
    for(uint lAmac = 0; lAmac <= 2; ++lAmac) {
      for(uint lImpr = 0; lImpr <= 1; ++lImpr) {
        for(uint lJoinOrder = 0; lJoinOrder < kNoJoinOrder; ++lJoinOrder) {
            os << ", ";
            os << "int";
            os << " mem"
               << '_' << gNst[lNested]
               << '_' << gImc[lAmac][lImpr]
               << '_' << gJor[lJoinOrder];
        }
      }
    }
  }
  // nomn...  // number of main nodes
  for(uint lNested = 0; lNested <= 1; ++lNested) {
    for(uint lAmac = 0; lAmac <= 2; ++lAmac) {
      for(uint lImpr = 0; lImpr <= 1; ++lImpr) {
        for(uint lJoinOrder = 0; lJoinOrder < kNoJoinOrder; ++lJoinOrder) {
            os << ", ";
            os << "int";
            os << " nomn"
               << '_' << gNst[lNested]
               << '_' << gImc[lAmac][lImpr]
               << '_' << gJor[lJoinOrder];
        }
      }
    }
  }
  // nosn...  // number of sub nodes (3d HT only)
  for(uint lNested = 0; lNested <= 1; ++lNested) {
    for(uint lAmac = 0; lAmac <= 2; ++lAmac) {
      for(uint lImpr = 0; lImpr <= 1; ++lImpr) {
        for(uint lJoinOrder = 0; lJoinOrder < kNoJoinOrder; ++lJoinOrder) {
            os << ", ";
            os << "int";
            os << " nosn"
               << '_' << gNst[lNested]
               << '_' << gImc[lAmac][lImpr]
               << '_' << gJor[lJoinOrder];
        }
      }
    }
  }
  // dir_fill...
  for(uint lNested = 0; lNested <= 1; ++lNested) {
    for(uint lAmac = 0; lAmac <= 2; ++lAmac) {
      for(uint lImpr = 0; lImpr <= 1; ++lImpr) {
        for(uint lJoinOrder = 0; lJoinOrder < kNoJoinOrder; ++lJoinOrder) {
            os << ", ";
            os << "double";
            os << " dir_fill"
               << '_' << gNst[lNested]
               << '_' << gImc[lAmac][lImpr]
               << '_' << gJor[lJoinOrder];
        }
      }
    }
  }
  // cnt_bld_...
  for(uint lNested = 0; lNested <= 1; ++lNested) {
    for(uint lAmac = 0; lAmac <= 2; ++lAmac) {
      for(uint lImpr = 0; lImpr <= 1; ++lImpr) {
        for(uint lJoinOrder = 0; lJoinOrder < kNoJoinOrder; ++lJoinOrder) {
            os << ", ";
            os << "int";
            os << " cnt_bld"
               << '_' << gNst[lNested]
               << '_' << gImc[lAmac][lImpr]
               << '_' << gJor[lJoinOrder];
        }
      }
    }
  }
  // cnt_prb...
  for(uint lNested = 0; lNested <= 1; ++lNested) {
    for(uint lAmac = 0; lAmac <= 2; ++lAmac) {
      for(uint lImpr = 0; lImpr <= 1; ++lImpr) {
        for(uint lJoinOrder = 0; lJoinOrder < kNoJoinOrder; ++lJoinOrder) {
            os << ", ";
            os << "int";
            os << " cnt_prb"
               << '_' << gNst[lNested]
               << '_' << gImc[lAmac][lImpr]
               << '_' << gJor[lJoinOrder];
        }
      }
    }
  }
  os << "};" << std::endl;
  return os;
}



std::ostream&
operator<<(std::ostream& os, const rt_stat_t& aRtStat) {
  return aRtStat.print(os, '|');
}

std::ostream& 
operator<<(std::ostream& os, const line_t& aLine) {
  return aLine.print(os, '|');
}


