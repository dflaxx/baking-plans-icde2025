#include "cb_glob.hh"

#include <fstream>
#include <string>
#include <algorithm>

CbGlob::CbGlob()
       : _min_log_card_key(0),
         _max_log_card_key(29),
         _min_log_card_fk(0),
         _max_log_card_fk(29),
         _min_log_dom_fk_reduction(0),
         _max_log_dom_fk_reduction(6),
         _min_log_domsize_fk(0),
         _max_log_domsize_fk(29),
         _min_skew(0),
         _max_skew(1),
         _card_key_factor(1.0),
         _card_fk_factor(1.0),
         _hwtno(0),
         _directory("."),
         _filename(""),
         _measure(""),
         _os(nullptr),
         _trace(0),
         _help(false) {
}

CbGlob::~CbGlob() {
  if(nullptr != _os) {
    delete _os;
  }
}


uint
CbGlob::get_no_rep(const uint x) const {
  if      (10 >= x) { return 17; }
  else if (12 >= x) { return 11; }
  else if (14 >= x) { return 7; }
  else if (16 >= x) { return 5; }
  else              { return 3; }
}

void CbGlob::min_log_card_key(const uint& x) { _min_log_card_key = x; }
void CbGlob::max_log_card_key(const uint& x) { _max_log_card_key = x; }
void CbGlob::min_log_card_fk(const uint& x)  { _min_log_card_fk = x; }
void CbGlob::max_log_card_fk(const uint& x)  { _max_log_card_fk = x; }
void CbGlob::min_log_dom_fk_reduction(const uint& x) { _min_log_dom_fk_reduction = x; }
void CbGlob::max_log_dom_fk_reduction(const uint& x) { _max_log_dom_fk_reduction = x; }
void CbGlob::min_log_domsize_fk(const uint& x) { _min_log_domsize_fk = x; }
void CbGlob::max_log_domsize_fk(const uint& x) { _max_log_domsize_fk = x; }
void CbGlob::min_skew(const uint& x) { _min_skew = x; }
void CbGlob::max_skew(const uint& x) { _max_skew = x; }
void CbGlob::card_key_factor(const double& x) { _card_key_factor = x; }
void CbGlob::card_fk_factor(const double& x) { _card_fk_factor = x; }
void CbGlob::hwtno(const uint& x) { _hwtno = x; }
void CbGlob::directory(const std::string& x) { _directory = x; }
void CbGlob::filename(const std::string& x) { _filename = x; }
void CbGlob::measure(const std::string& x) { _measure = x; }
void CbGlob::help(const bool& x) { _help = x; }

std::string
CbGlob::to_string() const {
  // don't have trailing zeros for factors
  std::stringstream lFactorsStream;
  lFactorsStream << std::noshowpoint
    << "k" << card_key_factor() << '_'
    << "f" << card_fk_factor();
  std::string lFactorsStr = lFactorsStream.str();
  lFactorsStr.erase(std::remove(lFactorsStr.begin(), lFactorsStr.end(), '.'), lFactorsStr.end());

  std::string lRes;
  lRes +=         std::to_string(hwtno()) + '_'
          + 'K' + std::to_string(min_log_card_key())
          + '_' + std::to_string(max_log_card_key()) + '_'
          + 'F' + std::to_string(min_log_card_fk())
          + '_' + std::to_string(max_log_card_fk()) + '_'
          + 'R' + std::to_string(min_log_dom_fk_reduction())
          + '_' + std::to_string(max_log_dom_fk_reduction()) + '_'
          + 'D' + std::to_string(min_log_domsize_fk())
          + '_' + std::to_string(max_log_domsize_fk()) + '_'
          + 'S' + std::to_string(min_skew())
          + '_' + std::to_string(max_skew()) + '_'
          + lFactorsStr;
  return lRes;
}

std::string
CbGlob::complete_filename() const {
  return  (directory() + '/' + measure() + '_' + filename() + "_X_" + to_string() + ".rel");
}

bool
CbGlob::alloc_os() {
  if(nullptr != _os) {
    return false;
  }
  const std::string lFilename = complete_filename();
  _os = new std::ofstream(lFilename);
  if(nullptr == _os) { 
    std::cout << "Can't open file '" << lFilename << "'." << std::endl;
    return false; 
  }
  if(!(*_os)) {
    std::cout << "Can't open file '" << lFilename << "'." << std::endl;
    delete _os;
    return false;
  }
  return true;
}

std::ostream&
CbGlob::print(std::ostream& os, const bool aPrintStrings) const {
  os << "to measure: " << measure() << std::endl
     << "min/max log card key: " << min_log_card_key() << ',' << max_log_card_key() << std::endl
     << "min/max log card fk : " << min_log_card_fk()  << ',' << max_log_card_fk()  << std::endl
     << "min/max dom red  fk : " << min_log_dom_fk_reduction()  << ',' << max_log_dom_fk_reduction() << std::endl
     << "min/max dom size fk : " << min_log_domsize_fk()  << ',' << max_log_domsize_fk() << std::endl
     << "min/max dom skew fk : " << min_skew()  << ',' << max_skew()  << std::endl
     << "key factor          : " << card_key_factor() << std::endl
     << "fk factor           : " << card_fk_factor() << std::endl;
  if(aPrintStrings) {
    os << "out        directory: " << directory() << std::endl
       << "out  filename prefix: " << filename() << std::endl
       << "out  filename       : " << complete_filename() << std::endl
       << "run  on      hwtno  : " << hwtno() << std::endl
       ;
  }

  return os;
}

