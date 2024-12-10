#ifndef CB_GLOB_HH
#define CB_GLOB_HH
#pragma once

#include "dfinfra/standard_includes.hh"


/*
 * Global Control Block for Experiment 8
 */
class CbGlob {
  private:
    CbGlob(const CbGlob&) = delete;
    CbGlob& operator=(const CbGlob&) = delete;
  public:
    CbGlob();
    ~CbGlob();
  public:
    inline uint max_no_rep() const { return 17; }
           uint get_no_rep(const uint) const;

    // const getters by value
    inline uint min_log_card_key() const { return _min_log_card_key; }
    inline uint max_log_card_key() const { return _max_log_card_key; }
    inline uint min_log_card_fk()  const { return _min_log_card_fk; }
    inline uint max_log_card_fk()  const { return _max_log_card_fk; }
    inline uint min_log_dom_fk_reduction()  const { return _min_log_dom_fk_reduction; }
    inline uint max_log_dom_fk_reduction()  const { return _max_log_dom_fk_reduction; }
    inline uint min_log_domsize_fk() const { return _min_log_domsize_fk; }
    inline uint max_log_domsize_fk() const { return _max_log_domsize_fk; }
    inline uint min_skew() const { return _min_skew; }
    inline uint max_skew() const { return _max_skew; }
    inline double card_key_factor() const { return _card_key_factor; }
    inline double card_fk_factor()  const { return _card_fk_factor; }
    inline uint hwtno()    const { return _hwtno; }
    inline const std::string& directory() const { return _directory; }
    inline const std::string& filename() const { return _filename; }
    inline const std::string& measure() const { return _measure; } 
    inline std::ostream& os() const { return (*_os); }
    inline std::ostream* os_ptr() const { return _os; }
    inline uint64_t trace() const { return _trace; }
    inline bool     help()  const { return _help; }
  public:
    // non-const getters by reference (to set value from outside, e.g., by CLI11)
    inline uint& min_log_card_key() { return _min_log_card_key; }
    inline uint& max_log_card_key() { return _max_log_card_key; }
    inline uint& min_log_card_fk()  { return _min_log_card_fk; }
    inline uint& max_log_card_fk()  { return _max_log_card_fk; }
    inline uint& min_log_dom_fk_reduction() { return _min_log_dom_fk_reduction; }
    inline uint& max_log_dom_fk_reduction() { return _max_log_dom_fk_reduction; }
    inline uint& min_log_domsize_fk() { return _min_log_domsize_fk; }
    inline uint& max_log_domsize_fk() { return _max_log_domsize_fk; }
    inline uint& min_skew() { return _min_skew; }
    inline uint& max_skew() { return _max_skew; }
    inline double& card_key_factor() { return _card_key_factor; }
    inline double& card_fk_factor()  { return _card_fk_factor; }
    inline uint& hwtno()    { return _hwtno; }
    inline std::string& directory() { return _directory; }
    inline std::string& filename() { return _filename; }
    inline std::string& measure() { return _measure; }
    inline std::ostream& os() { return (*_os); }
    inline std::ostream* os_ptr() { return _os; }
    inline uint64_t& trace() { return _trace; }
    inline bool&     help()  { return _help; }
  public:
    void min_log_card_key(const uint&);
    void max_log_card_key(const uint&);
    void min_log_card_fk(const uint&);
    void max_log_card_fk(const uint&);
    void min_log_dom_fk_reduction(const uint&);
    void max_log_dom_fk_reduction(const uint&);
    void min_log_domsize_fk(const uint&);
    void max_log_domsize_fk(const uint&);
    void min_skew(const uint&);
    void max_skew(const uint&);
    void card_key_factor(const double&);
    void card_fk_factor(const double&);
    void hwtno(const uint&); 
    void directory(const std::string&);
    void filename(const std::string&);
    void measure(const std::string&);
    void help(const bool&);
  public:
    std::string to_string() const;
    std::string complete_filename() const;
  public:
    bool alloc_os();
  public:
    std::ostream& print(std::ostream& os, const bool aPrintStrings = true) const;
  private:
    uint _min_log_card_key;  // minimum of log2 of cardinality of key relation  R
    uint _max_log_card_key;  // maximum of log2 of cardinality of key relation  R
    uint _min_log_card_fk;   // minimum of log2 of cardinality of foreign key relation S
    uint _max_log_card_fk;   // maximum of log2 of cardinality of foreign key relation S
    uint _min_log_dom_fk_reduction; // minimum of log2 join attribute domain reduction of S
    uint _max_log_dom_fk_reduction; // maximum of log2 join attribute domain reduction of S
    uint _min_log_domsize_fk; // minimum of log2 domain size S
    uint _max_log_domsize_fk; // maximum of log2 domain size S
    uint _min_skew;           // minimum skew considered for join attribute of foreign key relation S
    uint _max_skew;           // maximum skew considered for join attribute of foreign key relation S
    double _card_key_factor;  // scaling factor for key cardinality: card_key := 2^log_card_key * _card_key_factor
    double _card_fk_factor;   // scaling factor for fk cardinality
    uint _hwtno;              // hardware thread number where executable should run
    std::string   _directory; // output directory
    std::string   _filename;  // output filename
    std::string   _measure;   // what to measure: <build|probe>_<bun|bnu>
    std::ostream* _os;        // output stream
    uint64_t      _trace;     // trace
    bool          _help;      // help
};

#endif
