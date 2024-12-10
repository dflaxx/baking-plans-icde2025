#ifndef HYPERLOGLOG_ERTL_HH
#define HYPERLOGLOG_ERTL_HH

#include "dfinfra/standard_includes.hh"
#include "gminfra/bit_intrinsics.hh"
#include "gminfra/vector_ops.hh"
#include <cmath>
#include <limits>

template<typename Tuint>
class HllErtl {
  public:
    typedef Tuint   hval_t; // type of hash value to be inserted
    typedef uint8_t reg_t;  // type of register
    typedef std::vector<reg_t> reg_vt; // type of register vector [need better memory management]
  public:
    HllErtl();
    HllErtl(const uint aLog2NoReg);
    virtual ~HllErtl();
  public:
    virtual void init(const uint aLog2NoReg);
  public:
    virtual void insert(const hval_t aVal);
  public:
    hval_t estimate() const { return estimate_pois(); }
  public:
    void   union_with(const HllErtl& aHll);
    void   minim_with(const HllErtl& aHll);
  public:
    double estimate_orig() const; // original estimator
    double estimate_pois() const; // poisson estimator // best one
    hval_t estimate_maxl() const; // maximum likelihood estimator
    hval_t estimate_linc() const; // linear counting estimator
  public:
    hval_t no_reg_zero() const;
    hval_t no_reg_nonzero() const;
  public:
    void   extract_counts(uint_vt& C) const;
    uint   min_idx_gt_zero(const uint_vt& C) const;
    uint   max_idx_gt_zero(const uint_vt& C) const;
  public:
    double estimate_raw() const;
  public:
    inline uint64_t no_insert() const { return _no_insert; }
    inline uint     no_bit() const { return _no_bit; }
    inline uint     log2_noreg() const { return _log2_noreg; }
    inline uint     no_reg() const { return _regs.size(); }
    inline uint     get_p() const { return _log2_noreg; }
    inline uint     get_q() const { return no_bit() - get_p(); }
    inline uint     shift() const { return _shift; }
    inline hval_t   mask() const { return _mask; }
    inline reg_t    reg(const uint i) const { return _regs[i]; }
    inline const reg_vt& regs() const { return _regs; }
           double   alpha_m() const;
    inline double   alpha_inf() const { return (1/(2*std::log(2))); }
           double   sigma(double x) const;
           double   tau(double x) const;
  public:
    void clear(); // only sets _no_insert = 0, _regs[i] = 0
  public:
    std::ostream& print(std::ostream& os) const;
  private:
    uint64_t _no_insert;
    uint     _no_bit; // total number of bits = 32 or 64
    uint     _log2_noreg; // p
    uint     _shift;
    hval_t   _mask;
    reg_vt   _regs; 
};

template<typename Tuint>
HllErtl<Tuint>::HllErtl()
               : _no_insert(0),
                 _no_bit(sizeof(Tuint) * 8),
                 _log2_noreg(),
                 _shift(),
                 _mask(),
                 _regs() {
}



template<typename Tuint>
HllErtl<Tuint>::HllErtl(const uint aLog2NoReg)
               : _no_insert(0),
                 _no_bit(sizeof(Tuint) * 8),
                 _log2_noreg(aLog2NoReg),
                 _shift(aLog2NoReg),
                 _mask((((Tuint) 1) << aLog2NoReg) - 1),
                 _regs(1 << aLog2NoReg) {
  for(uint i = 0; i < no_reg(); ++i) {
    _regs[i] = 0;
  }
}

template<typename Tuint>
HllErtl<Tuint>::~HllErtl() {
}


template<typename Tuint>
void
HllErtl<Tuint>::init(const uint aLog2NoReg) {
  _no_insert  = 0;
  _no_bit     = sizeof(Tuint) * 8;
  _log2_noreg = aLog2NoReg;
  _shift      = aLog2NoReg;
  _mask       = (((hval_t) 1) << aLog2NoReg) - 1;
  _regs.resize(1 << aLog2NoReg);
  for(uint i = 0; i < no_reg(); ++i) {
    _regs[i] = 0;
  }
}

template<typename Tuint>
void
HllErtl<Tuint>::insert(const hval_t aVal) {
  ++_no_insert;
  const hval_t a = aVal & mask();
  const hval_t b = (aVal >> shift());
  const hval_t k = (0 == b) ? get_q() + 1 : 1 + idx_lowest_bit_set<hval_t>(aVal >> shift());
  // std::cout << "HllErtl::insert: " << aVal << " 0x" << std::hex << aVal << ' ';
  // std::cout << "a = " << std::dec << a << ", b = 0x" << std::hex << (aVal >> shift()) << ", k = " << k;
  // std::cout << std::endl;
  // std::cout << std::dec;
  assert( a < no_reg() );
  if(k > reg(a)) {
    _regs[a] = k;
  }
}

template<typename Tuint>
// typename HllErtl<Tuint>::hval_t
double
HllErtl<Tuint>::estimate_orig() const {
  if(1 >= no_insert()) {
    return no_insert();
  }
  const double m = no_reg();
  const double l2_32 = (double) (((uint64_t) 1) << 32);
  const double lEstRaw = estimate_raw();
  double lRes = 0;
  // std::cout << "HllErtl::estimate_orig: " << lEstRaw << " <?= " << ((5*m)/2) << std::endl;
  // std::cout << "c0 = " << c0 << ", est_raw = " << lEstRaw << std::endl;
  // if(c0 >= (m/2))
  if(lEstRaw <= ((5*m)/2)) {
    const double c0 = no_reg_zero();
    return no_reg() * std::log(no_reg() / c0);
  } else
  if(lEstRaw <= (l2_32 / 30)) {
    lRes = lEstRaw;
  } else {
    lRes = - l2_32 * std::log(1 - (lEstRaw / l2_32));
  }
  return lRes;
}

template<typename Tuint>
// typename HllErtl<Tuint>::hval_t
double
HllErtl<Tuint>::estimate_pois() const {
  const double m = no_reg();
  const uint   q = get_q();
  uint_vt C(q+2);
  extract_counts(C);
  // std::cout << "C = " << C << std::endl;
  double z = m * tau(1 - C[q+1]/m);
  for(uint k = q; 1 <= k; --k) {
    z = 0.5 * (z + C[k]);
  }
  z = z + m * sigma(C[0]/m);
  return alpha_inf() * m * m / z;
}

// TUTS NICHT [to elim or to fix]

template<typename Tuint>
typename HllErtl<Tuint>::hval_t
HllErtl<Tuint>::estimate_maxl() const {
  const double m = no_reg();
  const int    q = get_q();
  std::cout << "L000: p = " << get_p() << ", q = " << get_q() << std::endl; 
  uint_vt C(q+2);
  extract_counts(C);
  if(C[q+1]) {
    return std::numeric_limits<hval_t>::max();
  }
  const int K0_min = min_idx_gt_zero(C);
  const int K1_min = std::max<uint>(1, K0_min);
  const int K0_max = max_idx_gt_zero(C);
  const int K1_max = std::min(q, K0_max);
  std::cout << "L001: R ="; 
  for(uint i = 0; i < no_reg(); ++i) { 
    std::cout << ' ' << ((int) reg(i));
  }
  std::cout << std::endl;
  std::cout << "L002: C = " << C << std::endl;
  std::cout << "L003: K0_min/max = " << K0_min << '/' << K0_max << std::endl;
  double z = 0;
  for(int k = K1_max; K1_min >= k; --k) {
    z = 0.5 * z + C[k];
  }
  z = z / ((double) (1 << K1_min));
  double c = C[q+1];
  if(1 <= q) {
    c = c + C[K1_max];
  }
  double a  = z + C[0];
  double b  = z + C[q+1] / ((double) (1 << q));
  uint   mp = m - C[0];
  double x  = 0;
  if(b <= 1.5 * a) {
    x = mp / (0.5 * b + a);
  } else {
    x = (mp / b) * std::log(1 + b/a);
  }
  const double eps   = 0.01;
  const double delta = eps / std::sqrt(m);
  double xp, xpp, h, g, g_prev = 0;
  int    kappa;
  double x_delta = x;
  uint   lCount = 0;
  std::cout.precision(6);
  while(x_delta > (x * delta)) {
    std::cout << "L100: " << lCount << std::endl;
    std::cout << "  L101: x/log2(X): " << x << '/' << std::log2(x) << std::endl;
    kappa = 2 + std::floor(std::log2(x));
    xp    = x / ((double) (1 << (std::max<uint>(K1_max, kappa) - 1)));
    std::cout << "  L102: xp: " << xp << std::endl;
    xpp   = xp * xp;
    h     = xp - xpp/3 + (xpp * xpp) * (((double) 1.0 / (double) 45.0) - xpp / 472.5);
    std::cout << "  L103: " << (kappa - 1) << ' ' << K1_max << std::endl;
    for(int k = kappa - 1; k <= K1_max; ++k) {
      h  = (xp + h * (1-h)) / (xp + (1 - h));
      xp = 2 * xp;
    }
    g = c * h;
    for(int k = K1_max - 1; k >= K1_min; --k) {
      h  = (xp + h * (1 - h)) / (xp + (1 - h));
      g  = g + C[k] * h;
      xp = 2 * xp;
    }
    g = g + x * a;
    if((g > g_prev) && (mp >= g)) {
      x_delta = x_delta * ((mp - g) / (g - g_prev));
    } else {
      x_delta = 0;
    }
    x = x + x_delta;
    g_prev = g;
    ++lCount;
    if(10 < lCount) {
      break;
    }
  }

  return std::round(m * x);
}

template<typename Tuint>
typename HllErtl<Tuint>::hval_t
HllErtl<Tuint>::estimate_linc() const {
  return 1;
}



template<typename Tuint>
typename HllErtl<Tuint>::hval_t 
HllErtl<Tuint>::no_reg_zero() const {
  uint lRes = 0;
  for(uint i = 0; i < no_reg(); ++i) {
    lRes += (0 == reg(i));
  }
  return lRes;
}

template<typename Tuint>
typename HllErtl<Tuint>::hval_t 
HllErtl<Tuint>::no_reg_nonzero() const {
  uint lRes = 0;
  for(uint i = 0; i < no_reg(); ++i) {
    lRes += (0 != reg(i));
  }
  return lRes;
}

template<typename Tuint>
void
HllErtl<Tuint>::extract_counts(uint_vt& C) const {
  const uint   q = get_q();
  C.resize(q+2);
  for(uint i = 0; i < C.size(); ++i) {
    C[i] = 0;
  }
  for(uint i = 0; i < no_reg(); ++i) {
    ++C[reg(i)];
  }
}


template<typename Tuint>
uint
HllErtl<Tuint>::min_idx_gt_zero(const uint_vt& C) const {
  for(uint i = 0; i < C.size(); ++i) {
    if(0 < C[i]) {
      return i;
    }
  }
  return -1;
}

template<typename Tuint>
uint
HllErtl<Tuint>::max_idx_gt_zero(const uint_vt& C) const {
  for(int i = C.size() - 1; 0 <= i; --i) {
    if(0 < C[i]) {
      return i;
    }
  }
  return -1;
}


template<typename Tuint>
double
HllErtl<Tuint>::estimate_raw() const {
  double lSum = 0;
  for(uint i = 0; i < no_reg(); ++i) {
    lSum += 1 / ((double) ((hval_t) 1 << reg(i)));
  }
  return (1/lSum) * no_reg() * no_reg() * alpha_m(); 
}

template<typename Tuint>
double
HllErtl<Tuint>::alpha_m() const {
  switch(log2_noreg()) {
    case 4: return 0.673;
    case 5: return 0.697;
    default: break;
  }
  return alpha_inf();
}

template<typename Tuint>
double
HllErtl<Tuint>::sigma(double x) const {
  if(1 == x) {
    return std::numeric_limits<double>::infinity();
  }
  double y = 1;
  double z = x;
  double w = 0;
  do {
    x = x * x;
    w = z;
    z = z + x * y;
    y = 2 * y;
  } while(w != z);
  return z; 
}

template<typename Tuint>
double
HllErtl<Tuint>::tau(double x) const {
  if((0 == x) || (1 == x)) {
    return 0;
  }
  double y = 1;
  double z = 1 - x;
  double w = 0;
  do {
    x = std::sqrt(x);
    w = z;
    y = 0.5 * y;
    z = z - (1 - x) * (1 - x) * y;
  } while(w != z);
  return z/3;
}

template<typename Tuint>
void
HllErtl<Tuint>::union_with(const HllErtl& aHll) {
  if(log2_noreg() != aHll.log2_noreg()) {
    std::cout << "HllErtl<Tuint>::union_with:" << std::endl
              << "  log2_noreg()      = " << log2_noreg() << std::endl
              << "  aHll.log2_noreg() = " << aHll.log2_noreg() << std::endl
              << std::endl;
  }
  assert(no_bit() == aHll.no_bit());
  assert(log2_noreg() == aHll.log2_noreg());
  _no_insert += aHll.no_insert();
  for(uint i = 0; i < no_reg(); ++i) {
    if(_regs[i] < aHll.reg(i)) {
      _regs[i] = aHll.reg(i);
    }
  }  
}

// set's this.reg(i) to minimum of this.reg(i) and aHll.reg(i)
template<typename Tuint>
void
HllErtl<Tuint>::minim_with(const HllErtl& aHll) {
  if(log2_noreg() != aHll.log2_noreg()) {
    std::cout << "HllErtl<Tuint>::minim_with:" << std::endl
              << "  log2_noreg()      = " << log2_noreg() << std::endl
              << "  aHll.log2_noreg() = " << aHll.log2_noreg() << std::endl
              << std::endl;
  }
  assert(no_bit() == aHll.no_bit());
  assert(log2_noreg() == aHll.log2_noreg());
  _no_insert = std::min(no_insert(), aHll.no_insert()); // bad estimate
  for(uint i = 0; i < no_reg(); ++i) {
    if(_regs[i] > aHll.reg(i)) {
      _regs[i] = aHll.reg(i);
    }
  }
}

template<typename Tuint>
void
HllErtl<Tuint>::clear() {
  _no_insert = 0;
  for(uint i = 0; i < no_reg(); ++i) {
    _regs[i] = 0;
  }
}

template<typename Tuint>
std::ostream&
HllErtl<Tuint>::print(std::ostream& os) const {
  os << "Hll[" << no_insert() << "] = {";
  for(uint i = 0; i < no_reg(); ++i) {
    os << ' ' << ((uint) reg(i));
  }
  os << " }";
  return os;
}


#endif
