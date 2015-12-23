// copyright yanghao 2013

#pragma once

#include <array>
#include <cassert>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <vector>

namespace rsa {
using std::operator<<;
using std::operator>>;

static const uint64_t LONG_MASK = 0xffffffffL;
static const int DEFAULT_PRIME_CERTAINTY = 100;

namespace details {

template <class _iter>
struct is_uint32_iterator
    : std::is_same<typename std::iterator_traits<_iter>::value_type, uint32_t> {
};

template <class _itera>
_itera find_least_significant_nonzero(_itera abegin, _itera aend) {
  static_assert(is_uint32_iterator<_itera>::value, "invalid iterator type");
  while (abegin != aend) {
    if (*abegin != 0)
      return abegin;
    abegin++;
  }
  return aend;
}

template <class _itera>
_itera find_most_significant_nonzero(_itera abegin, _itera aend) {
  static_assert(is_uint32_iterator<_itera>::value, "invalid iterator type");
  std::reverse_iterator<_itera> arbegin(aend), arend(abegin);
  while (arbegin != arend) {
    if (*arbegin != 0)
      return --arbegin.base();
    arbegin++;
  }
  return aend;
}

template <class _itera, class _iterb>
int compare_range(_itera abegin, _itera aend, _iterb bbegin, _iterb bend) {
  static_assert(is_uint32_iterator<_itera>::value, "invalid iterator type");
  static_assert(is_uint32_iterator<_iterb>::value, "invalid iterator type");
  int result = 0;
  auto i = abegin;
  auto j = bbegin;
  for (; i != aend && j != bend; ++i, ++j) {
    if (*i != *j)
      result = (*i < *j) ? -1 : 1;
  }
  if (i != aend && j == bend) {
    for (; i != aend; ++i) {
      if (*i != 0)
        return 1;
    }
  }
  if (i == aend && j != bend) {
    for (; j != bend; ++j) {
      if (*j != 0)
        return -1;
    }
  }
  return result;
}

// ~a
template <class _itera> void twos_complement_range(_itera begin, _itera end) {
  static_assert(is_uint32_iterator<_itera>::value, "invalid iterator type");
  for (auto i = begin; i != end; ++i)
    *i = ~(*i);
}

// -a
template <class _itera> void negate_range(_itera begin, _itera end) {
  static_assert(is_uint32_iterator<_itera>::value, "invalid iterator type");
  twos_complement_range(begin, end);
  plus_one_range(begin, end);
}

// ++a
template <class _itera> void plus_one_range(_itera begin, _itera end) {
  static_assert(is_uint32_iterator<_itera>::value, "invalid iterator type");
  bool carry = true;
  for (auto i = begin; i != end; ++i) {
    uint64_t sumhere = (LONG_MASK & *i) + (carry ? 1 : 0);
    *i = uint32_t(sumhere);
    carry = (sumhere >> 32) != 0;
    if (!carry)
      break;
  }
}

// --a
template <class _itera> void subtract_one_range(_itera begin, _itera end) {
  static_assert(is_uint32_iterator<_itera>::value, "invalid iterator type");
  bool borrow = true;
  for (auto i = begin; i != end; ++i) {
    int64_t subhere = (LONG_MASK & *i) - (borrow ? 1 : 0);
    *i = uint32_t(subhere);
    borrow = subhere < 0;
    if (!borrow)
      break;
  }
}

// a += b
template <class _itera, class _iterb>
void plus_range(_itera abegin, _itera aend, _iterb bbegin, _iterb bend) {
  static_assert(is_uint32_iterator<_itera>::value, "invalid iterator type");
  static_assert(is_uint32_iterator<_iterb>::value, "invalid iterator type");
  bool carry = false;
  _itera ai;
  _iterb bi;
  for (ai = abegin, bi = bbegin; ai != aend && bi != bend; ++ai, ++bi) {
    uint64_t sumhere = (LONG_MASK & *ai) + (LONG_MASK & *bi) + (carry ? 1 : 0);
    *ai = uint32_t(sumhere);
    carry = (sumhere >> 32) != 0;
  }
  if (ai != aend && carry)
    *ai += 1;
  // assert(!(ai == aend && bi != bend)); // b should be shorter
}

// a -= b
template <class _itera, class _iterb>
void subtract_range(_itera abegin, _itera aend, _iterb bbegin, _iterb bend) {
  static_assert(is_uint32_iterator<_itera>::value, "invalid iterator type");
  static_assert(is_uint32_iterator<_iterb>::value, "invalid iterator type");
  bool borrow = false;
  _itera ai;
  _iterb bi;
  for (ai = abegin, bi = bbegin; ai != aend && bi != bend; ++ai, ++bi) {
    uint64_t along = LONG_MASK & *ai;
    uint64_t blong = LONG_MASK & *bi;
    if (along >= blong + (borrow ? 1 : 0)) {
      *ai = uint32_t(along - blong - (borrow ? 1 : 0));
      borrow = false;
    } else {
      *ai = uint32_t(along + (1UL << 32U) - blong - (borrow ? 1 : 0));
      borrow = true;
    }
  }
  // assert(!(ai == aend && bi != bend)); // b should be shorter
  if (borrow) { // still borrow
    while (ai != aend && *ai == 0) {
      *ai = *ai - 1; // 0 -> 0xffffffff
    }
    if (ai != aend)
      *ai = *ai - 1;
  }
}

// ab = a * b
template <class _itera, class _iterb, class _iterab>
void mult_range(_itera abegin, _itera aend, _iterb bbegin, _iterb bend,
                _iterab abbegin, _iterab abend) {
  static_assert(is_uint32_iterator<_itera>::value, "invalid iterator type");
  static_assert(is_uint32_iterator<_iterb>::value, "invalid iterator type");
  static_assert(is_uint32_iterator<_iterab>::value, "invalid iterator type");
  std::fill(abbegin, abend, 0u);

  // auto bend = bbegin + (aend - abegin);
  _iterab abihere = abbegin;
  for (auto bi = bbegin; bi != bend; ++bi, ++abihere) {
    uint64_t carry = 0;
    _itera ai = abegin;
    _iterab abi = abihere;
    for (; abi != abend && ai != aend; ++ai, ++abi) {
      uint64_t product =
          (*ai & LONG_MASK) * (*bi & LONG_MASK) + (*abi & LONG_MASK) + carry;
      *abi = (uint32_t)product;
      carry = product >> 32;
    }
    if (abi != abend)
      *abi = uint32_t(carry);
  }
}

// a <<= n
template <class _itera>
void left_shift_range(_itera abegin, _itera aend, uint32_t n) {
  static_assert(is_uint32_iterator<_itera>::value, "invalid iterator type");
  assert(n < 32);
  uint32_t carry = 0;
  for (auto i = abegin; i != aend; ++i) {
    uint64_t s = ((*i & LONG_MASK) << n) + carry;
    *i = uint32_t(s);
    carry = (s >> 32);
  }
}

// a >>= n
template <class _itera>
void right_shift_range(_itera abegin, _itera aend, uint32_t n) {
  static_assert(is_uint32_iterator<_itera>::value, "invalid iterator type");
  assert(n < 32);
  std::reverse_iterator<_itera> arbegin(aend), arend(abegin);
  uint32_t carry = 0;
  for (auto i = arbegin; i != arend; ++i) {
    uint64_t s = ((*i & LONG_MASK) + ((carry & LONG_MASK) << 32));
    *i = uint32_t(s >> n);
    carry = uint32_t(s - ((*i & LONG_MASK) << n));
  }
}

// a & (1 << n)
template <class _itera>
bool test_bit_range(_itera abegin, _itera aend, uint32_t n) {
  static_assert(is_uint32_iterator<_itera>::value, "invalid iterator type");
  uint32_t unitid = n / 32;
  uint32_t offset = n % 32;
  return (*(abegin + unitid) & (1ull << offset)) != 0;
}

// number of bits
inline int number_of_bits(uint32_t i) {
  i = i - ((i >> 1) & 0x55555555);
  i = (i & 0x33333333) + ((i >> 2) & 0x33333333);
  return (((i + (i >> 4)) & 0x0F0F0F0F) * 0x01010101) >> 24;
}

template <class _itera> int number_of_bits_range(_itera abegin, _itera aend) {
  static_assert(is_uint32_iterator<_itera>::value, "invalid iterator type");
  int n = 0;
  for (auto i = abegin; i != aend; ++i) {
    n += number_of_bits(*i);
  }
  return n;
}

// quotient = dividend / divisor, didivend -= quotient * divisor (remainder)
template <class _iterDividend, class _iterDivisor, class _iterQuotient>
void divide_and_remain_range(_iterDividend ddbegin, _iterDividend ddend,
                             _iterDivisor dsbegin, _iterDivisor dsend,
                             _iterQuotient qbegin, _iterQuotient qend) {
  static_assert(is_uint32_iterator<_iterDividend>::value,
                "invalid iterator type");
  static_assert(is_uint32_iterator<_iterDivisor>::value,
                "invalid iterator type");
  static_assert(is_uint32_iterator<_iterQuotient>::value,
                "invalid iterator type");

  int cmp = compare_range(ddbegin, ddend, dsbegin, dsend);
  if (cmp < 0) { // dd < ds
    std::fill(qbegin, qend, 0u);
    return;
  }

  _iterDivisor dsbig = find_most_significant_nonzero(dsbegin, dsend);
  _iterDividend ddbig = find_most_significant_nonzero(ddbegin, ddend);
  if (dsbig == dsend) { // ds = 0
    throw std::runtime_error("divided by zero");
  }
  if (ddbig == ddend) { // dd = 0
    std::fill(qbegin, qend, 0u);
    return;
  }
  uint64_t dshead = *dsbig & LONG_MASK;

  typedef std::reverse_iterator<_iterDividend> _riterDividend;

  // 0 ... [0? ddtoclear ... ddcurrenton] ddcurrenton-1 ddcurrenton-2 ...
  // ddbegin
  _iterDividend ddcurrenton = ddbig - (dsbig - dsbegin);
  _iterDividend ddtoclear = ddbig; // will be cleared
  _riterDividend ddcurrenton_rev(ddcurrenton + 1);
  _riterDividend ddtoclear_rev(ddtoclear + 1);
  _riterDividend ddrevend(ddbegin);

  std::vector<uint32_t> dsmultgess(dsend - dsbegin + 1,
                                   0u); // cache for multiplication
  while (ddcurrenton_rev != ddrevend) {
    ddcurrenton = ddcurrenton_rev.base() - 1;
    ddtoclear = ddtoclear_rev.base() - 1;

    int cmp = compare_range(ddcurrenton, ddend, dsbegin, dsend);
    if (cmp >= 0) {
      assert(ddtoclear - ddcurrenton == dsbig - dsbegin);
      uint64_t ddtoclearhead = *ddtoclear & LONG_MASK;
      if (ddtoclear + 1 != ddend) {
        ddtoclearhead += ((*(ddtoclear + 1) & LONG_MASK) << 32);
      }
      assert(ddtoclearhead >= dshead);
      uint64_t qgess =
          (dshead + 1ul == 0) ? 1ul : ddtoclearhead / (dshead + 1ul);
      qgess = std::max<uint64_t>(1ul, qgess);
      assert(qgess >= 1ul && uint32_t(qgess >> 32) == 0);
      if (qgess > 1) {
        uint32_t qgess_[] = {uint32_t(qgess)};
        std::fill(dsmultgess.begin(), dsmultgess.end(), 0u);
        mult_range(dsbegin, dsend, qgess_, qgess_ + 1, dsmultgess.begin(),
                   dsmultgess.end());
        subtract_range(ddcurrenton, ddend, dsmultgess.begin(),
                       dsmultgess.end());
        *(qbegin + (ddcurrenton - ddbegin)) += uint32_t(qgess);
      } else {
        subtract_range(ddcurrenton, ddend, dsbegin, dsbig + 1);
        *(qbegin + (ddcurrenton - ddbegin)) += 1;
      }
    } else {
      ddcurrenton_rev++;
      ddtoclear_rev++;
    }
  }
}
}

// fixed_uint
template <int num> class fixed_uint {
  static_assert(num > 1, "num should be over 2");

public:
  typedef uint32_t unit_t;
  typedef uint64_t duint_t;
  static const int unit_num = num;

  inline fixed_uint() { std::fill(mag_.begin(), mag_.end(), 0); }
  fixed_uint(uint64_t val) {
    std::fill(mag_.begin(), mag_.end(), 0);
    uint32_t highword = (uint32_t)(uint64_t(val) >> 32);
    if (highword == 0) {
      mag_[0] = uint32_t(val);
    } else {
      mag_[1] = highword;
      mag_[0] = uint32_t(val);
    }
  }
  explicit fixed_uint(const std::array<uint32_t, unit_num> &val) : mag_(val) {}
  fixed_uint(const fixed_uint &v) : mag_(v.mag_) {}
  template <int num2> explicit fixed_uint(const fixed_uint<num2> &v) {
    std::fill(mag_.begin(), mag_.end(), 0);
    if (num2 <= num)
      std::copy(v.mag_.begin(), v.mag_.end(), mag_.begin());
    else
      std::copy(v.mag_.begin(), v.mag_.begin() + num, mag_.begin());
  }
  template <class _it> explicit fixed_uint(_it begin, _it end) {
    std::fill(mag_.begin(), mag_.end(), 0);
    std::copy(begin, end, mag_.begin());
  }
  fixed_uint(const std::string &s, int radix = 10) {
    std::fill(mag_.begin(), mag_.end(), 0);
    from_string(s, radix);
  }

  inline const std::array<uint32_t, unit_num> &magnitude() const {
    return mag_;
  }
  inline uint32_t mag_at(int i) const { return mag_[i]; }
  inline uint32_t &mag_at(int i) { return mag_[i]; }
  inline uint32_t operator[](int i) const { return mag_[i]; }
  inline uint32_t &operator[](int i) { return mag_[i]; }

  typedef uint32_t value_type;
  typedef typename std::array<uint32_t, num>::iterator iterator;
  typedef typename std::array<uint32_t, num>::const_iterator const_iterator;

  inline iterator begin() { return mag_.begin(); }
  inline iterator end() { return mag_.end(); }
  inline const_iterator begin() const { return mag_.begin(); }
  inline const_iterator end() const { return mag_.end(); }
  inline const_iterator cbegin() const { return mag_.cbegin(); }
  inline const_iterator cend() const { return mag_.cend(); }

  inline uint32_t to_uint32() const { return mag_[0]; }
  inline uint64_t to_uint64() const {
    return (mag_[0] & LONG_MASK) + ((mag_[1] & LONG_MASK) << 32);
  }
  std::string to_string(int radix = 10) const;
  void from_string(const std::string &s, int radix = 10);

  template <class randgen_t> void randomize(randgen_t &rnd);
  template <class randgen_t>
  void set_to_probable_prime(randgen_t &rnd, int maxbitoffset = 0,
                             int certainty = DEFAULT_PRIME_CERTAINTY);

public:
  // comparison
  int compare_to(const fixed_uint &v) const {
    return details::compare_range(mag_.begin(), mag_.end(), v.mag_.begin(),
                                  v.mag_.end());
  }
  inline bool operator==(const fixed_uint &v) const {
    return compare_to(v) == 0;
  }
  inline bool operator!=(const fixed_uint &v) const {
    return compare_to(v) != 0;
  }
  inline bool operator<(const fixed_uint &v) const {
    return compare_to(v) == -1;
  }
  inline bool operator>(const fixed_uint &v) const {
    return compare_to(v) == 1;
  }
  inline bool operator<=(const fixed_uint &v) const {
    return compare_to(v) <= 0;
  }
  inline bool operator>=(const fixed_uint &v) const {
    return compare_to(v) >= 0;
  }

  bool is_zero() const;
  bool is_one() const;
  inline bool is_odd() const { return test_bit(0); }
  inline bool is_even() const { return !test_bit(0); }
  inline bool is_power_of_2() const { return number_of_bits() == 1; }

  int most_significant_nonzero_unit() const;
  int most_significant_nonzero_bit() const;
  int least_significant_nonzero_unit() const;
  int least_significant_nonzero_bit() const;
  inline int number_of_bits() const {
    return details::number_of_bits_range(mag_.begin(), mag_.end());
  }

  void left_shift(uint32_t n);
  void right_shift(uint32_t n);
  inline bool test_bit(uint32_t n) const {
    return details::test_bit_range(mag_.begin(), mag_.end(), n);
  }

  inline fixed_uint left_shifted(uint32_t n) const {
    fixed_uint s(*this);
    s.left_shift(n);
    return s;
  }
  inline fixed_uint right_shifted(uint32_t n) const {
    fixed_uint s(*this);
    s.right_shift(n);
    return s;
  }

  inline void twos_complement() {
    details::twos_complement_range(mag_.begin(), mag_.end());
  }
  void negate();
  inline void plus_one() { details::plus_one_range(mag_.begin(), mag_.end()); }
  inline void subtract_one() {
    details::subtract_one_range(mag_.begin(), mag_.end());
  }

  fixed_uint operator++(int) {
    fixed_uint o(*this);
    plus_one();
    return o;
  }
  fixed_uint &operator++() {
    plus_one();
    return *this;
  }
  fixed_uint operator--(int) {
    fixed_uint o(*this);
    subtract_one();
    return o;
  }
  fixed_uint &operator--() {
    subtract_one();
    return *this;
  }

  fixed_uint &operator+=(const fixed_uint &v);
  fixed_uint &operator-=(const fixed_uint &v);
  fixed_uint &operator*=(const fixed_uint &v);
  fixed_uint &operator/=(const fixed_uint &v);

  // this = remainder, returns the quotient
  fixed_uint divide_and_remain(const fixed_uint &ds);
  fixed_uint mod(const fixed_uint &ds) const;

  fixed_uint mod_plus(const fixed_uint &b, const fixed_uint &m) const;
  fixed_uint mod_negate(const fixed_uint &m) const { return (m - mod(m)); }
  fixed_uint mod_subtract(const fixed_uint &b, const fixed_uint &m) const {
    return mod_plus(b.mod_negate(m), m);
  }
  fixed_uint mod_mult(const fixed_uint &b, const fixed_uint &m) const;
  fixed_uint mod_inverse(const fixed_uint &b) const;
  fixed_uint mod_pow(const fixed_uint &exp, const fixed_uint &m) const;

  fixed_uint gcd(const fixed_uint &b) const;
  fixed_uint binary_gcd(const fixed_uint &b) const;
  bool is_coprime_with(const fixed_uint &b) const { return gcd(b).is_one(); }

  template <class randgen_t>
  bool passes_miller_rabin_test(int iterations, randgen_t &rnd) const;

  template <class randgen_t>
  bool is_prime_with_certainty(randgen_t &rnd, int certainty) const;

private:
  std::array<uint32_t, unit_num> mag_;

  template <int num> friend class fixed_uint;
};

/************************************************************************/
/* fixed_uint implementation                                    */
/************************************************************************/

template <int num>
std::string fixed_uint<num>::to_string(int radix /*= 10*/) const {
  if (radix <= 1)
    throw std::runtime_error("invalid radix");
  if (is_zero())
    return "0";
  fixed_uint base = radix;
  fixed_uint self(*this);
  std::string s;
  while (!self.is_zero()) {
    int b = self.mod(base).to_uint32();
    s = char(b < 10 ? (char(b) + '0') : (char(b - 10) + 'A')) + s;
    self = self / base;
  }
  return s;
}

template <int num>
void fixed_uint<num>::from_string(const std::string &s, int radix /*= 10*/) {
  if (radix <= 1)
    throw std::runtime_error("invalid radix");
  *this = fixed_uint();
  fixed_uint base = radix;
  for (int i = 0; i < s.size(); i++) {
    assert(isdigit(s[i]) || isalpha(s[i]));
    fixed_uint b =
        isdigit(s[i]) ? (s[i] - '0')
                      : (islower(s[i]) ? (s[i] - 'a' + 10) : (s[i] - 'A' + 10));
    *this = *this * base + b;
  }
}

template <int num>
template <class randgen_t>
void fixed_uint<num>::randomize(randgen_t &rnd) {
  std::uniform_int_distribution<> dist;
  for (int i = 0; i < num; i++) {
    mag_[i] = dist(rnd);
  }
}

static const uint64_t SMALL_PRIME_PRODUCT(3ull * 5 * 7 * 11 * 13 * 17 * 19 *
                                          23 * 29 * 31 * 37 * 41);

template <int num>
template <class randgen_t>
void fixed_uint<num>::set_to_probable_prime(
    randgen_t &rnd, int maxbitoffset,
    int certainty /*= DEFAULT_PRIME_CERTAINTY*/) {
  while (true) {
    // candidate
    randomize<randgen_t>(rnd);
    mag_.back() |= (1u << 31); // make the bit full
    right_shift(maxbitoffset); // offset to lower it
    mag_[0] |= 1;              // make it odd

    uint64_t r = mod(SMALL_PRIME_PRODUCT).to_uint64();
    if ((r % 3 == 0) || (r % 5 == 0) || (r % 7 == 0) || (r % 11 == 0) ||
        (r % 13 == 0) || (r % 17 == 0) || (r % 19 == 0) || (r % 23 == 0) ||
        (r % 29 == 0) || (r % 31 == 0) || (r % 37 == 0) || (r % 41 == 0))
      continue; // Candidate is composite; try another
    if (is_prime_with_certainty<randgen_t>(rnd, certainty))
      break;
  }
}

template <int num> bool fixed_uint<num>::is_zero() const {
  for (int i = 0; i < num; i++)
    if (mag_[i] != 0)
      return false;
  return true;
}

template <int num> bool fixed_uint<num>::is_one() const {
  if (mag_[0] != 1)
    return false;
  for (int i = 1; i < num; i++)
    if (mag_[i] != 0)
      return false;
  return true;
}

template <int num> int fixed_uint<num>::most_significant_nonzero_unit() const {
  auto i = details::find_most_significant_nonzero(mag_.begin(), mag_.end());
  return i == mag_.end() ? -1 : (i - mag_.begin());
}

template <int num> int fixed_uint<num>::most_significant_nonzero_bit() const {
  int msn = most_significant_nonzero_unit();
  if (msn == -1)
    return -1;
  uint32_t b = mag_[msn];
  int a = -1;
  for (; b != 0; a++)
    b >>= 1;
  assert(a >= 0);
  return msn * 32 + a;
}

template <int num> int fixed_uint<num>::least_significant_nonzero_unit() const {
  auto i = details::find_least_significant_nonzero(mag_.begin(), mag_.end());
  return i == mag_.end() ? -1 : (i - mag_.begin());
}

template <int num> int fixed_uint<num>::least_significant_nonzero_bit() const {
  int lsn = least_significant_nonzero_unit();
  if (lsn == -1)
    return -1;
  uint32_t b = mag_[lsn];
  int a = 32;
  for (; b != 0; a--)
    b <<= 1;
  assert(a >= 0);
  return lsn * 32 + a;
}

template <int num> void fixed_uint<num>::left_shift(uint32_t n) {
  while (n != 0) {
    uint32_t nn = std::min(32u, n);
    details::left_shift_range(mag_.begin(), mag_.end(), nn);
    n -= nn;
  }
}

template <int num> void fixed_uint<num>::right_shift(uint32_t n) {
  while (n != 0) {
    uint32_t nn = std::min(32u, n);
    details::right_shift_range(mag_.begin(), mag_.end(), nn);
    n -= nn;
  }
}

template <int num> void fixed_uint<num>::negate() {
  twos_complement();
  plus_one();
}

template <int num>
fixed_uint<num> &fixed_uint<num>::operator+=(const fixed_uint &v) {
  details::plus_range(mag_.begin(), mag_.end(), v.mag_.begin(), v.mag_.end());
  return *this;
}

template <int num>
fixed_uint<num> &fixed_uint<num>::operator-=(const fixed_uint &v) {
  details::subtract_range(mag_.begin(), mag_.end(), v.mag_.begin(),
                          v.mag_.end());
  return *this;
}

template <int num>
fixed_uint<num> &fixed_uint<num>::operator*=(const fixed_uint &v) {
  *this = (*this * v);
  return *this;
}

template <int num>
fixed_uint<num> &fixed_uint<num>::operator/=(const fixed_uint &v) {
  *this = (*this / v);
  return *this;
}

template <int num>
fixed_uint<num> fixed_uint<num>::divide_and_remain(const fixed_uint &ds) {
  fixed_uint<num> q;
  details::divide_and_remain_range(mag_.begin(), mag_.end(), ds.mag_.begin(),
                                   ds.mag_.end(), q.mag_.begin(), q.mag_.end());
  return q;
}

template <int num>
fixed_uint<num> fixed_uint<num>::mod(const fixed_uint &ds) const {
  fixed_uint<num> n(*this), q;
  details::divide_and_remain_range(n.mag_.begin(), n.mag_.end(),
                                   ds.mag_.begin(), ds.mag_.end(),
                                   q.mag_.begin(), q.mag_.end());
  return n;
}

template <int num>
fixed_uint<num> fixed_uint<num>::mod_plus(const fixed_uint &b,
                                          const fixed_uint &m) const {
  if (m.is_zero()) {
    throw std::runtime_error("m is zero!");
  }
  if (m.is_one()) {
    return 0;
  }
  if (b.is_zero())
    return mod(m);

  fixed_uint<num + 1> self(*this);
  fixed_uint<num + 1> temp;
  details::plus_range(self.mag_.begin(), self.mag_.end(), b.mag_.begin(),
                      b.mag_.end());
  details::divide_and_remain_range(
      self.mag_.begin(), self.mag_.end(), m.mag_.begin(), m.mag_.end(),
      temp.mag_.begin(), temp.mag_.end()); // self %= m
  return fixed_uint(self);
}

template <int num>
fixed_uint<num> fixed_uint<num>::mod_mult(const fixed_uint &b,
                                          const fixed_uint &m) const {
  if (m.is_zero())
    throw std::runtime_error("m is zero!");
  if (m.is_one())
    return 0;
  if (b.is_zero())
    return 0;

  fixed_uint<num * 2> prod = mult_unlimit(*this, b), temp;
  details::divide_and_remain_range(
      prod.mag_.begin(), prod.mag_.end(), m.mag_.begin(), m.mag_.end(),
      temp.mag_.begin(), temp.mag_.end()); // prod %= m
  assert(prod < fixed_uint<num * 2>(m));
  return fixed_uint(prod);
}

template <int num>
fixed_uint<num> fixed_uint<num>::mod_inverse(const fixed_uint &m) const {
  if (m.is_zero())
    throw std::runtime_error("modulus is zero!");
  if (m.is_one())
    return 0;

  typedef fixed_uint<2 * num> fixed_unsigned_int2;

  fixed_unsigned_int2 mm(m);
  fixed_unsigned_int2 s = 0, old_s = 1;
  fixed_unsigned_int2 t = 1, old_t = 0;
  fixed_unsigned_int2 r(m), old_r(*this);
  fixed_unsigned_int2 q;
  fixed_unsigned_int2 tmp;

  while (!r.is_zero()) {
    q = old_r / r;
    tmp = old_r;
    old_r = r;
    r = tmp.mod_subtract(q * r, mm);
    tmp = old_s;
    old_s = s;
    s = tmp.mod_subtract(q * s, mm);
    tmp = old_t;
    old_t = t;
    t = tmp.mod_subtract(q * t, mm);
  }
  return fixed_uint(old_s).mod(m);
}

template <int num>
fixed_uint<num> fixed_uint<num>::mod_pow(const fixed_uint &exp,
                                         const fixed_uint &md) const {
  if (exp.is_zero())
    return 1ull;
  if (exp.is_one())
    return mod(md);
  if (is_zero())
    return 0;

  int msb = exp.most_significant_nonzero_bit();
  assert(msb >= 0);
  fixed_uint<2 * num> prod(1ull);
  fixed_uint<2 * num> self(mod(md));
  fixed_uint<2 * num> mmd(md);
  for (int i = msb; i >= 0; --i) {
    prod = prod * prod;
    prod = prod.mod(mmd);
    if (exp.test_bit(i)) {
      prod *= self;
      prod = prod.mod(mmd);
    }
  }

  return fixed_uint(prod);
}

template <int num>
fixed_uint<num> fixed_uint<num>::gcd(const fixed_uint &b) const {
  fixed_uint<num> aa(*this), bb(b);
  if (aa.is_zero())
    return bb;
  if (bb.is_zero())
    return aa;
  if (aa.is_one() || bb.is_one())
    return 1;

  while (!bb.is_zero()) {
    auto tt = bb;
    bb = aa.mod(bb);
    aa = tt;
  }
  return aa;
}

template <int num>
fixed_uint<num> fixed_uint<num>::binary_gcd(const fixed_uint &b) const {
  fixed_uint aa(*this), bb(b);
  int i = 0;
  if (aa.is_zero())
    return bb;
  if (bb.is_zero())
    return aa;
  if (aa.is_one() || bb.is_one())
    return 1;

  int common_low_zero_nums = std::min(aa.least_significant_nonzero_bit(),
                                      bb.least_significant_nonzero_bit());
  aa.right_shift(common_low_zero_nums);
  bb.right_shift(common_low_zero_nums);

  aa.right_shift(aa.least_significant_nonzero_bit());
  bb.right_shift(bb.least_significant_nonzero_bit());
  int cmp;
  while ((cmp = aa.compare_to(bb)) != 0) {
    fixed_uint temp = aa;
    if (cmp < 0) { // aa < bb
      aa = bb - temp;
      bb = temp;
    } else { // bb < aa
      aa = temp - bb;
    }
    aa.right_shift(aa.least_significant_nonzero_bit());
  }
  return aa.left_shifted(common_low_zero_nums);
}

template <int num>
template <class randgen_t>
bool fixed_uint<num>::passes_miller_rabin_test(int iterations,
                                               randgen_t &rnd) const {
  if (is_zero()) {
    return false;
  }
  fixed_uint this_1 = *this;
  --this_1;
  fixed_uint m = this_1;
  if (m.is_zero()) {
    return false;
  }

  int a = m.least_significant_nonzero_bit();
  assert(a != -1);
  m.right_shift(a);

  for (int i = 0; i < iterations; i++) {
    fixed_uint b;
    do {
      b.randomize(rnd);
    } while (b <= 1 || b >= *this);

    int j = 0;
    fixed_uint z = b.mod_pow(m, *this);
    while (!((j == 0 && z.is_one()) || z == this_1)) {
      if (j > 0 && z.is_one() || ++j == a)
        return false;
      z = z.mod_pow(2, *this);
    }
  }
  return true;
}

template <int num>
template <class randgen_t>
bool fixed_uint<num>::is_prime_with_certainty(randgen_t &rnd,
                                              int certainty) const {
  int rounds = 0;
  int n = (std::min(certainty, std::numeric_limits<int>::max() - 1) + 1) / 2;
  int sizeInBits = most_significant_nonzero_bit() + 1;
  if (sizeInBits == 0) // zero
    return false;

  if (sizeInBits < 256) {
    rounds = 27;
  } else if (sizeInBits < 512) {
    rounds = 15;
  } else if (sizeInBits < 768) {
    rounds = 8;
  } else if (sizeInBits < 1024) {
    rounds = 4;
  } else {
    rounds = 2;
  }
  rounds = n < rounds ? n : rounds;
  return passes_miller_rabin_test<randgen_t>(rounds, rnd);
}

/************************************************************************/
/* global functions implementation                                      */
/************************************************************************/

template <int num> fixed_uint<num> operator~(const fixed_uint<num> &v) {
  fixed_uint<num> s(v);
  s.twos_complement();
  return s;
}

template <int num> fixed_uint<num> operator-(const fixed_uint<num> &v) {
  fixed_uint<num> s(v);
  s.negate();
  return s;
}

template <int num>
fixed_uint<num> operator+(const fixed_uint<num> &a, const fixed_uint<num> &b) {
  fixed_uint<num> s(a);
  s += b;
  return s;
}

template <int num>
fixed_uint<num> operator-(const fixed_uint<num> &a, const fixed_uint<num> &b) {
  fixed_uint<num> s(a);
  s -= b;
  return s;
}

// unwrapped
template <int num1, int num2>
fixed_uint<num1 + num2> mult_unlimit(const fixed_uint<num1> &xx,
                                     const fixed_uint<num2> &yy) {
  std::array<uint32_t, num1 + num2> p;
  std::fill(p.begin(), p.end(), 0);
  details::mult_range(xx.magnitude().begin(), xx.magnitude().end(),
                      yy.magnitude().begin(), yy.magnitude().end(), p.begin(),
                      p.end());
  return fixed_uint<num1 + num2>(p);
}

// wrapped
template <int num>
fixed_uint<num> operator*(const fixed_uint<num> &xx,
                          const fixed_uint<num> &yy) {
  std::array<uint32_t, num> p;
  std::fill(p.begin(), p.end(), 0);
  details::mult_range(xx.magnitude().begin(), xx.magnitude().end(),
                      yy.magnitude().begin(), yy.magnitude().end(), p.begin(),
                      p.end());
  return fixed_uint<num>(p);
}

template <int num>
fixed_uint<num> operator/(const fixed_uint<num> &xx,
                          const fixed_uint<num> &yy) {
  fixed_uint<num> xxx(xx);
  return xxx.divide_and_remain(yy);
}

// streamer
template <int num>
std::ostream &operator<<(std::ostream &s, const fixed_uint<num> &v) {
  if (s.flags() & s.hex) {
    s << "[hex " << v.to_string(16) << "]";
  } else if (s.flags() & s.dec) {
    s << "[" << v.to_string(10) << "]";
  } else if (s.flags() & s.oct) {
    s << "[oct " << v.to_string(8) << "]";
  }
  return s;
}

template <int num>
std::istream &operator>>(std::istream &s, fixed_uint<num> &v) {
  std::string txt;
  s >> txt;
  if (s.flags() & s.hex) {
    v.from_string(txt, 16);
  } else if (s.flags() & s.dec) {
    v.from_string(txt, 10);
  } else if (s.flags() & s.oct) {
    v.from_string(txt, 8);
  }
  return s;
}
}
