#ifndef ARBITRARY_PRECISION_FLOATS_HPP
#define ARBITRARY_PRECISION_FLOATS_HPP

/*
  apf
  Arbitrary Precision Floats
  A wrapper on GNU MPFR

  Kat Matheson

  Usage: compile with "-lmpfr -lgmp" tags for MPFR
*/

#include "mpfr.h"
#include <cmath>
#include <concepts>
#include <functional>
#include <iostream>

class apf;

// Note: integer values larger than long are cast to long in
// constructor/assignment operators
template <class T>
concept num_no_apf =
    (std::integral<T> && !std::is_same_v<T, bool> && !std::is_same_v<T, char> &&
     !std::is_same_v<T, signed char> && !std::is_same_v<T, unsigned char> &&
     !std::is_same_v<T, wchar_t> && !std::is_same_v<T, char16_t> &&
     !std::is_same_v<T, char32_t>) ||
    std::floating_point<T>;

template <class T>
concept num = num_no_apf<T> || std::same_as<T, apf>;

class apf {
public:
  mpfr_t value;
  // Precision (in bits, not digits)
  static const mpfr_prec_t DefaultPrecision = 1000;
  inline static mpfr_prec_t Precision = DefaultPrecision;
  static void ResetPrecision();

  inline static mpfr_rnd_t Rounding = MPFR_RNDN;
  // numbers with significant digits exceeding this value will be printed as
  // []e+/-[]
  inline static int PrintExpThreshold = 10;
  // hashing parameters
  inline static size_t hashMax = 200;
  inline static size_t hashDelta = 20;

  apf();
  apf(const apf &);
  apf(const apf &&) noexcept;
  apf(const std::string &);
  static apf from_string(const std::string &);
  template <num T> apf(T x);
  ~apf();

  // assignment operators
  apf &operator=(const apf &);
  apf &operator=(const apf &&) noexcept;
  template <num T> apf &operator=(T x);

  // arithmetic operators
  apf operator-() const;
  template <num T> apf operator+(const T &rhs) const;
  template <num T> apf operator-(const T &rhs) const;
  template <num T> apf operator*(const T &rhs) const;
  template <num T> apf operator/(const T &rhs) const;

  template <num_no_apf T> friend apf operator+(const T &lhs, const apf &rhs);
  template <num_no_apf T> friend apf operator-(const T &lhs, const apf &rhs);
  template <num_no_apf T> friend apf operator*(const T &lhs, const apf &rhs);
  template <num_no_apf T> friend apf operator/(const T &lhs, const apf &rhs);

  template <num T> apf &operator+=(const T &rhs);
  template <num T> apf &operator-=(const T &rhs);
  template <num T> apf &operator*=(const T &rhs);
  template <num T> apf &operator/=(const T &rhs);

  template <num T> bool operator<(const T &rhs) const;
  template <num T> bool operator>(const T &rhs) const;
  template <num T> bool operator<=(const T &rhs) const;
  template <num T> bool operator>=(const T &rhs) const;
  template <num T> bool operator==(const T &rhs) const;
  template <num T> bool operator!=(const T &rhs) const;

  template <num_no_apf T> friend bool operator<(const T &lhs, const apf &rhs);
  template <num_no_apf T> friend bool operator>(const T &lhs, const apf &rhs);
  template <num_no_apf T> friend bool operator<=(const T &lhs, const apf &rhs);
  template <num_no_apf T> friend bool operator>=(const T &lhs, const apf &rhs);
  template <num_no_apf T> friend bool operator==(const T &lhs, const apf &rhs);
  template <num_no_apf T> friend bool operator!=(const T &lhs, const apf &rhs);

  // functions
  static apf exp(const apf &);
  static apf log(const apf &);
  static apf log10(const apf &);
  static apf cos(const apf &);
  static apf sin(const apf &);
  static apf sqrt(const apf &);
  static apf abs(const apf &);
  static apf erf(const apf &);
  static apf normalCDF(const apf &);
  template <num T> static apf pow(const apf &base, const T &exponent);

  static double trim(const apf &);
  std::string to_string() const;
  friend std::ostream &operator<<(std::ostream &, const apf &);
  static const apf Inf, NegInf;
};

namespace std {
template <> struct hash<apf> {
  size_t operator()(const apf &a) const;
};
}; // namespace std

#include "./apf.tpp"

#endif
