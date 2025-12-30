#ifndef ARBITRARY_PRECISION_FLOATS_TPP
#define ARBITRARY_PRECISION_FLOATS_TPP

template <num T> apf::apf(T x) {
  mpfr_init2(value, Precision);
  if constexpr (std::is_same_v<T, apf>)
    mpfr_set(value, x.value, Rounding);
  else if constexpr (std::floating_point<T>)
    mpfr_set_d(value, static_cast<double>(x), Rounding);
  else
    mpfr_set_si(value, static_cast<long>(x), Rounding);
}

template <num T> apf &apf::operator=(T x) {
  if constexpr (std::is_same_v<T, apf>)
    mpfr_set(value, x.value, Rounding);
  else if constexpr (std::floating_point<T>)
    mpfr_set_d(value, static_cast<double>(x), Rounding);
  else
    mpfr_set_si(value, static_cast<long>(x), Rounding);
  return *this;
}

template <num T> apf apf::operator+(const T &rhs) const {
  apf result;
  if constexpr (std::is_same_v<T, apf>)
    mpfr_add(result.value, value, rhs.value, Rounding);
  else if constexpr (std::floating_point<T>)
    mpfr_add_d(result.value, value, rhs, Rounding);
  else {
    mpz_t tmp;
    mpz_init(tmp);
    mpz_set_si(tmp, static_cast<long>(rhs));
    mpfr_add_z(result.value, value, tmp, Rounding);
    mpz_clear(tmp);
  }
  return result;
}

template <num T> apf apf::operator-(const T &rhs) const {
  apf result;
  if constexpr (std::is_same_v<T, apf>)
    mpfr_sub(result.value, value, rhs.value, Rounding);
  else if constexpr (std::floating_point<T>)
    mpfr_sub_d(result.value, value, rhs, Rounding);
  else {
    mpz_t tmp;
    mpz_init(tmp);
    mpz_set_si(tmp, static_cast<long>(rhs));
    mpfr_sub_z(result.value, value, tmp, Rounding);
    mpz_clear(tmp);
  }
  return result;
}

template <num T> apf apf::operator*(const T &rhs) const {
  apf result;
  if constexpr (std::is_same_v<T, apf>)
    mpfr_mul(result.value, value, rhs.value, Rounding);
  else if constexpr (std::floating_point<T>)
    mpfr_mul_d(result.value, value, rhs, Rounding);
  else {
    mpz_t tmp;
    mpz_init(tmp);
    mpz_set_si(tmp, static_cast<long>(rhs));
    mpfr_mul_z(result.value, value, tmp, Rounding);
    mpz_clear(tmp);
  }
  return result;
}

// returns +/- INF is apf lhs != 0 and rhs = 0
template <num T> apf apf::operator/(const T &rhs) const {
  apf result;
  if constexpr (std::is_same_v<T, apf>)
    mpfr_div(result.value, value, rhs.value, Rounding);
  else if constexpr (std::floating_point<T>)
    mpfr_div_d(result.value, value, rhs, Rounding);
  else {
    mpz_t tmp;
    mpz_init(tmp);
    mpz_set_si(tmp, static_cast<long>(rhs));
    mpfr_div_z(result.value, value, tmp, Rounding);
    mpz_clear(tmp);
  }
  return result;
}

template <num_no_apf T> apf operator+(const T &lhs, const apf &rhs) {
  return apf(lhs) + rhs;
}

template <num_no_apf T> apf operator-(const T &lhs, const apf &rhs) {
  return apf(lhs) - rhs;
}

template <num_no_apf T> apf operator*(const T &lhs, const apf &rhs) {
  return apf(lhs) * rhs;
}

template <num_no_apf T> apf operator/(const T &lhs, const apf &rhs) {
  return apf(lhs) / rhs;
}

template <num T> apf &apf::operator+=(const T &rhs) {
  *this = *this + rhs;
  return *this;
}

template <num T> apf &apf::operator-=(const T &rhs) {
  *this = *this - rhs;
  return *this;
}

template <num T> apf &apf::operator*=(const T &rhs) {
  *this = *this * rhs;
  return *this;
}

template <num T> apf &apf::operator/=(const T &rhs) {
  *this = *this / rhs;
  return *this;
}

template <num T> bool apf::operator<(const T &rhs) const {
  if constexpr (std::is_same_v<T, apf>)
    return mpfr_cmp(value, rhs.value) < 0;
  else if constexpr (std::floating_point<T>)
    return mpfr_cmp_d(value, static_cast<double>(rhs)) < 0;
  else {
    mpz_t tmp;
    mpz_init(tmp);
    mpz_set_si(tmp, static_cast<long>(rhs));
    bool result = mpfr_cmp_z(value, tmp) < 0;
    mpz_clear(tmp);
    return result;
  }
}

template <num T> bool apf::operator>(const T &rhs) const {
  if constexpr (std::is_same_v<T, apf>)
    return mpfr_cmp(value, rhs.value) > 0;
  else if constexpr (std::floating_point<T>)
    return mpfr_cmp_d(value, static_cast<double>(rhs)) > 0;
  else {
    mpz_t tmp;
    mpz_init(tmp);
    mpz_set_si(tmp, static_cast<long>(rhs));
    bool result = mpfr_cmp_z(value, tmp) > 0;
    mpz_clear(tmp);
    return result;
  }
}

template <num T> bool apf::operator<=(const T &rhs) const {
  return !(*this > rhs);
}

template <num T> bool apf::operator>=(const T &rhs) const {
  return !(*this < rhs);
}

template <num T> bool apf::operator==(const T &rhs) const {
  if constexpr (std::is_same_v<T, apf>)
    return mpfr_equal_p(value, rhs.value) != 0;
  else if constexpr (std::floating_point<T>)
    return mpfr_cmp_d(value, static_cast<double>(rhs)) == 0;
  else {
    mpz_t tmp;
    mpz_init(tmp);
    mpz_set_si(tmp, static_cast<long>(rhs));
    bool result = mpfr_cmp_z(value, tmp) == 0;
    mpz_clear(tmp);
    return result;
  }
}

template <num T> bool apf::operator!=(const T &rhs) const {
  return !(*this == rhs);
}

template <num_no_apf T> bool operator<(const T &lhs, const apf &rhs) {
  return apf(lhs) < rhs;
}

template <num_no_apf T> bool operator>(const T &lhs, const apf &rhs) {
  return apf(lhs) > rhs;
}

template <num_no_apf T> bool operator<=(const T &lhs, const apf &rhs) {
  return apf(lhs) <= rhs;
}

template <num_no_apf T> bool operator>=(const T &lhs, const apf &rhs) {
  return apf(lhs) >= rhs;
}

template <num_no_apf T> bool operator==(const T &lhs, const apf &rhs) {
  return apf(lhs) == rhs;
}

template <num_no_apf T> bool operator!=(const T &lhs, const apf &rhs) {
  return apf(lhs) != rhs;
}

template <num T> apf apf::pow(const apf &base, const T &exponent) {
  apf result;
  if constexpr (std::is_same_v<T, apf>)
    mpfr_pow(result.value, base.value, exponent.value, Rounding);
  else if constexpr (std::floating_point<T>) {
    apf aexp = exponent;
    mpfr_pow(result.value, base.value, aexp.value, Rounding);
  } else
    mpfr_pow_si(result.value, base.value, static_cast<long>(exponent),
                Rounding);

  return result;
}

#endif
