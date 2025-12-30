#include "./apf.hpp"
#include <sstream>

void apf::ResetPrecision() { Precision = DefaultPrecision; }

apf::apf() {
  mpfr_init2(value, Precision);
  mpfr_set_d(value, 0.0, Rounding);
}

apf::apf(const apf &other) {
  mpfr_init2(value, Precision);
  mpfr_set(value, other.value, Rounding);
}

apf::apf(const apf &&other) noexcept {
  mpfr_init2(value, Precision);
  mpfr_set(value, other.value, Rounding);
}

apf::apf(const std::string &s) {
  mpfr_init2(value, Precision);
  mpfr_set_str(value, s.c_str(), 10, Rounding);
}

apf apf::from_string(const std::string &s) {
  apf result;
  mpfr_set_str(result.value, s.c_str(), 10, Rounding);
  return result;
}

apf::~apf() { mpfr_clear(value); }

apf &apf::operator=(const apf &other) {
  if (this != &other)
    mpfr_set(value, other.value, Rounding);
  return *this;
}

// mpfr doesn't have an explicit move assignment
apf &apf::operator=(const apf &&other) noexcept {
  if (this != &other)
    mpfr_set(value, other.value, Rounding);
  return *this;
}

apf apf::operator-() const {
  apf result;
  mpfr_neg(result.value, value, Rounding);
  return result;
}

apf apf::exp(const apf &x) {
  apf result;
  mpfr_exp(result.value, x.value, Rounding);
  return result;
}

apf apf::log(const apf &x) {
  apf result;
  mpfr_log(result.value, x.value, Rounding);
  return result;
}

apf apf::log10(const apf &x) {
  apf result;
  mpfr_log10(result.value, x.value, Rounding);
  return result;
}

apf apf::cos(const apf &x) {
  apf result;
  mpfr_cos(result.value, x.value, Rounding);
  return result;
}

apf apf::sin(const apf &x) {
  apf result;
  mpfr_sin(result.value, x.value, Rounding);
  return result;
}

apf apf::sqrt(const apf &x) {
  apf result;
  mpfr_sqrt(result.value, x.value, Rounding);
  return result;
}

apf apf::abs(const apf &x) {
  apf result;
  mpfr_set(result.value, x.value, Rounding);
  if (result < apf(0))
    mpfr_neg(result.value, result.value, Rounding);
  return result;
}

apf apf::erf(const apf &x) {
  apf result;
  mpfr_erf(result.value, x.value, Rounding);
  return result;
}

apf apf::normalCDF(const apf &x) {
  return (1 + apf::erf(x / apf::sqrt(2))) / 2;
}

double apf::trim(const apf &x) { return mpfr_get_d(x.value, Rounding); }

std::string apf::to_string() const {
  std::stringstream ss;

  if (mpfr_sgn(value) == 0) {
    ss << "0";
    return ss.str();
  }

  mp_exp_t exp;
  // 0 digits = all significant digits
  char *str = mpfr_get_str(nullptr, &exp, 10, 0, value, Rounding);
  std::string s(str);
  free(str);

  bool neg = s[0] == '-';
  if (neg)
    s = s.substr(1);

  exp--; // adjust to match normalized scientific notation

  if (std::abs(exp) > apf::PrintExpThreshold) {
    // scientific notation
    std::string mantissa;
    mantissa += s[0];
    if (s.size() > 1)
      mantissa += "." + s.substr(1);

    if (neg)
      mantissa = "-" + mantissa;

    ss << mantissa << " e" << exp;
  } else {
    // normal decimal: insert decimal point at correct place
    std::string result;
    if (exp >= 0) {
      result = s.substr(0, exp + 1); // integer part
      if (exp + 1 < (mp_exp_t)s.size())
        result += "." + s.substr(exp + 1); // fractional part
    } else {
      result = "0." + std::string(-exp, '0') + s;
    }

    if (neg)
      result = "-" + result;
    ss << result;
  }

  return ss.str();
}

std::ostream &operator<<(std::ostream &os, const apf &a) {
  os << a.to_string();
  return os;
}

const apf apf::Inf = []() -> apf {
  apf x;
  mpfr_set_inf(x.value, 1);
  return x;
}();

const apf apf::NegInf = []() -> apf {
  apf x;
  mpfr_set_inf(x.value, -1);
  return x;
}();

namespace std {
size_t hash<apf>::operator()(const apf &a) const {
  hash<double> dhash;
  size_t _hash = dhash(apf::trim(a));
  mp_exp_t *exp = new mp_exp_t(1);

  // https://machinecognitis.github.io/Math.Mpfr.Native/html/b706c297-237f-4a39-79f5-763763069e04.htm
  // hashMax is number of digits in the string
  char *astr =
      mpfr_get_str(nullptr, exp, 10, apf::hashMax, a.value, apf::Rounding);

  string prehash(astr);
  delete (astr);
  size_t n = prehash.length();
  if (n < apf::hashDelta) {
    return _hash;
  }

  for (size_t i{}; i < n; i += apf::hashDelta) {
    string hashStep = prehash.substr(n - i, apf::hashDelta);
    hash<string> strhash;
    _hash = _hash ^ (strhash(hashStep) << 1);
  }

  return _hash;
}
} // namespace std
