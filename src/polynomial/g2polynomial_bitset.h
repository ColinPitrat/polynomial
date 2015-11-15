#pragma once

/*
 * Polynomials on GF(2) implemented using a bitset.
 * This is slow on large sparse polynomials.
 */
#include <iostream>
#include <vector>
#include <utility>
#include <algorithm>
#include <cmath>
#include <bitset>

// If changing MAX_SIZE then derivationMask must only be changed
//#define MAX_SIZE 512
#define MAX_SIZE 1024
//#define MAX_SIZE 131072

static std::bitset<MAX_SIZE> derivationMask("1010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010", 1024, '1', '0');

class G2Poly {
  public:
    G2Poly();
    G2Poly(std::vector<uint64_t> coeffs);

    int64_t degree() const;
    bool null() const;

    void simplify(uint64_t maxdeg = MAX_SIZE-1);
    G2Poly power(uint64_t n) const;
    G2Poly unpower(uint64_t n) const;
    G2Poly derivate() const;
    G2Poly squareFreePart(uint64_t p) const;
    std::vector<std::pair<G2Poly, uint64_t>> squareFreeFactors(uint64_t p) const;
    std::vector<G2Poly> distinctDegreeFactors(uint64_t p) const;
    G2Poly cantorZassenhaus(uint64_t p) const;
    G2Poly knuth(uint64_t p) const;
    G2Poly mceliece(uint64_t p) const;
    G2Poly equalDegreeFactorize(uint64_t p, uint64_t i) const;

    G2Poly& operator-=(const G2Poly &p);
    G2Poly& operator+=(const G2Poly &p);
    G2Poly& operator*=(const G2Poly &p);

    static G2Poly Xn(uint64_t n);
    static G2Poly Rand(uint64_t n);

  // TODO: make this protected by removing all direct usage in non friend non member methods
  //protected:
    int64_t degree_;
		std::bitset<MAX_SIZE> coeffs_;
};

bool operator==(const G2Poly &a, const G2Poly &b) {
  return a.coeffs_ == b.coeffs_;
}

bool operator!=(const G2Poly &a, const G2Poly &b) {
  return a.coeffs_ != b.coeffs_;
}

bool operator<(const G2Poly &a, const G2Poly &b) {
  if(a.degree() < b.degree()) return true;
  if(a.degree() > b.degree()) return false;
  for(uint64_t i = a.degree(); i >= 0; i--) {
    if (b.coeffs_[i] && !a.coeffs_[i]) return true;
    if (!b.coeffs_[i] && a.coeffs_[i]) return false;
  }
  return false;
}

std::ostream& operator<<(std::ostream &s, const G2Poly &p) {
  bool first = true;
  if(p.null()) s << "0";
  for(int i = p.degree(); i >= 0; i--) {
    if(p.coeffs_[i]) {
      if (!first) {
        s << " + ";
      }
      first = false;
      if(i > 0) {
        s << "X";
        if(i > 1) {
          s << "^" << i;
        }
      } else {
        s << (p.coeffs_[i] ? "1" : "0");
      }
    }
  }
	return s;
}

G2Poly operator-(const G2Poly &a, const G2Poly &b) {
  auto r = a;
  r -= b;
  return r;
}

G2Poly operator*(const G2Poly &a, const G2Poly &b) {
  auto r = a;
  r *= b;
  return r;
}

G2Poly operator+(const G2Poly &a, const G2Poly &b) {
  auto r = a;
  r += b;
  return r;
}

std::pair<G2Poly, G2Poly> euclidDivide(const G2Poly &a, const G2Poly &b) {
  //std::cerr << "euclidDivide(\n   - " << a << "\n   - " << b << "\n) = ";
  G2Poly r(a);
  G2Poly q;
  while(r.degree() >= b.degree()) {
    int n = r.degree() - b.degree();
    auto bn = b;
    bn.coeffs_ <<= n;
    r -= bn;
    q += G2Poly::Xn(n);
  }
  r.simplify(std::max(b.degree(), 0ll));
  q.simplify(std::max(a.degree() - b.degree(), 0ll));
  //std::cerr << "(" << q << ", " << r << ")" << std::endl;
  return std::make_pair(q, r);
}

G2Poly operator%(const G2Poly &a, const G2Poly &b) {
  return euclidDivide(a,b).second;
}

G2Poly operator/(const G2Poly &a, const G2Poly &b) {
  return euclidDivide(a,b).first;
}

G2Poly::G2Poly() {
  degree_ = -1;
}

/*
G2Poly::G2Poly(std::vector<int> coeffs) {
  coeffs_ = coeffs;
}
*/

G2Poly::G2Poly(std::vector<uint64_t> coeffs) {
  size_t i = 0;
  for(auto it = coeffs.begin(); it != coeffs.end(); it++) {
    coeffs_[i] = (*it == 1);
    i++;
  }
  simplify(coeffs.size());
}

int64_t G2Poly::degree() const {
  return degree_;
}

G2Poly G2Poly::Xn(uint64_t n) {
  G2Poly p;
  p.coeffs_[n] = true;
  p.degree_ = n;
  return p;
}

G2Poly G2Poly::Rand(uint64_t n) {
  G2Poly p;
  for(uint64_t i = 0; i < n; i++) {
    if(rand() & 1) {
      p.coeffs_[i] = true;
      p.degree_ = i;
    }
  }
  return p;
}

void G2Poly::simplify(uint64_t maxdeg) {
  degree_ = -1;
	for(uint64_t i = maxdeg+1; i > 0; i--) {
    if(coeffs_[i-1]) {
      degree_ = i-1;
      break;
    }
  }
}

G2Poly& G2Poly::operator-=(const G2Poly &other) {
  return operator+=(other);
}

G2Poly& G2Poly::operator+=(const G2Poly &other) {
  coeffs_ ^= other.coeffs_;
  simplify(std::max(degree_, other.degree_));
  return *this;
}

G2Poly& G2Poly::operator*=(const G2Poly &other) {
  std::bitset<MAX_SIZE> result;
  for(int i = 0; i <= degree(); i++) {
    if(coeffs_[i]) {
      result ^= (other.coeffs_ << i);
    }
  }
  auto max_size = std::max(degree_, 0ll) + std::max(other.degree_, 0ll);
  coeffs_ = result;
  simplify(max_size);
  return *this;
}

bool G2Poly::null() const {
  return coeffs_.none();
}

G2Poly G2Poly::derivate() const {
  G2Poly p(*this);
  p.coeffs_ >>= 1;
  p.coeffs_ &= derivationMask;
  p.simplify(std::max(degree_, 0ll));
  return p;
}

G2Poly gcd(G2Poly a, G2Poly b) {
  if(a.degree() < b.degree()) {
    G2Poly tmp = a;
    a = b;
    b = tmp;
  }
  //std::cerr << "gcd(\n   - " << a << "\n   - " << b << "\n): " << std::endl;
  G2Poly r = a % b;
  while(!r.null()) {
    //std::cerr << "gcd(\n   - " << b << "\n   - " << r << "\n): " << std::endl;
    a = b;
    b = r;
    r = a % b;
  }
  //std::cerr << "gcd = " << b << std::endl;
  b.simplify(std::max(a.degree(), 0ll));
  return b;
}

G2Poly G2Poly::squareFreePart(uint64_t p) const {
  //std::cerr << "squareFreePart(" << *this << ", " << p << ")" << std::endl;
  G2Poly dp = this->derivate();
  if(dp.null() && this->degree() > 0) {
    G2Poly newp = this->unpower(p);
    return newp.squareFreePart(p);
  }
  G2Poly g = gcd(*this, dp);
  G2Poly r = (*this) / g;
  //std::cerr << "squareFreePart(" << *this << ") = " << (*this) << " / " << g << " = " << r << " (dp = " << dp << ")" << std::endl;
  return r;
}

G2Poly G2Poly::power(uint64_t n) const {
  auto result = (*this);
  result.degree_ = 0;
  for(int i = 0; i <= degree(); i++) {
    if(coeffs_[i]) {
      result.coeffs_[i*n] = coeffs_[i];
      result.degree_ = i*n;
    }
  }
  return result;
}

G2Poly G2Poly::unpower(uint64_t n) const {
  auto result = (*this);
  result.degree_ = 0;
  for(int i = 0; i <= degree(); i++) {
    result.coeffs_[i] = 0;
    if (i%n == 0) {
      if(coeffs_[i]) {
        result.coeffs_[i/n] = true;
        result.degree_ = i/n;
      }
    }
  }
  return result;
}

std::vector<std::pair<G2Poly, uint64_t>> G2Poly::squareFreeFactors(uint64_t p) const {
  uint64_t i = 1;
  auto unit = G2Poly::Xn(0);
  auto f = *this;
  //std::cerr << "squareFreeFactors on " << f << std::endl;
  auto g = this->derivate();
  //std::cerr << "derivate: " << g << std::endl;
  std::vector<std::pair<G2Poly, uint64_t>> result;
  if (!g.null()) {
    auto c = gcd(f, g);
    //std::cerr << "gcd(1): " << c << std::endl;
    auto w = f/c;
    while (w != unit) {
      //std::cerr << "w: " << w << std::endl;
      auto y = gcd(w, c);
      //std::cerr << "gcd(2): " << y << std::endl;
      auto z = w / y;
      if (z != unit) {
        //std::cerr << "Adding(1) " << z << " with degree " << i << std::endl;
        result.push_back(std::make_pair(z, i));
      } else {
        //std::cerr << "Reserving: " << w << " (" << z << ") -> " << i << std::endl;
      }
      i++;
      w = y;
      c = c / y;
    }
    if (c != unit) {
      //std::cerr << "Unpowering " << c << " with degree " << p << std::endl;
      c = c.unpower(p);
      //std::cerr << "Adding(2) " << c << " with degree " << p << std::endl;
      result.push_back(std::make_pair(c, p));
      return result;
    }
    else {
      return result;
    }
  } else {
    f = f.unpower(p);
    //std::cerr << "Recursively calling squareFreeFactors on " << f << std::endl;
    result = f.squareFreeFactors(p);
    for(auto it = result.begin(); it != result.end(); it++) {
      it->second *= p;
    }
    return result;
  }
}

/* Return distinct degree factors for a polynom on a finite field of characteristic p */
std::vector<G2Poly> G2Poly::distinctDegreeFactors(uint64_t p) const {
  std::vector<G2Poly> g;
  auto f = this->squareFreePart(p);
  g.resize(f.degree());
  //auto f = *this;
  //std::cerr << "distinctDegreeFactors(" << *this << ") => " << f << ": " << std::endl;
  for(int64_t i = 0; i < f.degree(); i++) {
    g[i] = gcd(f, G2Poly::Xn(pow(p, i+1)) - G2Poly::Xn(1));
    f = f / g[i];
    //std::cerr << "  factor " << i << ": " << g[i] << std::endl;
  }
  //std::cerr << "squareFreePart(" << *this << ") = " << (*this) << " / " << g << " = " << r << " (dp = " << dp << ")" << std::endl;
  return g;
}

G2Poly G2Poly::equalDegreeFactorize(uint64_t p, uint64_t d) const {
  while(true) {
    auto a = G2Poly::Rand(degree()-1);
    //std::cerr << "EDF(1): " << a << std::endl;
    auto g = gcd(*this, a);
    //std::cerr << "EDF(2): " << g << std::endl;
    auto unit = G2Poly::Xn(0);
    if(g != unit) {
      return g;
    }
    uint64_t m = (std::pow(p, d)-1) / 2;
    auto a2 = a.power(m) + unit;
    //std::cerr << "EDF(3): " << a2 << std::endl;
    g = gcd(*this, a2);
    //std::cerr << "EDF(4): " << g << std::endl;
    if(!a2.null() && g != unit) {
      return g;
    }
  }
}

G2Poly G2Poly::cantorZassenhaus(uint64_t p) const {
  auto sfp = squareFreePart(p);
  if(sfp != (*this)) {
    return sfp;
  }
  auto ddf = distinctDegreeFactors(p);
  auto unit = G2Poly::Xn(0);
  if(ddf[degree()] == (*this)) {
    return *this;
  }
  for(int64_t i = 0; i < degree(); i++) {
    if(ddf[i] != unit) {
      return ddf[i].equalDegreeFactorize(p, i+1);
    }
  }
  std::cerr << "Oups, Cantor-Zassenhaus terminated without finding anything: shouldn't reach this line !!!" << std::endl;
  return *this;
}

G2Poly G2Poly::knuth(uint64_t p) const {
  auto sfp = squareFreePart(p);
  if(sfp != (*this)) {
    return sfp;
  }
  auto ddf = distinctDegreeFactors(p);
  if(ddf[degree()] == (*this)) {
    return *this;
  }
  auto unit = G2Poly::Xn(0);
  for(uint64_t d = 0; d < ddf.size(); d++) {
    if(ddf[d] != unit) {
      std::cerr << "Knuth for " << ddf[d] << std::endl;
      //while(true)
      //  auto t = G2Poly::Rand(degree()-1);
      for(int i = 0; i < ddf[d].degree(); i++) {
        auto t = G2Poly::Xn(2*i+1);
        std::cerr << "  With " << t << std::endl;
        auto u = t;
        for(int j = 2; j < ddf[d].degree(); j++) {
          u = gcd((t+u*u), ddf[d]);
        }
        //std::cerr << "  U = " << u << std::endl;
        if(!u.null() && u != unit) {
          return u;
        }
      }
    }
  }

  std::cerr << "Oups, Knuth didn't work !!!" << std::endl;
  return *this;
}

G2Poly G2Poly::mceliece(uint64_t p) const {
  /* Assume square free polynomial
  auto sfp = squareFreePart(p);
  if(sfp != (*this)) {
    return sfp;
  }
  auto ddf = distinctDegreeFactors(p);
  if(ddf[degree()] == (*this)) {
    return *this;
  }
  */
  //std::cerr << "McEliece for " << (*this) << std::endl;
  auto x = G2Poly::Xn(1);
  uint64_t N = 2;
  auto r = G2Poly::Xn(2);
  while(true) {
    //auto r = G2Poly::Xn(std::pow(2,N)) % (*this);
    r = (r*r) % (*this);
    //std::cerr << "X^(2^" << N << ") mod f = " << r << std::endl;
    if(r == x) {
      break;
    }
    if(r.null()) {
      std::cerr << "ERROR: Bitset size too small ..." << std::endl;
      return (*this);
    }
    N++;
  }
  //std::cerr << "N = " << N << std::endl;
  auto unit = G2Poly::Xn(0);
  for(int64_t i = 1; i < degree(); i++) {
    auto dT = G2Poly::Xn(i) % (*this);
    G2Poly Ti;
    //auto Ti = G2Poly::Xn(i) % (*this);
    for(uint64_t j = 0; j < N; j++) {
      //auto dT = G2Poly::Xn(i << j) % (*this);
      Ti += dT;
      dT = (dT*dT) % (*this);
    }
    //std::cerr << "T[" << i << "] = " << Ti << std::endl;
    if(!Ti.null()) {
      auto Ci = gcd((*this), Ti);
      auto Qi = (*this) / Ci;
      if( (Ci != unit) && (Ci * Qi == (*this)) ) {
        //std::cerr << "C[" << i << "] = " << Ci << " * Q[" << i << "] = " << Qi << " works !!!" << std::endl;
        return Ci;
      }
    }
  }
  std::cerr << "Oups, McEliece didn't work !!!" << std::endl;
  return *this;
}
