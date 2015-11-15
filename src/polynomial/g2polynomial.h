#pragma once

/*
 * Polynomials on GF(2) implemented using a vector of non-null indices.
 * This is fast on large sparse polynomials.
 */
#include <iostream>
#include <vector>
#include <algorithm>
#include <utility>
#include <cmath>
#include <bitset>

class G2Poly {
  public:
    G2Poly() = default;
    G2Poly(std::vector<uint64_t> coeffs);

    int64_t degree() const;
    bool null() const;

    G2Poly timesXn(uint64_t n) const;
    void minusTimesXn(const G2Poly &b, uint64_t n);
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
		std::vector<uint64_t> coeffs_;
};

bool operator==(const G2Poly &a, const G2Poly &b) {
  return a.coeffs_ == b.coeffs_;
}

bool operator!=(const G2Poly &a, const G2Poly &b) {
  return a.coeffs_ != b.coeffs_;
}

std::ostream& operator<<(std::ostream &s, const G2Poly &p) {
  bool first = true;
  if(p.coeffs_.empty()) s << "0";
  for(auto it = p.coeffs_.begin(); it != p.coeffs_.end(); it++) {
    if(!first) {
      s << " + ";
    }
    first = false;
    if(*it == 0) {
      s << "1";
    } else {
      s << "X";
      if(*it > 1) {
        s << "^" << (*it);
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

G2Poly G2Poly::timesXn(uint64_t n) const {
  auto r = (*this);
  auto size = r.coeffs_.size();
  for(size_t i = 0; i < size; i++) {
    r.coeffs_[i] += n;
  }
  /*
  auto itEnd = r.coeffs_.end();
  for(auto it = r.coeffs_.begin(); it != itEnd; it++) {
    (*it) += n;
  }
  */
  return r;
}

void G2Poly::minusTimesXn(const G2Poly &b, uint64_t n) {
  std::vector <uint64_t> result;
  result.reserve(std::max(b.coeffs_.size(), coeffs_.size()));
  size_t i(0), j(0);
  auto size = coeffs_.size();
  auto size2 = b.coeffs_.size();
  auto ci = coeffs_[i];
  auto cj = b.coeffs_[j] + n;
  while(i != size && j != size2) {
    if(ci > cj) {
      result.push_back(ci);
      ci = coeffs_[++i];
    } else if(ci < cj) {
      result.push_back(cj);
      cj = b.coeffs_[++j] + n;
    } else {
      ci = coeffs_[++i];
      cj = b.coeffs_[++j] + n;
    }
  }
  // Insert remaining coeffs from other as they are obviously not in this
  while(i != size) {
    result.push_back(coeffs_[i]);
    i++;
  }
  while(j != size2) {
    result.push_back(b.coeffs_[j] + n);
    j++;
  }
  this->coeffs_.swap(result);
}

std::pair<G2Poly, G2Poly> euclidDivide(const G2Poly &a, const G2Poly &b) {
  //std::cout << "euclidDivide(\n   - " << a << "\n   - " << b << "\n) = ";
  G2Poly r(a);
  G2Poly q; 
  auto db = b.degree();
  auto dr = r.degree();
  while(dr >= db) {
    uint64_t n = dr - db;
    q.coeffs_.push_back(n);
    //r -= b.timesXn(n);
    r.minusTimesXn(b, n);
    dr = r.degree();
    //std::cout << "    n = " << n << " (" << q << ", " << r << ")" << std::endl;
  }
  //std::cout << "(" << q << ", " << r << ")" << std::endl;
  return std::make_pair(q, r);
}

G2Poly operator%(const G2Poly &a, const G2Poly &b) {
  G2Poly r(a);
  while(r.degree() >= b.degree()) {
    int64_t n = r.degree() - b.degree();
    //r -= b.timesXn(n);
    r.minusTimesXn(b, n);
  }
  return r;
}

G2Poly operator/(const G2Poly &a, const G2Poly &b) {
  return euclidDivide(a,b).first;
}

G2Poly::G2Poly(std::vector<uint64_t> coeffs) {
  uint64_t i = coeffs.size();
  for(auto rit = coeffs.rbegin(); rit != coeffs.rend(); rit++) {
    i--;
    if(*rit) {
      coeffs_.push_back(i);
    }
  }
}

int64_t G2Poly::degree() const { 
  if(coeffs_.empty()) {
    return -1;
  }
  return coeffs_[0]; 
}

G2Poly G2Poly::Xn(uint64_t n) {
  G2Poly p;
  p.coeffs_.push_back(n);
  return p;
}

G2Poly G2Poly::Rand(uint64_t n) {
  G2Poly p;
  for(uint64_t i = n; i > 0; i--) {
    if(rand() & 1) {
      p.coeffs_.push_back(i-1);
    }
  }
  return p;
}

G2Poly& G2Poly::operator-=(const G2Poly &other) {
  return operator+=(other);
}

/* Old code: (working but slow because of insert / erase)
G2Poly& G2Poly::operator+=(const G2Poly &other) {
  auto it = coeffs_.begin(); 
  auto it2 = other.coeffs_.begin(); 
  while(it != coeffs_.end() && it2 != other.coeffs_.end()) {
    // Coeffs in both polynomials => will be zero in sum
    if(*it == *it2) {
      it = coeffs_.erase(it);
      it2++;
    } else if(*it > *it2) {
      it++;
    } else if(*it < *it2) {
      it = coeffs_.insert(it, *it2);
      it2++;
    }
  }
  // Insert remaining coeffs from other as they are obviously not in this
  while(it != other.coeffs_.end()) {
    coeffs_.push_back(*it);
    it++;
  }
  while(it2 != other.coeffs_.end()) {
    coeffs_.push_back(*it2);
    it2++;
  }
  return *this;
}
*/

G2Poly& G2Poly::operator+=(const G2Poly &other) {
  std::vector <uint64_t> result;
  result.reserve(std::max(other.coeffs_.size(), coeffs_.size()));
  size_t i(0), j(0);
  auto size = coeffs_.size();
  auto size2 = other.coeffs_.size();
  while(i != size && j != size2) {
    if(coeffs_[i] > other.coeffs_[j]) {
      result.push_back(coeffs_[i]);
      i++;
    } else if(coeffs_[i] < other.coeffs_[j]) {
      result.push_back(other.coeffs_[j]);
      j++;
    } else {
      i++;
      j++;
    }
  }
  // Insert remaining coeffs from other as they are obviously not in this
  while(i != size) {
    result.push_back(coeffs_[i]);
    i++;
  }
  while(j != size2) {
    result.push_back(other.coeffs_[j]);
    j++;
  }
  this->coeffs_.swap(result);
  return *this;
}

G2Poly& G2Poly::operator*=(const G2Poly &other) {
  std::vector<uint64_t> result;
  //std::cout << "Product of " << *this << " - " << other << std::endl;
  for(auto it = coeffs_.begin(); it != coeffs_.end(); ++it) {
    auto it3 = result.begin();
    for(auto it2 = other.coeffs_.begin(); it2 != other.coeffs_.end(); ++it2) {
      uint64_t toToggle = (*it) + (*it2);
      //std::cout << "  To toggle: " << toToggle << std::endl;
      while(it3 != result.end() && (*it3) > toToggle) {
        it3++;
      }
      if(it3 != result.end() && (*it3) == toToggle) {
        //std::cout << "    removing: " << toToggle << std::endl;
        it3 = result.erase(it3);
      } else {
        //std::cout << "    setting: " << toToggle << std::endl;
        it3 = result.insert(it3, toToggle);
      }
    }
  }
  this->coeffs_.swap(result);
  return *this;
}

bool G2Poly::null() const {
  return coeffs_.empty();
}

G2Poly G2Poly::derivate() const {
  G2Poly p;
  for(auto it = coeffs_.begin(); it != coeffs_.end(); it++) {
    if((*it) & 1) {
      p.coeffs_.push_back((*it)-1);
    }
  }
  return p;
}

G2Poly gcd(G2Poly a, G2Poly b) {
  if(a.degree() < b.degree()) {
    G2Poly tmp = a;
    a = b;
    b = tmp;
  }
  //std::cout << "gcd(\n   - " << a << "\n   - " << b << "\n): " << std::endl;
  G2Poly r = a % b;
  while(!r.null()) {
    //std::cout << "gcd(\n   - " << b << "\n   - " << r << "\n): " << std::endl;
    a = b;
    b = r;
    r = a % b; 
  }
  //std::cout << "gcd = " << b << std::endl;
  return b;
}

G2Poly G2Poly::squareFreePart(uint64_t p) const {
  //std::cout << "squareFreePart(" << *this << ", " << p << ")" << std::endl;
  G2Poly dp = this->derivate();
  if(dp.null() && this->degree() > 0) { 
    G2Poly newp = this->unpower(p);
    return newp.squareFreePart(p);
  }
  G2Poly g = gcd(*this, dp);
  G2Poly r = (*this) / g;
  //std::cout << "squareFreePart(" << *this << ") = " << (*this) << " / " << g << " = " << r << " (dp = " << dp << ")" << std::endl;
  return r;
}

G2Poly G2Poly::power(uint64_t n) const {
  G2Poly result;
  for(auto it = coeffs_.begin(); it != coeffs_.end(); ++it) {
    result.coeffs_.push_back((*it)*n);
  }
  return result;
}

G2Poly G2Poly::unpower(uint64_t n) const {
  G2Poly result;
  for(auto it = coeffs_.begin(); it != coeffs_.end(); ++it) {
    if ((*it)%n == 0) {
      result.coeffs_.push_back((*it)/n);
    }
  }
  return result;
}

std::vector<std::pair<G2Poly, uint64_t>> G2Poly::squareFreeFactors(uint64_t p) const {
  uint64_t i = 1;
  auto unit = G2Poly::Xn(0);
  auto f = *this;
  //std::cout << "squareFreeFactors on " << f << std::endl;
  auto g = this->derivate();
  //std::cout << "derivate: " << g << std::endl;
  std::vector<std::pair<G2Poly, uint64_t>> result;
  if (!g.null()) {
    auto c = gcd(f, g);
    //std::cout << "gcd(1): " << c << std::endl;
    auto w = f/c;
    while (w != unit) {
      //std::cout << "w: " << w << std::endl;
      auto y = gcd(w, c);
      //std::cout << "gcd(2): " << y << std::endl;
      auto z = w / y;
      if (z != unit) {
        //std::cout << "Adding(1) " << z << " with degree " << i << std::endl;
        result.push_back(std::make_pair(z, i));
      } else {
        //std::cout << "Reserving: " << w << " (" << z << ") -> " << i << std::endl;
      }
      i++;
      w = y;
      c = c / y;
    }
    if (c != unit) {
      //std::cout << "Unpowering " << c << " with degree " << p << std::endl;
      c = c.unpower(p);
      //std::cout << "Adding(2) " << c << " with degree " << p << std::endl;
      result.push_back(std::make_pair(c, p));
      return result;
    }
    else {
      return result;
    }
  } else {
    f = f.unpower(p);
    //std::cout << "Recursively calling squareFreeFactors on " << f << std::endl;
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
  //std::cout << "distinctDegreeFactors(" << *this << ") => " << f << ": " << std::endl;
  for(int64_t i = 0; i < f.degree(); i++) {
    g[i] = gcd(f, G2Poly::Xn(pow(p, i+1)) - G2Poly::Xn(1));
    f = f / g[i];
    //std::cout << "  factor " << i << ": " << g[i] << std::endl;
  }
  //std::cout << "squareFreePart(" << *this << ") = " << (*this) << " / " << g << " = " << r << " (dp = " << dp << ")" << std::endl;
  return g;
}

G2Poly G2Poly::equalDegreeFactorize(uint64_t p, uint64_t d) const {
  while(true) {
    auto a = G2Poly::Rand(degree()-1);
    //std::cout << "EDF(1): " << a << std::endl;
    auto g = gcd(*this, a);
    //std::cout << "EDF(2): " << g << std::endl;
    auto unit = G2Poly::Xn(0);
    if(g != unit) {
      return g;
    }
    uint64_t m = (std::pow(p, d)-1) / 2;
    auto a2 = a.power(m) + unit;
    //std::cout << "EDF(3): " << a2 << std::endl;
    g = gcd(*this, a2);
    //std::cout << "EDF(4): " << g << std::endl;
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
      std::cout << "Knuth for " << ddf[d] << std::endl;
      //while(true)
      //  auto t = G2Poly::Rand(degree()-1);
      for(int i = 0; i < ddf[d].degree(); i++) {
        auto t = G2Poly::Xn(2*i+1);
        std::cout << "  With " << t << std::endl;
        auto u = t;
        for(int j = 2; j < ddf[d].degree(); j++) {
          u = gcd((t+u*u), ddf[d]);
        }
        //std::cout << "  U = " << u << std::endl;
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
  std::cout << "McEliece for " << (*this) << std::endl;
  auto x = G2Poly::Xn(1);
  uint64_t N = 2;
  auto r = G2Poly::Xn(2);
  while(true) {
    //auto r = G2Poly::Xn(std::pow(2,N)) % (*this);
    r = (r*r) % (*this);
    //std::cout << "X^(2^" << N << ") mod f = " << r << std::endl;
    if(r == x) {
      break;
    }
    if(r.null()) {
      std::cout << "ERROR: Bitset size too small ..." << std::endl;
      return (*this);
    }
    N++;
  }
  std::cout << "N = " << N << std::endl;
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
    std::cout << "T[" << i << "] = " << Ti << std::endl;
    if(!Ti.null()) {
      auto Ci = gcd((*this), Ti);
      auto Qi = (*this) / Ci;
      if( (Ci != unit) && (Ci * Qi == (*this)) ) {
        std::cout << "C[" << i << "] = " << Ci << " * Q[" << i << "] = " << Qi << " works !!!" << std::endl;
        return Ci;
      }
    }
  }
  std::cerr << "Oups, McEliece didn't work !!!" << std::endl;
  return *this;
}
