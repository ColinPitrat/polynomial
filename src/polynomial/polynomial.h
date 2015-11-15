#pragma once

#include <iostream>
#include <vector>
#include <algorithm>
#include <utility>
#include <cmath>

template <typename T>
class Poly {
/*
  friend std::ostream& operator<< <>(std::ostream &s, const Poly<T> &p);
  friend Poly<T> operator*(const Poly<T> &a, const Poly<T> &b);
  friend Poly<T> operator*<>(const Poly<T> &a, const Poly<T> &b);
  friend Poly<T> operator%<>(const Poly<T> &a, const Poly<T> &b);
  */

  public:
    Poly<T>() = default;
    Poly<T>(std::vector<T> coeffs);

    int degree() const;
    bool null() const;

    void simplify();
    Poly<T> power(int n) const;
    Poly<T> unpower(int n) const;
    Poly<T> derivate() const;
    Poly<T> squareFreePart(int p) const;
    std::vector<std::pair<Poly<T>, int>> squareFreeFactors(int p) const;
    std::vector<Poly<T>> distinctDegreeFactors(int p) const;
    Poly<T> cantorZassenhaus(int p) const;
    Poly<T> equalDegreeFactorize(int p, int i) const;

    Poly<T>& operator-=(const Poly<T> &p);
    Poly<T>& operator+=(const Poly<T> &p);
    Poly<T>& operator*=(const Poly<T> &p);

    static Poly<T> Xn(int n);
    static Poly<T> Rand(int n);

  //protected:
		std::vector<T> coeffs_;
};

template <typename T>
bool operator==(const Poly<T> &a, const Poly<T> &b) {
  if(a.degree() != b.degree()) {
    return false;
  }
	for(size_t i = 0; i < a.coeffs_.size(); i++) {
    if(a.coeffs_[i] != b.coeffs_[i]) {
      return false;
    }
  }
  return true;
}

template <typename T>
bool operator!=(const Poly<T> &a, const Poly<T> &b) {
  return !(a == b);
}

template <typename T>
std::ostream& operator<<(std::ostream &s, const Poly<T> &p) {
	int i = p.degree();
  bool first = true;
  if(p.coeffs_.empty()) s << "0";
	for(auto rit = p.coeffs_.rbegin(); rit != p.coeffs_.rend(); ++rit, --i) {
    if(*rit != T(0)) {
      if(*rit < T(0)) {
        if(!first) {
          s << " - ";
        } else {
          s << "-";
        }
      } else if (!first) {
        s << " + ";
      }
      first = false;
      if((*rit != T(1) && (*rit > T(0) || *rit != T(-1))) || i == 0) {
        if((*rit) < T(0)) {
          s << -(*rit);
        } else {
          s << (*rit);
        }
        if(i > 0) {
          s << "*";
        }
      }
      if(i > 0) {
        s << "X";
        if(i > 1) {
          s << "^" << i;
        }
      }
    }
  }
	return s;
}

template <typename T>
Poly<T> operator-(const Poly<T> &a, const Poly<T> &b) {
  auto p = a;
  p -= b;
  return p;
}

template <typename T>
Poly<T> operator*(const Poly<T> &a, const Poly<T> &b) {
  auto r = a;
  r *= b;
  return r;
}

template <typename T>
Poly<T> operator+(const Poly<T> &a, const Poly<T> &b) {
  auto result = a;
  result += b;
  return result;
}

template <typename T>
Poly<T> operator*(T m, const Poly<T> &p) {
  Poly<T> r(p);
  for(auto it = r.coeffs_.begin(); it != r.coeffs_.end(); it++) {
    *it *= m;
  }
  r.simplify();
  return r;
}

template <typename T>
Poly<T> operator/(const Poly<T> &p, T d) {
  Poly<T> r(p);
  for(auto it = r.coeffs_.begin(); it != r.coeffs_.end(); it++) {
    *it /= d;
  }
  r.simplify();
  return r;
}

template <typename T>
std::pair<Poly<T>, Poly<T>> euclidDivide(const Poly<T> &a, const Poly<T> &b) {
  //std::cout << "euclidDivide(\n   - " << a << "\n   - " << b << "\n) = ";
  Poly<T> r(a);
  Poly<T> q; 
  auto d = b.degree();
  int qsize = std::max(a.degree() - b.degree(), 0)+1;
  q.coeffs_.resize(qsize);
  while(r.degree() >= d) {
    int n = r.degree() - b.degree();
    //std::cout << "N=" << n << ", deg(r)=" << r.degree() << ", deg(b)=" << b.degree() << std::endl;
    // TODO / BUG: divide by 0 when b is not simplified
    q.coeffs_[n] = r.coeffs_[r.degree()] / b.coeffs_[b.degree()];
    r -= q.coeffs_[n] * b * Poly<T>::Xn(n);
    //std::cout << " -> " << q.coeffs_[n] << "x^" << n << "\n    q = " << q << "\n    r = " << r << std::endl;
  }
  r.simplify();
  q.simplify();
  //std::cout << "(" << q << ", " << r << ")" << std::endl;
  return std::make_pair(q, r);
}

template <typename T>
Poly<T> operator%(const Poly<T> &a, const Poly<T> &b) {
  return euclidDivide(a,b).second;
}

template <typename T>
Poly<T> operator/(const Poly<T> &a, const Poly<T> &b) {
  return euclidDivide(a,b).first;
}

template <typename T>
Poly<T>::Poly(std::vector<T> coeffs) {
  coeffs_ = coeffs;
}

template <typename T>
int Poly<T>::degree() const { 
  return coeffs_.size()-1; 
}

template <typename T>
Poly<T> Poly<T>::Xn(int n) {
  Poly<T> p;
  for(int i = 0; i < n; i++) {
    p.coeffs_.push_back(static_cast<T>(0));
  }
  p.coeffs_.push_back(static_cast<T>(1));
  return p;
}

template <typename T>
Poly<T> Poly<T>::Rand(int n) {
  Poly<T> p;
  for(int i = 0; i < n; i++) {
    p.coeffs_.push_back(static_cast<T>(rand()));
  }
  p.simplify();
  return p;
}

template <typename T>
void Poly<T>::simplify() {
  int to_remove = 0;
	for(auto rit = this->coeffs_.rbegin(); rit != this->coeffs_.rend(); ++rit) {
    if(*rit != T(0)) {
      break;
    }
    to_remove++;
  }
  for(int i = 0; i < to_remove; i++) {
    this->coeffs_.pop_back();
  }
}

template <typename T>
Poly<T>& Poly<T>::operator-=(const Poly<T> &other) {
  size_t i = 0;
  for(auto it = other.coeffs_.begin(); it != other.coeffs_.end(); it++) {
    if(i < this->coeffs_.size()) {
      this->coeffs_[i] -= *it;
    } else {
      this->coeffs_.push_back(-(*it));
    }
    i++;
  }
  this->simplify();
  return *this;
}

template <typename T>
Poly<T>& Poly<T>::operator*=(const Poly<T> &other) {
  auto deg = degree() + other.degree();
  Poly<T> result;
  result.coeffs_.resize(deg+1);
  int i(0);
  for(auto itA = coeffs_.begin(); itA != coeffs_.end(); itA++) {
    int j(0);
    for(auto itB = other.coeffs_.begin(); itB != other.coeffs_.end(); itB++) {
      result.coeffs_[i+j] += (*itB)*(*itA);
      j++;
    }
    i++;
  }
  *this = result;
  return *this;
}

template <typename T>
Poly<T>& Poly<T>::operator+=(const Poly<T> &other) {
  size_t i = 0;
  for(auto it = other.coeffs_.begin(); it != other.coeffs_.end(); it++) {
    if(i < this->coeffs_.size()) {
      this->coeffs_[i] += *it;
    } else {
      this->coeffs_.push_back(*it);
    }
    i++;
  }
  this->simplify();
  return *this;
}

template <typename T>
bool Poly<T>::null() const {
  for(auto it = this->coeffs_.begin(); it != this->coeffs_.end(); ++it) {
    if(*it != T(0)) {
      return false;
    }
  }
  return true;
}

template <typename T>
Poly<T> Poly<T>::derivate() const {
  Poly<T> p;
  int i = 0;
  for(auto it = this->coeffs_.begin(); it != this->coeffs_.end(); ++it) {
    if(i > 0) {
      p.coeffs_.push_back(T(i)*(*it));
    }
    i++;
  }
  p.simplify();
  return p;
}

template <typename T>
Poly<T> gcd(Poly<T> a, Poly<T> b) {
  if(a.degree() < b.degree()) {
    Poly<T> tmp = a;
    a = b;
    b = tmp;
  }
  //std::cout << "gcd(\n   - " << a << "\n   - " << b << "\n): " << std::endl;
  Poly<T> r = a % b;
  while(!r.null()) {
    //std::cout << "gcd(\n   - " << b << "\n   - " << r << "\n): " << std::endl;
    a = b;
    b = r;
    r = a % b; 
  }
  //std::cout << "gcd = " << b;
  // Normalize the polynom:
  b = b / b.coeffs_[b.degree()];
  //std::cout << " = " << b << std::endl;
  return b;
}

template <typename T>
Poly<T> Poly<T>::squareFreePart(int p) const {
  //std::cout << "squareFreePart(" << *this << ", " << p << ")" << std::endl;
  Poly<T> dp = this->derivate();
  if(dp.null() && this->degree() > 0) { 
    Poly<T> newp = this->unpower(p);
    return newp.squareFreePart(p);
  }
  Poly<T> g = gcd(*this, dp);
  Poly<T> r = (*this) / g;
  //std::cout << "squareFreePart(" << *this << ") = " << (*this) << " / " << g << " = " << r << " (dp = " << dp << ")" << std::endl;
  return r;
}

template <typename T>
Poly<T> Poly<T>::power(int n) const {
  auto result = (*this);
  for(int i = 1; i < n; i++) {
    result = result * (*this);
  }
  return result;
}

template <typename T>
Poly<T> Poly<T>::unpower(int n) const {
  auto result = (*this);
  for(int i = 0; i <= degree(); i++) {
    result.coeffs_[i] = 0;
    if (i%n == 0) {
      result.coeffs_[i/n] = coeffs_[i];
    }
  }
  result.simplify();
  return result;
}

template <typename T>
std::vector<std::pair<Poly<T>, int>> Poly<T>::squareFreeFactors(int p) const {
  int i = 1;
  auto unit = Poly<T>::Xn(0);
  auto f = *this;
  //std::cout << "squareFreeFactors on " << f << std::endl;
  auto g = this->derivate();
  //std::cout << "derivate: " << g << std::endl;
  std::vector<std::pair<Poly<T>, int>> result;
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
template <typename T>
std::vector<Poly<T>> Poly<T>::distinctDegreeFactors(int p) const {
  std::vector<Poly<T>> g;
  auto f = this->squareFreePart(p);
  //auto f = *this;
  g.resize(f.degree());
  //std::cout << "distinctDegreeFactors(" << *this << ") => " << f << ": " << std::endl;
  for(int i = 0; i < f.degree(); i++) {
    g[i] = gcd(f, Poly<T>::Xn(pow(p, i+1)) - Poly<T>::Xn(1));
    g[i].simplify();
    f = f / g[i];
    //std::cout << "  factor " << i << ": " << g[i] << std::endl;
  }
  //std::cout << "squareFreePart(" << *this << ") = " << (*this) << " / " << g << " = " << r << " (dp = " << dp << ")" << std::endl;
  return g;
}

template <typename T>
Poly<T> Poly<T>::equalDegreeFactorize(int p, int d) const {
  while(true) {
    auto a = Poly<T>::Rand(degree()-1);
    //std::cout << "EDF(1): " << a << std::endl;
    auto g = gcd(*this, a);
    //std::cout << "EDF(2): " << g << std::endl;
    auto unit = Poly<T>::Xn(0);
    if(g != unit) {
      return g;
    }
    int m = (std::pow(p, d)-1) / 2;
    auto a2 = a.power(m) + unit;
    //std::cout << "EDF(3): " << a2 << std::endl;
    g = gcd(*this, a2);
    //std::cout << "EDF(4): " << g << std::endl;
    if(!a2.null() && g != unit) {
      return g;
    }
  }
}

template <typename T>
Poly<T> Poly<T>::cantorZassenhaus(int p) const {
  auto sfp = squareFreePart(p);
  if(sfp != (*this)) {
    return sfp;
  }
  auto ddf = distinctDegreeFactors(p);
  auto unit = Poly<T>::Xn(0);
  if(ddf[degree()] == (*this)) {
    return *this;
  }
  for(int i = 0; i < degree(); i++) {
    if(ddf[i] != unit) {
      return ddf[i].equalDegreeFactorize(p, i+1);
    }
  }
  std::cerr << "Oups, Cantor-Zassenhaus terminated without finding anything: shouldn't reach this line !!!" << std::endl;
  return *this;
}
