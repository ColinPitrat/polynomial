#pragma once

#include <iostream>

template <unsigned int n>
class FFElem {
  public:
    static const unsigned int characteristic;

    FFElem<n>();
    FFElem<n>(unsigned int v);
    FFElem<n>(int v);

    void operator+=(const FFElem<n> &other);
    void operator-=(const FFElem<n> &other);
    void operator*=(const FFElem<n> &other);
    void operator/=(const FFElem<n> &other);

    int compare(const FFElem<n> &other) const;
    void toStream(std::ostream &s) const;
    FFElem<n> inverse() const;

  protected:
    void normalize();

		unsigned int val_;
};

template <unsigned int n>
FFElem<n>::FFElem() : val_(0) {};

template <unsigned int n>
FFElem<n>::FFElem(unsigned int v) : val_(v) {
  normalize();
};

template <unsigned int n>
FFElem<n>::FFElem(int v) {
  if(v < 0) { 
    v += (-v/n + 1)*n; 
  } 
  val_ = v; 
  normalize();
};

template <unsigned int n>
const unsigned int FFElem<n>::characteristic = n;

template <unsigned int n>
void FFElem<n>::normalize() {
  if(this->val_ >= n) {
    this->val_ = this->val_ % n;
  }
};

template <unsigned int n>
void FFElem<n>::operator+=(const FFElem<n> &other) {
  val_ += other.val_;
  normalize();
}

template <unsigned int n>
void FFElem<n>::operator-=(const FFElem<n> &other) {
  while(other.val_ > val_) { val_ += n; }
  val_ -= other.val_;
  normalize();
}

template <unsigned int n>
void FFElem<n>::operator*=(const FFElem<n> &other) {
  val_ *= other.val_;
  normalize();
}

template <unsigned int n>
void FFElem<n>::operator/=(const FFElem<n> &other) {
  operator*=(other.inverse());
  normalize();
}

template <unsigned int n>
int FFElem<n>::compare(const FFElem<n> &other) const {
  if(val_ < other.val_) {
    return -1;
  } else if(val_ > other.val_) {
    return 1;
  }
  return 0;
}

template <unsigned int n>
bool operator<(const FFElem<n> &a, const FFElem<n> &b) {
  return a.compare(b) < 0;
}

template <unsigned int n>
bool operator<=(const FFElem<n> &a, const FFElem<n> &b) {
  return a.compare(b) <= 0;
}

template <unsigned int n>
bool operator>(const FFElem<n> &a, const FFElem<n> &b) {
  return a.compare(b) > 0;
}

template <unsigned int n>
bool operator>=(const FFElem<n> &a, const FFElem<n> &b) {
  return a.compare(b) >= 0;
}

template <unsigned int n>
bool operator==(const FFElem<n> &a, const FFElem<n> &b) {
  return a.compare(b) == 0;
}

template <unsigned int n>
bool operator!=(const FFElem<n> &a, const FFElem<n> &b) {
  return a.compare(b) != 0;
}

template <unsigned int n>
FFElem<n> operator+(const FFElem<n> &a, const FFElem<n> &b) {
  auto r = a;
  r += b;
  return r;
}

template <unsigned int n>
FFElem<n> operator-(const FFElem<n> &a) {
  FFElem<n> r(0);
  r -= a;
  return r;
}

template <unsigned int n>
FFElem<n> operator-(const FFElem<n> &a, const FFElem<n> &b) {
  auto r = a;
  r -= b;
  return r;
}

template <unsigned int n>
FFElem<n> operator*(const FFElem<n> &a, const FFElem<n> &b) {
  auto r = a;
  r *= b;
  return r;
}

template <unsigned int n>
FFElem<n> operator/(const FFElem<n> &a, const FFElem<n> &b) {
  auto r = a;
  r /= b;
  return r;
}

template <unsigned int n>
FFElem<n> FFElem<n>::inverse() const {
  FFElem<n> self(*this);
  for(unsigned int i = 0; i < n; i++) {
    FFElem<n> inv(i);
    FFElem<n> unit(1);
    if(inv*self == unit) {
      return inv;
    }
  }
  return FFElem<n>(0);
}

template <unsigned int n>
void FFElem<n>::toStream(std::ostream &s) const {
  s << val_ << "[" << n << "]";
}

template <unsigned int n>
std::ostream &operator<<(std::ostream &s, const FFElem<n> &v) {
  v.toStream(s);
  return s;
}
