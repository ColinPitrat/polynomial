#define BOOST_TEST_MODULE TestPolynomial
#include <boost/test/included/unit_test.hpp>
#include "polynomial.h"
#include "../finite_fields/finite_fields.h"
#include <sstream>
#include <boost/rational.hpp>
#include <boost/multiprecision/gmp.hpp>
#include <boost/multiprecision/cpp_int.hpp>
#include <algorithm>
#include <utility>
#include <cstdlib>

struct PolynomialFixture
{
};

BOOST_FIXTURE_TEST_SUITE(TestPolynomial, PolynomialFixture)

BOOST_AUTO_TEST_CASE( test_default_constructed_is_null_and_has_degree_minus_1 ) 
{
  Poly<int> p;

  BOOST_CHECK_EQUAL(true, p.null());
  BOOST_CHECK_EQUAL(-1, p.degree());
}

BOOST_AUTO_TEST_CASE( test_constructed_X2_has_degree_2 ) 
{
  Poly<int> p({0, 0, 1});

  BOOST_CHECK_EQUAL(2, p.degree());
}

BOOST_AUTO_TEST_CASE( test_constant_has_degree_0 ) 
{
  Poly<int> p({42});

  BOOST_CHECK_EQUAL(0, p.degree());
}

BOOST_AUTO_TEST_CASE( test_X0_has_degree_0 ) 
{
  Poly<int> x0 = Poly<int>::Xn(0);

  BOOST_CHECK_EQUAL(0, x0.degree());
}

BOOST_AUTO_TEST_CASE( test_X5_has_degree_5 ) 
{
  Poly<int> x5 = Poly<int>::Xn(5);

  BOOST_CHECK_EQUAL(5, x5.degree());
}

BOOST_AUTO_TEST_CASE( test_constructed_X3_equal_X3 ) 
{
  Poly<int> x3 = Poly<int>::Xn(3);
  Poly<int> cx3({0, 0, 0, 1});

  BOOST_CHECK_EQUAL(x3, cx3);
}

BOOST_AUTO_TEST_CASE( test_X5_not_equal_X3 ) 
{
  Poly<int> x3 = Poly<int>::Xn(3);
  Poly<int> x5 = Poly<int>::Xn(5);

  BOOST_CHECK(x3 != x5);
}

BOOST_AUTO_TEST_CASE( test_constructed_3X3_equal_3X3 ) 
{
  Poly<int> c3x3({0, 0, 0, 3});
  Poly<int> x3 = Poly<int>::Xn(3);

  Poly<int> m3x3 = 3*x3;

  BOOST_CHECK_EQUAL(m3x3, c3x3);
}

BOOST_AUTO_TEST_CASE( test_usual_product ) 
{
  Poly<int> p1({1, 1});
  Poly<int> p2({1, 2});
  Poly<int> result({1, 3, 2});

  BOOST_CHECK_EQUAL(result, p1*p2);
}

BOOST_AUTO_TEST_CASE( test_usual_product2 ) 
{
  Poly<int> p1({5, 2, 3});
  Poly<int> p2({-6, -5, 1});
  Poly<int> result({-30, -37, -23, -13, 3});

  BOOST_CHECK_EQUAL(result, p1*p2);
}

BOOST_AUTO_TEST_CASE( test_toStream_of_null ) 
{
  Poly<int> p;
  std::ostringstream oss;

  oss << p;

  BOOST_CHECK_EQUAL(oss.str(), "0");
}

BOOST_AUTO_TEST_CASE( test_toStream_of_1 ) 
{
  Poly<int> p({1});
  std::ostringstream oss;

  oss << p;

  BOOST_CHECK_EQUAL(oss.str(), "1");
}

BOOST_AUTO_TEST_CASE( test_toStream_of_minus_1 ) 
{
  Poly<int> p({-1});
  std::ostringstream oss;

  oss << p;

  BOOST_CHECK_EQUAL(oss.str(), "-1");
}

BOOST_AUTO_TEST_CASE( test_toStream_of_X ) 
{
  Poly<int> p({0, 1});
  std::ostringstream oss;

  oss << p;

  BOOST_CHECK_EQUAL(oss.str(), "X");
}

BOOST_AUTO_TEST_CASE( test_toStream_of_minus_X ) 
{
  Poly<int> p({0, -1});
  std::ostringstream oss;

  oss << p;

  BOOST_CHECK_EQUAL(oss.str(), "-X");
}

BOOST_AUTO_TEST_CASE( test_toStream_of_minus_2X ) 
{
  Poly<int> p({0, -2});
  std::ostringstream oss;

  oss << p;

  BOOST_CHECK_EQUAL(oss.str(), "-2*X");
}

BOOST_AUTO_TEST_CASE( test_toStream_of_X2_plus_2X_minus_1 ) 
{
  Poly<int> p({-1, 2, 1});
  std::ostringstream oss;

  oss << p;

  BOOST_CHECK_EQUAL(oss.str(), "X^2 + 2*X - 1");
}

BOOST_AUTO_TEST_CASE( test_toStream_of_some_polynomial ) 
{
  Poly<int> p({0, 1, -2, 3, 0, -4});
  std::ostringstream oss;

  oss << p;

  BOOST_CHECK_EQUAL(oss.str(), "-4*X^5 + 3*X^3 - 2*X^2 + X");
}

BOOST_AUTO_TEST_CASE( test_toStream_on_finite_field_with_minus_1 ) 
{
  Poly<FFElem<3>> p({2, 1, 2, 0, 1, 2});
  std::ostringstream oss;

  oss << p;

  BOOST_CHECK_EQUAL(oss.str(), "2[3]*X^5 + X^4 + 2[3]*X^2 + X + 2[3]");
}

BOOST_AUTO_TEST_CASE( test_derivate_of_polynomial ) 
{
  Poly<int> p({6, 7, 1, 3});
  Poly<int> dp({7, 2, 9});

  BOOST_CHECK_EQUAL(p.derivate(), dp);
}

BOOST_AUTO_TEST_CASE( test_derivate_of_null ) 
{
  Poly<int> p;
  Poly<int> dp;

  BOOST_CHECK_EQUAL(p.derivate(), dp);
}

BOOST_AUTO_TEST_CASE( test_remainder ) 
{
  Poly<int> a({0, -2, 3, -1, -1, 1});
  Poly<int> b({1, -1, 1});
  Poly<int> r({-1, 1});

  BOOST_CHECK_EQUAL(r, a%b);
}

BOOST_AUTO_TEST_CASE( test_divide_no_rem ) 
{
  Poly<int> a({0, 0, 6, 3});
  Poly<int> b({2, 1});
  Poly<int> r({0, 0, 3});

  BOOST_CHECK_EQUAL(r, a/b);
}

BOOST_AUTO_TEST_CASE( test_divide_with_rem ) 
{
  Poly<int> a({1, 0, 6, 3});
  Poly<int> b({2, 1});
  Poly<int> r({0, 0, 3});

  BOOST_CHECK_EQUAL(r, a/b);
}

BOOST_AUTO_TEST_CASE( test_gcd ) 
{
  Poly<int> a({0, -2, 3, -1, -1, 1});
  Poly<int> b({1, -1, 1});
  Poly<int> r({-1, 1});

  BOOST_CHECK_EQUAL(r, a%b);
}

BOOST_AUTO_TEST_CASE( test_gcd_rational ) 
{
  // A(x) = x^4 - 2x^2 + 1 = (x - 1)(x - 1)(x + 1)(x + 1)
  // B(x) = 4x^3 - 4x = 4x(x - 1)(x + 1)
  // pgcd(A, B) = (x - 1)(x + 1)
  Poly<boost::rational<int>> a({1, 0, -2, 0, 1});
  Poly<boost::rational<int>> b = a.derivate();
  Poly<boost::rational<int>> r({-1, 0, 1});

  BOOST_CHECK_EQUAL(r, gcd(a,b));
}

BOOST_AUTO_TEST_CASE( test_square_free ) 
{
  Poly<boost::rational<int>> a({1, 0, -2, 0, 1});
  Poly<boost::rational<int>> r({-1, 0, 1});

  BOOST_CHECK_EQUAL(a.squareFreePart(0), r);
}

BOOST_AUTO_TEST_CASE( test_square_free2 ) 
{
  Poly<boost::rational<int64_t>> a({-96, 16, 304, -24, -358, -7, 192, 23, -46, -9, 4, 1});
  Poly<boost::rational<int64_t>> r({12, 4, -15, -5, 3, 1});

  BOOST_CHECK_EQUAL(a.squareFreePart(0), r);
}

BOOST_AUTO_TEST_CASE( test_square_free_null_derivative ) 
{
  FFElem<2> zero(0);
  FFElem<2> one(1);
  Poly<FFElem<2>> a({one, zero, one, zero, one});
  Poly<FFElem<2>> r({one, one, one});

  BOOST_CHECK_EQUAL(a.squareFreePart(2), r);
}

// Exemples with large precision integers
typedef boost::multiprecision::int512_t int512_t;

BOOST_AUTO_TEST_CASE( test_int512_t ) 
{
  Poly<boost::rational<int512_t>> a({int512_t(1), int512_t(0), int512_t(2)});

  BOOST_CHECK_EQUAL(a.degree(), 2);
}

typedef boost::multiprecision::mpz_int mpint;
BOOST_AUTO_TEST_CASE( test_mpz_int ) 
{
  mpint zero = 0;
  mpint one = 1;
  mpint two = 2;
  Poly<boost::rational<mpint>> a({one, zero, two});

  BOOST_CHECK_EQUAL(a.degree(), 2);
}

/* Doesn't terminate 
BOOST_AUTO_TEST_CASE( test_distinct_degree_factors_finite_field_GF2 ) 
{
  FFElem<2> zero(0);
  FFElem<2> one(1);
  Poly<FFElem<2>> f({one, zero, zero, one, zero, zero, one});

  auto result = f.distinctDegreeFactors(2);

  BOOST_CHECK_EQUAL(result.size(), 6);
}
*/

BOOST_AUTO_TEST_CASE( test_square_free_factors ) 
{
  // Example from wikipedia: 
  // https://en.wikipedia.org/wiki/Factorization_of_polynomials_over_finite_fields#Example_of_a_square-free_factorization
  FFElem<3> zero(0);
  FFElem<3> one(1);
  FFElem<3> two(2);
  Poly<FFElem<3>> a({one, zero, two, two, zero, one, one, zero, two, two, zero, one});

  Poly<FFElem<3>> f1({one, one});
  Poly<FFElem<3>> f2({one, zero, one});
  Poly<FFElem<3>> f3({two, one});

  auto result = a.squareFreeFactors(3);

  // Sort by exponents because we don't care of the order of the result in the vector
  std::sort(result.begin(), result.end(), [](const std::pair<Poly<FFElem<3>>, int>& lhs, const std::pair<Poly<FFElem<3>>, int>& rhs) { return lhs.second < rhs.second; } );
  BOOST_CHECK_EQUAL(result.size(), 3);
  BOOST_CHECK_EQUAL(result[0].first, f1);
  BOOST_CHECK_EQUAL(result[0].second, 1);
  BOOST_CHECK_EQUAL(result[1].first, f2);
  BOOST_CHECK_EQUAL(result[1].second, 3);
  BOOST_CHECK_EQUAL(result[2].first, f3);
  BOOST_CHECK_EQUAL(result[2].second, 4);
}

BOOST_AUTO_TEST_CASE( test_square_free_factors_on_GF2 ) 
{
  FFElem<2> zero(0);
  FFElem<2> one(1);
  // 11x11 = 101
  Poly<FFElem<2>> a({one, zero, one});

  Poly<FFElem<2>> f1({one, one});

  auto result = a.squareFreeFactors(2);

  BOOST_CHECK_EQUAL(result.size(), 1);
  BOOST_CHECK_EQUAL(result[0].first, f1);
  BOOST_CHECK_EQUAL(result[0].second, 2);
}

BOOST_AUTO_TEST_CASE( test_distinct_degree_factors_on_nintendo_1 )
{
  FFElem<2> zero(0);
  FFElem<2> one(1);
  // 73AF = 111001110101111 
  Poly<FFElem<2>> a({one, one, one, one, zero, one, zero, one, one, one, zero, zero, one, one, one});
  Poly<FFElem<2>> unit({one});
  // 83 = 10000011
  Poly<FFElem<2>> f1({one, one, zero, zero, zero, zero, zero, one});
  // E5 = 11100101
  Poly<FFElem<2>> f2({one, zero, one, zero, zero, one, one, one});

  auto result = a.distinctDegreeFactors(2);

  BOOST_CHECK_EQUAL(result.size(), 14);
  BOOST_CHECK_EQUAL(result[0], unit);
  BOOST_CHECK_EQUAL(result[1], unit);
  BOOST_CHECK_EQUAL(result[2], unit);
  BOOST_CHECK_EQUAL(result[3], unit);
  BOOST_CHECK_EQUAL(result[4], unit);
  BOOST_CHECK_EQUAL(result[5], unit);
  BOOST_CHECK_EQUAL(result[6], a);
  BOOST_CHECK(result[7].null());
  BOOST_CHECK(result[8].null());
  BOOST_CHECK(result[9].null());
  BOOST_CHECK(result[10].null());
  BOOST_CHECK(result[11].null());
  BOOST_CHECK(result[12].null());
  BOOST_CHECK(result[13].null());
}

BOOST_AUTO_TEST_CASE( test_gcd_GF2_infinite_loop )
{
  FFElem<2> zero(0);
  FFElem<2> one(1);
  Poly<FFElem<2>> unit({one});
  // 73AF = 111001110101111 
  Poly<FFElem<2>> a({one, one, one, one, zero, one, zero, one, one, one, zero, zero, one, one, one});
  // EDF(1): X^11 + X^9 + X^7 + X^5 + X^4 + X^3 + X
  Poly<FFElem<2>> b({one, zero, zero, one, one, one, zero, one, zero, one, zero, one});

  BOOST_CHECK_EQUAL(unit, gcd(a,b));
}

BOOST_AUTO_TEST_CASE( test_cantor_zassenhaus_on_nintendo_1 )
{
  srand(time(NULL));
  FFElem<2> zero(0);
  FFElem<2> one(1);
  // 73AF = 111001110101111 
  // EDF(1): X^11 + X^9 + X^7 + X^5 + X^4 + X^3 + X
  // Poly<FFElem<2>> a({one, zero, zero, one, one, one, zero, one, zero, one, zero, one});
  Poly<FFElem<2>> a({one, one, one, one, zero, one, zero, one, one, one, zero, zero, one, one, one});
  // 83 = 10000011
  Poly<FFElem<2>> f1({one, one, zero, zero, zero, zero, zero, one});
  // E5 = 11100101
  Poly<FFElem<2>> f2({one, zero, one, zero, zero, one, one, one});

  auto result = a.cantorZassenhaus(2);

  //BOOST_CHECK_EQUAL(result.size(), 2);
  // result[0] and result[1] could be interverted ...
  BOOST_CHECK(result == f1 || result == f2);
  //BOOST_CHECK_EQUAL(result[0].second, 1);
  //BOOST_CHECK_EQUAL(result[1].first, f2);
  //BOOST_CHECK_EQUAL(result[1].second, 1);
}

// Example from: http://blog.fkraiem.org/2013/11/30/polynomial-factorisation-over-finite-fields-part-2-distinct-degree-factorisation/
BOOST_AUTO_TEST_CASE( test_distinct_degree_factors ) 
{
  FFElem<5> zero(0);
  FFElem<5> one(1);
  FFElem<5> two(2);
  FFElem<5> three(3);
  FFElem<5> four(4);
  Poly<FFElem<5>> f1({one, one});
  Poly<FFElem<5>> f2({two, one});
  Poly<FFElem<5>> f3({one, one, one});
  Poly<FFElem<5>> f4({two, one, one});

  Poly<FFElem<5>> f = f1*f2*f3*f4;
  Poly<FFElem<5>> g1({two, three, one});
  Poly<FFElem<5>> g2({two, three, four, two, one});
  Poly<FFElem<5>> nullPoly;

  auto result = f.distinctDegreeFactors(5);

  BOOST_CHECK_EQUAL(result.size(), 6);
  BOOST_CHECK_EQUAL(result[0], g1);
  BOOST_CHECK_EQUAL(result[1], g2);
  BOOST_CHECK_EQUAL(result[2], nullPoly);
  BOOST_CHECK_EQUAL(result[3], nullPoly);
  BOOST_CHECK_EQUAL(result[4], nullPoly);
  BOOST_CHECK_EQUAL(result[5], nullPoly);
}

BOOST_AUTO_TEST_SUITE_END()
