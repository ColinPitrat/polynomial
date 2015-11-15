#define BOOST_TEST_MODULE TestG2Polynomial
#include <boost/test/included/unit_test.hpp>

// g2polynomial and g2polynomial_bitset provide the exact same class with the same interface but different implementations
// This test can be used to test both
//#include "g2polynomial.h"
#include "g2polynomial_bitset.h"

#include <sstream>
#include <algorithm>
#include <utility>
#include <cstdlib>
#include <iostream>

struct G2PolynomialFixture
{
};

BOOST_FIXTURE_TEST_SUITE(TestG2Polynomial, G2PolynomialFixture)

BOOST_AUTO_TEST_CASE( test_default_constructed_is_null_and_has_degree_minus_1 ) 
{
  G2Poly p;

  BOOST_CHECK_EQUAL(true, p.null());
  BOOST_CHECK_EQUAL(-1, p.degree());
}

BOOST_AUTO_TEST_CASE( test_constructed_X2_has_degree_2 ) 
{
  G2Poly p({0, 0, 1});

  BOOST_CHECK_EQUAL(2, p.degree());
}

BOOST_AUTO_TEST_CASE( test_one_has_degree_0 ) 
{
  G2Poly p({1});

  BOOST_CHECK_EQUAL(0, p.degree());
}

BOOST_AUTO_TEST_CASE( test_X0_has_degree_0 ) 
{
  G2Poly x0 = G2Poly::Xn(0);

  BOOST_CHECK_EQUAL(0, x0.degree());
}

BOOST_AUTO_TEST_CASE( test_X5_has_degree_5 ) 
{
  G2Poly x5 = G2Poly::Xn(5);

  BOOST_CHECK_EQUAL(5, x5.degree());
}

BOOST_AUTO_TEST_CASE( test_constructed_X3_equal_X3 ) 
{
  G2Poly x3 = G2Poly::Xn(3);
  G2Poly cx3({0, 0, 0, 1});

  BOOST_CHECK_EQUAL(x3, cx3);
}

BOOST_AUTO_TEST_CASE( test_X5_not_equal_X3 ) 
{
  G2Poly x3 = G2Poly::Xn(3);
  G2Poly x5 = G2Poly::Xn(5);

  BOOST_CHECK(x3 != x5);
}

BOOST_AUTO_TEST_CASE( test_usual_sum ) 
{
  G2Poly p1({1, 1});
  G2Poly p2({1, 1});
  G2Poly result({0});

  BOOST_CHECK_EQUAL(result, p1+p2);
}

BOOST_AUTO_TEST_CASE( test_usual_sum2 ) 
{
  G2Poly p1({1, 0, 0, 1});
  G2Poly p2({0, 0, 1, 1});
  G2Poly result({1, 0, 1});

  BOOST_CHECK_EQUAL(result, p1+p2);
}

BOOST_AUTO_TEST_CASE( test_usual_diff ) 
{
  G2Poly p1({1, 1});
  G2Poly p2({1, 1});
  G2Poly result({0});

  BOOST_CHECK_EQUAL(result, p1-p2);
}

BOOST_AUTO_TEST_CASE( test_usual_diff2 ) 
{
  G2Poly p1({1, 0, 0, 1});
  G2Poly p2({0, 0, 1, 1});
  G2Poly result({1, 0, 1});

  BOOST_CHECK_EQUAL(result, p1-p2);
}

BOOST_AUTO_TEST_CASE( test_usual_product ) 
{
  G2Poly p1({1, 1});
  G2Poly p2({1, 1});
  G2Poly result({1, 0, 1});

  BOOST_CHECK_EQUAL(result, p1*p2);
}

BOOST_AUTO_TEST_CASE( test_usual_product2 ) 
{
  G2Poly p1({1, 0, 1});
  G2Poly p2({0, 1, 1});
  G2Poly result({0, 1, 1, 1, 1});

  BOOST_CHECK_EQUAL(result, p1*p2);
}

BOOST_AUTO_TEST_CASE( test_toStream_of_null ) 
{
  G2Poly p;
  std::ostringstream oss;

  oss << p;

  BOOST_CHECK_EQUAL(oss.str(), "0");
}

BOOST_AUTO_TEST_CASE( test_toStream_of_1 ) 
{
  G2Poly p({1});
  std::ostringstream oss;

  oss << p;

  BOOST_CHECK_EQUAL(oss.str(), "1");
}

BOOST_AUTO_TEST_CASE( test_toStream_of_X ) 
{
  G2Poly p({0, 1});
  std::ostringstream oss;

  oss << p;

  BOOST_CHECK_EQUAL(oss.str(), "X");
}

BOOST_AUTO_TEST_CASE( test_toStream_of_X2_plus_1 ) 
{
  G2Poly p({1, 0, 1});
  std::ostringstream oss;

  oss << p;

  BOOST_CHECK_EQUAL(oss.str(), "X^2 + 1");
}

BOOST_AUTO_TEST_CASE( test_toStream_of_some_polynomial ) 
{
  G2Poly p({0, 1, 1, 1, 0, 1});
  std::ostringstream oss;

  oss << p;

  BOOST_CHECK_EQUAL(oss.str(), "X^5 + X^3 + X^2 + X");
}

BOOST_AUTO_TEST_CASE( test_derivate_of_polynomial ) 
{
  G2Poly p({1, 1, 1, 1});
  G2Poly dp({1, 0, 1});

  BOOST_CHECK_EQUAL(p.derivate(), dp);
}

BOOST_AUTO_TEST_CASE( test_derivate_of_one ) 
{
  G2Poly p({1});
  G2Poly dp;

  BOOST_CHECK_EQUAL(p.derivate(), dp);
}

BOOST_AUTO_TEST_CASE( test_derivate_of_null ) 
{
  G2Poly p;
  G2Poly dp;

  BOOST_CHECK_EQUAL(p.derivate(), dp);
}

BOOST_AUTO_TEST_CASE( test_remainder ) 
{
  G2Poly a({1, 0, 1, 0, 1});
  G2Poly b({1, 0, 1});
  G2Poly r({1});

  BOOST_CHECK_EQUAL(r, a%b);
}

BOOST_AUTO_TEST_CASE( test_remainder2 ) 
{
  G2Poly a({1, 0, 0, 0, 1});
  G2Poly b({1, 1, 0, 1, 1});
  G2Poly r({0, 1, 0, 1});

  BOOST_CHECK_EQUAL(r, a%b);
}

BOOST_AUTO_TEST_CASE( test_divide_no_rem ) 
{
  G2Poly a({0, 0, 1, 1});
  G2Poly b({1, 1});
  G2Poly r({0, 0, 1});

  BOOST_CHECK_EQUAL(r, a/b);
}

BOOST_AUTO_TEST_CASE( test_divide_with_rem ) 
{
  G2Poly a({1, 0, 1, 1});
  G2Poly b({1, 1});
  G2Poly r({0, 0, 1});

  BOOST_CHECK_EQUAL(r, a/b);
}

BOOST_AUTO_TEST_CASE( test_gcd ) 
{
  G2Poly a({1, 0, 0, 0, 1});
  G2Poly b({1, 1, 0, 1, 1});
  G2Poly r({1, 0, 1});

  BOOST_CHECK_EQUAL(r, gcd(a,b));
}

BOOST_AUTO_TEST_CASE( test_square_free ) 
{
  G2Poly a({1, 0, 0, 0, 1});
  G2Poly r({1, 1});

  BOOST_CHECK_EQUAL(a.squareFreePart(2), r);
}

BOOST_AUTO_TEST_CASE( test_square_free_null_derivative ) 
{
  G2Poly a({1, 0, 1, 0, 1});
  G2Poly r({1, 1, 1});

  BOOST_CHECK_EQUAL(a.squareFreePart(2), r);
}

BOOST_AUTO_TEST_CASE( test_distinct_degree_factors_finite_field_GF2 ) 
{
  G2Poly f({1, 0, 0, 1, 0, 0, 1});

  auto result = f.distinctDegreeFactors(2);

  BOOST_CHECK_EQUAL(result.size(), 6);
}

BOOST_AUTO_TEST_CASE( test_square_free_factors_on_GF2 ) 
{
  G2Poly a({1, 0, 1});
  G2Poly f1({1, 1});

  auto result = a.squareFreeFactors(2);

  BOOST_CHECK_EQUAL(result.size(), 1);
  BOOST_CHECK_EQUAL(result[0].first, f1);
  BOOST_CHECK_EQUAL(result[0].second, 2);
}

/*
BOOST_AUTO_TEST_CASE( test_distinct_degree_factors_on_nintendo_1 )
{
  // 73AF = 111001110101111 
  G2Poly a({1, 1, 1, 1, 0, 1, 0, 1, 1, 1, 0, 0, 1, 1, 1});
  G2Poly unit({1});
  // 83 = 10000011
  G2Poly f1({1, 1, 0, 0, 0, 0, 0, 1});
  // E5 = 11100101
  G2Poly f2({1, 0, 1, 0, 0, 1, 1, 1});

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
*/

/*
BOOST_AUTO_TEST_CASE( test_distinct_degree_factors_on_nintendo_2 )
{
  // 738377C1 = 1110011100000110111011111000001
  G2Poly a({1, 1, 1, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 1, 1, 0, 1, 1, 1, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 1});
  G2Poly unit({1});
  // CD55 = 1100110101010101
  G2Poly f1({1, 1, 0, 0, 1, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1});
  // B0C5 = 1011000011000101
  G2Poly f2({1, 0, 1, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1, 0, 1});

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
*/

BOOST_AUTO_TEST_CASE( test_gcd_GF2_infinite_loop )
{
  G2Poly unit({1});
  // 73AF = 111001110101111 
  G2Poly a({1, 1, 1, 1, 0, 1, 0, 1, 1, 1, 0, 0, 1, 1, 1});
  // EDF(1): X^11 + X^9 + X^7 + X^5 + X^4 + X^3 + X
  G2Poly b({1, 0, 0, 1, 1, 1, 0, 1, 0, 1, 0, 1});

  BOOST_CHECK_EQUAL(unit, gcd(a,b));
}

BOOST_AUTO_TEST_CASE( test_mc_eliece )
{
  // Example from: http://www.ams.org/journals/mcom/1969-23-108/S0025-5718-1969-0257039-X/S0025-5718-1969-0257039-X.pdf
  G2Poly a({1, 1, 0, 0, 1, 1, 0, 1});
  G2Poly f1({1, 0, 1, 1, 1, 1});
  G2Poly f2({1, 1, 1});

  auto result = a.mceliece(2);

  //BOOST_CHECK_EQUAL(result.size(), 2);
  // result[0] and result[1] could be interverted ...
  BOOST_CHECK(result == f1 || result == f2);
  //BOOST_CHECK_EQUAL(result[0].second, 1);
  //BOOST_CHECK_EQUAL(result[1].first, f2);
  //BOOST_CHECK_EQUAL(result[1].second, 1);
}

/*
BOOST_AUTO_TEST_CASE( test_knuth_on_nintendo_1 )
{
  srand(time(NULL));
  // E5 = 11100101
  G2Poly a({1, 1, 1, 1, 0, 1, 0, 1, 1, 1, 0, 0, 1, 1, 1});
  // 83 = 10000011
  G2Poly f1({1, 1, 0, 0, 0, 0, 0, 1});
  // E5 = 11100101
  G2Poly f2({1, 0, 1, 0, 0, 1, 1, 1});

  auto result = a.knuth(2);

  //BOOST_CHECK_EQUAL(result.size(), 2);
  // result[0] and result[1] could be interverted ...
  BOOST_CHECK(result == f1 || result == f2);
  //BOOST_CHECK_EQUAL(result[0].second, 1);
  //BOOST_CHECK_EQUAL(result[1].first, f2);
  //BOOST_CHECK_EQUAL(result[1].second, 1);
}

BOOST_AUTO_TEST_CASE( test_knuth_on_nintendo_2 )
{
  srand(time(NULL));
  // 738377C1 = 1110011100000110111011111000001
  G2Poly a({1, 1, 1, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 1, 1, 0, 1, 1, 1, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 1});
  // CD55 = 1100110101010101
  G2Poly f1({1, 1, 0, 0, 1, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1});
  // B0C5 = 1011000011000101
  G2Poly f2({1, 0, 1, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1, 0, 1});

  auto result = a.knuth(2);

  //BOOST_CHECK_EQUAL(result.size(), 2);
  // result[0] and result[1] could be interverted ...
  BOOST_CHECK(result == f1 || result == f2);
  //BOOST_CHECK_EQUAL(result[0].second, 1);
  //BOOST_CHECK_EQUAL(result[1].first, f2);
  //BOOST_CHECK_EQUAL(result[1].second, 1);
}

BOOST_AUTO_TEST_CASE( test_knuth_on_nintendo_3 )
{
  srand(time(NULL));
  std::cout << "Start test_knuth_on_nintendo_3" << std::endl;
  // 6677E20146508FB7 = 110011001110111111000100000000101000110010100001000111110110111
  G2Poly a({1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 1, 1, 0, 0, 1, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 1, 1, 1, 1, 1, 0, 1, 1, 0, 1, 1, 1});
  // B0C152F9 = 10110000110000010101001011111001
  G2Poly f1({1, 0, 1, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0, 0, 1, 0, 1, 1, 1, 1, 1, 0, 0, 1});
  // EBF2831F = 11101011111100101000001100011111
  G2Poly f2({1, 1, 1, 0, 1, 0, 1, 1, 1, 1, 1, 1, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1, 1, 1, 1, 1});

  auto result = a.knuth(2);

  //BOOST_CHECK_EQUAL(result.size(), 2);
  // result[0] and result[1] could be interverted ...
  BOOST_CHECK(result == f1 || result == f2);
  //BOOST_CHECK_EQUAL(result[0].second, 1);
  //BOOST_CHECK_EQUAL(result[1].first, f2);
  //BOOST_CHECK_EQUAL(result[1].second, 1);
  std::cout << "End test_knuth_on_nintendo_3" << std::endl;
}
*/

BOOST_AUTO_TEST_CASE( test_mc_eliece_on_nintendo_1 )
{
  // E5 = 11100101
  G2Poly a({1, 1, 1, 1, 0, 1, 0, 1, 1, 1, 0, 0, 1, 1, 1});
  // 83 = 10000011
  G2Poly f1({1, 1, 0, 0, 0, 0, 0, 1});
  // E5 = 11100101
  G2Poly f2({1, 0, 1, 0, 0, 1, 1, 1});

  auto result = a.mceliece(2);

  //BOOST_CHECK_EQUAL(result.size(), 2);
  // result[0] and result[1] could be interverted ...
  BOOST_CHECK(result == f1 || result == f2);
  //BOOST_CHECK_EQUAL(result[0].second, 1);
  //BOOST_CHECK_EQUAL(result[1].first, f2);
  //BOOST_CHECK_EQUAL(result[1].second, 1);
}

BOOST_AUTO_TEST_CASE( test_mc_eliece_on_nintendo_2 )
{
  // 738377C1 = 1110011100000110111011111000001
  G2Poly a({1, 1, 1, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 1, 1, 0, 1, 1, 1, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 1});
  // CD55 = 1100110101010101
  G2Poly f1({1, 1, 0, 0, 1, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1});
  // B0C5 = 1011000011000101
  G2Poly f2({1, 0, 1, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1, 0, 1});

  auto result = a.mceliece(2);

  //BOOST_CHECK_EQUAL(result.size(), 2);
  // result[0] and result[1] could be interverted ...
  BOOST_CHECK(result == f1 || result == f2);
  //BOOST_CHECK_EQUAL(result[0].second, 1);
  //BOOST_CHECK_EQUAL(result[1].first, f2);
  //BOOST_CHECK_EQUAL(result[1].second, 1);
}

BOOST_AUTO_TEST_CASE( test_mc_eliece_on_nintendo_3 )
{
  std::cout << "Start test_mc_eliece_on_nintendo_3" << std::endl;
  // 6677E20146508FB7 = 110011001110111111000100000000101000110010100001000111110110111
  G2Poly a({1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 1, 1, 0, 0, 1, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 1, 1, 1, 1, 1, 0, 1, 1, 0, 1, 1, 1});
  // B0C152F9 = 10110000110000010101001011111001
  G2Poly f1({1, 0, 1, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0, 0, 1, 0, 1, 1, 1, 1, 1, 0, 0, 1});
  // EBF2831F = 11101011111100101000001100011111
  G2Poly f2({1, 1, 1, 0, 1, 0, 1, 1, 1, 1, 1, 1, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1, 1, 1, 1, 1});

  auto result = a.mceliece(2);

  //BOOST_CHECK_EQUAL(result.size(), 2);
  // result[0] and result[1] could be interverted ...
  BOOST_CHECK(result == f1 || result == f2);
  //BOOST_CHECK_EQUAL(result[0].second, 1);
  //BOOST_CHECK_EQUAL(result[1].first, f2);
  //BOOST_CHECK_EQUAL(result[1].second, 1);
  std::cout << "End test_mc_eliece_on_nintendo_3" << std::endl;
}

BOOST_AUTO_TEST_CASE( test_mc_eliece_on_nintendo_6 )
{
  std::cout << "Start test_mc_eliece_on_nintendo_6" << std::endl;
  // 6C38C43D39E141ADE270BA243FADE25323C685535EA02B9CBB74360A4A3BB44642AD10FB14A5E47CCA60550F0906E6D0C72F6A8B465C5267390293804AF6FC33 
  G2Poly a({1, 1, 0, 1, 1, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 1, 1, 1, 1, 0, 1, 0, 0, 1, 1, 1, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 1, 0, 1, 0, 1, 1, 0, 1, 1, 1, 1, 0, 0, 0, 1, 0, 0, 1, 1, 1, 0, 0, 0, 0, 1, 0, 1, 1, 1, 0, 1, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 0, 1, 0, 1, 1, 0, 1, 1, 1, 1, 0, 0, 0, 1, 0, 0, 1, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 1, 1, 0, 1, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0, 1, 0, 0, 1, 1, 0, 1, 0, 1, 1, 1, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 1, 1, 0, 0, 1, 1, 1, 0, 0, 1, 0, 1, 1, 1, 0, 1, 1, 0, 1, 1, 1, 0, 1, 0, 0, 0, 0, 1, 1, 0, 1, 1, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 1, 0, 0, 1, 0, 1, 0, 0, 0, 1, 1, 1, 0, 1, 1, 1, 0, 1, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0, 1, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 1, 1, 0, 0, 0, 1, 0, 1, 0, 0, 1, 0, 1, 0, 0, 1, 0, 1, 1, 1, 1, 0, 0, 1, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 1, 1, 0, 0, 1, 0, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 1, 1, 0, 1, 1, 1, 0, 0, 1, 1, 0, 1, 1, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1, 1, 1, 0, 0, 1, 0, 1, 1, 1, 1, 0, 1, 1, 0, 1, 0, 1, 0, 1, 0, 0, 0, 1, 0, 1, 1, 0, 1, 0, 0, 0, 1, 1, 0, 0, 1, 0, 1, 1, 1, 0, 0, 0, 1, 0, 1, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 1, 1, 0, 0, 1, 1, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 1, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 1, 0, 1, 1, 1, 1, 0, 1, 1, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 0, 0, 1, 1});
  // B21251C47616E9D6B8FB36D373B746C5EF4D0EF5DC5035696CC9659D33749A7F
  G2Poly f1({1, 0, 1, 1, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 1, 0, 0, 0, 1, 1, 1, 0, 0, 0, 1, 0, 0, 0, 1, 1, 1, 0, 1, 1, 0, 0, 0, 0, 1, 0, 1, 1, 0, 1, 1, 1, 0, 1, 0, 0, 1, 1, 1, 0, 1, 0, 1, 1, 0, 1, 0, 1, 1, 1, 0, 0, 0, 1, 1, 1, 1, 1, 0, 1, 1, 0, 0, 1, 1, 0, 1, 1, 0, 1, 1, 0, 1, 0, 0, 1, 1, 0, 1, 1, 1, 0, 0, 1, 1, 1, 0, 1, 1, 0, 1, 1, 1, 0, 1, 0, 0, 0, 1, 1, 0, 1, 1, 0, 0, 0, 1, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 0, 1, 0, 0, 1, 1, 0, 1, 0, 0, 0, 0, 1, 1, 1, 0, 1, 1, 1, 1, 0, 1, 0, 1, 1, 1, 0, 1, 1, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 1, 0, 1, 0, 1, 1, 0, 1, 0, 0, 1, 0, 1, 1, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 0, 0, 1, 0, 1, 1, 0, 0, 1, 0, 1, 1, 0, 0, 1, 1, 1, 0, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 1, 1, 1, 0, 1, 0, 0, 1, 0, 0, 1, 1, 0, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1});
  // FA595299956FF31D84CC2EFC1A4FF1070AFE2A411AAAA61CB61B13F6320A18D5 
  G2Poly f2({1, 1, 1, 1, 1, 0, 1, 0, 0, 1, 0, 1, 1, 0, 0, 1, 0, 1, 0, 1, 0, 0, 1, 0, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 0, 1, 0, 1, 0, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 0, 0, 0, 1, 1, 1, 0, 1, 1, 0, 0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 1, 0, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 1, 1, 0, 1, 0, 0, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 1, 0, 1, 0, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 1, 0, 1, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 0, 1, 1, 0, 0, 0, 0, 1, 1, 1, 0, 0, 1, 0, 1, 1, 0, 1, 1, 0, 0, 0, 0, 1, 1, 0, 1, 1, 0, 0, 0, 1, 0, 0, 1, 1, 1, 1, 1, 1, 0, 1, 1, 0, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1, 1, 0, 1, 0, 1, 0, 1});

  auto result = a.mceliece(2);

  //BOOST_CHECK_EQUAL(result.size(), 2);
  // result[0] and result[1] could be interverted ...
  BOOST_CHECK(result == f1 || result == f2);
  //BOOST_CHECK_EQUAL(result[0].second, 1);
  //BOOST_CHECK_EQUAL(result[1].first, f2);
  //BOOST_CHECK_EQUAL(result[1].second, 1);
  std::cout << "End test_mc_eliece_on_nintendo_6" << std::endl;
}

/*
BOOST_AUTO_TEST_CASE( test_cantor_zassenhaus_on_nintendo_1 )
{
  srand(time(NULL));
  // 73AF = 111001110101111 
  // EDF(1): X^11 + X^9 + X^7 + X^5 + X^4 + X^3 + X
  // Poly<FFElem<2>> a({1, 0, 0, 1, 1, 1, 0, 1, 0, 1, 0, 1});
  G2Poly a({1, 1, 1, 1, 0, 1, 0, 1, 1, 1, 0, 0, 1, 1, 1});
  // 83 = 10000011
  G2Poly f1({1, 1, 0, 0, 0, 0, 0, 1});
  // E5 = 11100101
  G2Poly f2({1, 0, 1, 0, 0, 1, 1, 1});

  auto result = a.cantorZassenhaus(2);

  //BOOST_CHECK_EQUAL(result.size(), 2);
  // result[0] and result[1] could be interverted ...
  BOOST_CHECK(result == f1 || result == f2);
  //BOOST_CHECK_EQUAL(result[0].second, 1);
  //BOOST_CHECK_EQUAL(result[1].first, f2);
  //BOOST_CHECK_EQUAL(result[1].second, 1);
}
*/

/*
BOOST_AUTO_TEST_CASE( test_cantor_zassenhaus_on_nintendo_2 )
{
  srand(time(NULL));
  // 738377C1 = 1110011100000110111011111000001
  G2Poly a({1, 1, 1, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 1, 1, 0, 1, 1, 1, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 1});
  // CD55 = 1100110101010101
  G2Poly f1({1, 1, 0, 0, 1, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1});
  // B0C5 = 1011000011000101
  G2Poly f2({1, 0, 1, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1, 0, 1});

  auto result = a.cantorZassenhaus(2);

  //BOOST_CHECK_EQUAL(result.size(), 2);
  // result[0] and result[1] could be interverted ...
  BOOST_CHECK(result == f1 || result == f2);
  //BOOST_CHECK_EQUAL(result[0].second, 1);
  //BOOST_CHECK_EQUAL(result[1].first, f2);
  //BOOST_CHECK_EQUAL(result[1].second, 1);
}

BOOST_AUTO_TEST_CASE( test_cantor_zassenhaus_on_nintendo_3 )
{
  srand(time(NULL));
  // 6677E20146508FB7 = 110011001110111111000100000000101000110010100001000111110110111
  G2Poly a({1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 1, 1, 0, 0, 1, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 1, 1, 1, 1, 1, 0, 1, 1, 0, 1, 1, 1});
  // B0C152F9 = 10110000110000010101001011111001
  G2Poly f1({1, 0, 1, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0, 0, 1, 0, 1, 1, 1, 1, 1, 0, 0, 1});
  // EBF2831F = 11101011111100101000001100011111
  G2Poly f2({1, 1, 1, 0, 1, 0, 1, 1, 1, 1, 1, 1, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1, 1, 1, 1, 1});

  auto result = a.cantorZassenhaus(2);

  //BOOST_CHECK_EQUAL(result.size(), 2);
  // result[0] and result[1] could be interverted ...
  BOOST_CHECK(result == f1 || result == f2);
  //BOOST_CHECK_EQUAL(result[0].second, 1);
  //BOOST_CHECK_EQUAL(result[1].first, f2);
  //BOOST_CHECK_EQUAL(result[1].second, 1);
}
*/

BOOST_AUTO_TEST_SUITE_END()
