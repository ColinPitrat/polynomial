#define BOOST_TEST_MODULE TestFiniteFields
#include <boost/test/included/unit_test.hpp>
#include "finite_fields.h"
#include <sstream>

struct FiniteFieldsFixture
{
};

BOOST_FIXTURE_TEST_SUITE(TestFiniteFields, FiniteFieldsFixture)

BOOST_AUTO_TEST_CASE( test_default_constructed_is_0 ) 
{
  FFElem<5> p;

  BOOST_CHECK_EQUAL(p, FFElem<5>(0));
}

BOOST_AUTO_TEST_CASE( test_constructed_7_mod_5_is_2 ) 
{
  FFElem<5> p(7);

  BOOST_CHECK_EQUAL(p, FFElem<5>(2));
}

BOOST_AUTO_TEST_CASE( test_1_less_than_2 ) 
{
  FFElem<5> one(1);
  FFElem<5> two(2);

  BOOST_CHECK(one < two);
  BOOST_CHECK(one <= two);
  BOOST_CHECK(!(one >= two));
  BOOST_CHECK(!(one > two));
  BOOST_CHECK(!(one == two));
  BOOST_CHECK(one != two);
}

BOOST_AUTO_TEST_CASE( test_2_not_less_than_2 ) 
{
  FFElem<5> two(2);

  BOOST_CHECK(two == two);
  BOOST_CHECK(two <= two);
  BOOST_CHECK(two >= two);
  BOOST_CHECK(!(two < two));
  BOOST_CHECK(!(two > two));
}

BOOST_AUTO_TEST_CASE( test_3_not_less_than_2 ) 
{
  FFElem<5> three(3);
  FFElem<5> two(2);

  BOOST_CHECK(three != two);
  BOOST_CHECK(!(three == two));
  BOOST_CHECK(!(three < two));
  BOOST_CHECK(!(three <= two));
  BOOST_CHECK(three >= two);
  BOOST_CHECK(three > two);
}

BOOST_AUTO_TEST_CASE( test_2_plus_2_mod_3 ) 
{
  FFElem<3> a(2);
  FFElem<3> b(1);

  BOOST_CHECK_EQUAL(a+a, b);
}

BOOST_AUTO_TEST_CASE( test_2_minus_6_mod_7 ) 
{
  FFElem<7> a(2);
  FFElem<7> b(6);
  FFElem<7> c(3);

  BOOST_CHECK_EQUAL(a-b, c);
}

BOOST_AUTO_TEST_CASE( test_4_times_3_mod_5 ) 
{
  FFElem<5> a(4);
  FFElem<5> b(3);

  BOOST_CHECK_EQUAL(a*b, FFElem<5>(2));
}

BOOST_AUTO_TEST_CASE( test_4_over_3_mod_5 ) 
{
  FFElem<5> a(4);
  FFElem<5> b(3);

  BOOST_CHECK_EQUAL(a/b, FFElem<5>(3));
}

BOOST_AUTO_TEST_CASE( test_to_stream ) 
{
  std::ostringstream oss;
  FFElem<5> a(4);

  oss << a;

  BOOST_CHECK_EQUAL(oss.str(), "4[5]");
}

BOOST_AUTO_TEST_CASE( test_inverse_1_is_1 ) 
{
  FFElem<5> a(1);

  BOOST_CHECK_EQUAL(a.inverse(), a);
}

BOOST_AUTO_TEST_CASE( test_inverse_2_is_3 ) 
{
  FFElem<5> a(2);
  FFElem<5> b(3);

  BOOST_CHECK_EQUAL(a.inverse(), b);
}

BOOST_AUTO_TEST_CASE( test_inverse_3_is_2 ) 
{
  FFElem<5> a(3);
  FFElem<5> b(2);

  BOOST_CHECK_EQUAL(a.inverse(), b);
}

BOOST_AUTO_TEST_CASE( test_inverse_4_is_4 ) 
{
  FFElem<5> a(4);

  BOOST_CHECK_EQUAL(a.inverse(), a);
}

BOOST_AUTO_TEST_SUITE_END()
