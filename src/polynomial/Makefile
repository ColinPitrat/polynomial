all:
	g++ -O2 -Wall -Werror -std=c++11 test_polynomial.cc -o test_polynomial -lgmp
	g++ -O2 -Wall -Werror -std=c++11 test_g2polynomial.cc -o test_g2polynomial
	./test_polynomial --log_level=test_suite
	./test_g2polynomial --log_level=test_suite

clean:
	rm -f test_polynomial test_g2polynomial
