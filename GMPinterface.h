#include <iostream>
#include <cmath>
#include <vector>
#include <gmp.h>
#include <gmpxx.h>
#include <stdlib.h>
#include <time.h>

#define mpz mpz_class
#define mpf mpf_class

using namespace std;

struct primepower {
    mpz a;
    unsigned int n;

    primepower(mpz aa, unsigned int nn) : a { aa }, n { nn } {}
};

int is_prop_prime(const mpz& a);
void printprimefactors(const vector<primepower>& v);
vector<mpz> primefactor(mpz a);
mpz round(const mpf_class& a);
mpz pow(const mpz& a, const mpz& b);

vector<primepower> primepowerpairs(const vector<mpz>& v);
mpz gcd(mpz a, mpz b);
mpz Eulerphi(mpz n);
string decToBin(mpz a);
mpz modularExponentiation(mpz a, mpz b, mpz n);
int JacobiSymbol(mpz a, mpz b);
mpz modularSqrt(mpz a, mpz p);
