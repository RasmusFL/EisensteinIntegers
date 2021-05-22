#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>
#include "../GMPinterface.h"

using namespace std;

struct EI;

void printEIvector(const vector<EI>& v);

class EI
    // Eisenstein integer
{
    // coefficients a and b in a + bw
    mpz a;
    mpz b;
    mpz N = a*a + b*b - a*b;
public:

    // constructors
    EI() : a{ 0 }, b{ 0 } {}
    EI(mpz aa) : a{ aa }, b{ 0 } {}
    EI(mpz aa, mpz bb) : a{ aa }, b{ bb } {}

    // get member values, real and imaginary part, add get for a and b!
    const mpz get_a() const { return a; }
    const mpz get_b() const { return b; }
    const mpz norm() const { return N; }

    // set member values
    void set_a(const mpz& aa) { a = aa; N = a*a + b*b - a*b; }
    void set_b(const mpz& bb) { b = bb; N = a*a + b*b - a*b; }

    // operator overloads
    EI operator + (const EI& alpha) const { return EI(a + alpha.a, b + alpha.b); }
    EI operator - (const EI& alpha) const { return EI(a - alpha.a, b - alpha.b); }
    EI operator * (const EI& alpha) const { return EI(a*alpha.a - b*alpha.b, a*alpha.b + b*alpha.a - b*alpha.b); }
    EI operator / (const EI& alpha) const;      // calculates the quotient from the Euclidean division
    EI operator % (const EI& alpha) const;
    EI& operator += (const EI& alpha) { *this = *this + alpha; return *this; }
    EI& operator -= (const EI& alpha) { *this = *this - alpha; return *this; }
    EI& operator *= (const EI& alpha) { *this = *this * alpha; return *this; }
    EI& operator /= (const EI& alpha) { *this = *this / alpha; return *this; }
    EI& operator %= (const EI& alpha) { *this = *this % alpha; return *this; }
    bool operator < (const EI& alpha) const { return N < alpha.norm(); }     // for sorting by norm
};

// operator functions
EI conj(const EI& alpha);
EI quo(const EI& alpha, const EI& beta);
EI rem(const EI& alpha, const EI& beta);
vector<EI> associates(const EI& alpha);
bool operator == (const EI& alpha, const EI& beta);
bool operator != (const EI& alpha, const EI& beta);
ostream& operator << (ostream& os, const EI& alpha);

// utilities
string to_string(const EI& alpha);
struct EIprimepower { EI alpha; unsigned int n; };
void printEIproduct(const vector<EI>& v);
EI EIproduct(const vector<EI>& v);
EI EIpow(const EI& alpha, mpz n);
EI EImodularExponentiation(EI alpha, mpz n, EI beta);
EI EIprimary(EI alpha);

// primes and algorithms
bool isEIprime(const EI& alpha);
EI EIgcd(const EI& alpha, const EI& beta);
EI EIgcd(const vector<EI>& EIs);
bool EIpairwisecoprime(const vector<EI>& EIs);
vector<EI> EIexgcd(const EI& alpha, const EI& beta);
EI EIlcm(const EI& alpha, const EI& beta);
EI EIlcm(const vector<EI>& EIs);
EI EImodlinearsolve(const EI& alpha, const EI& beta, const EI& gamma);
EI EIChRem(const vector<EI> alpha, const vector<EI> moduli);
vector<EI> EIprimefactor(EI alpha);
EI CubicRes(EI alpha, EI pi);
EI CubicResSymbol(EI alpha, EI beta);
