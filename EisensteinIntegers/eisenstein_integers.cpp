#include "eisenstein_integers.h"

void printEIvector(const vector<EI>& v)
// prints a vector of EIs in '{}'
{
    cout << "{";
    for (int i = 0; i < v.size(); ++i) {
        if (i != v.size() - 1)
            cout << v[i] << ", ";
        else
            cout << v[i];
    }
    cout << "}";
}

// operator and member functions
// -----------------------------------------------------------------------------------

EI conj(const EI& alpha)
{
    return EI(alpha.get_a() - alpha.get_b(), -alpha.get_b());
}

EI quo(const EI& alpha, const EI& beta)
{
    EI prod = alpha * conj(beta);
    mpf p = mpf(prod.get_a())/beta.norm();
    mpf q = mpf(prod.get_b())/beta.norm();
    EI gamma {round(p), round(q)};
    return gamma;
}

EI rem(const EI& alpha, const EI& beta)
{
    return alpha % beta;
}

EI EI::operator / (const EI& alpha) const
{
    return quo(*this, alpha);
}

EI EI::operator % (const EI& alpha) const
{
    EI gamma = quo(*this, alpha);       // determine the quotient gamma
    EI rho = *this - alpha*gamma;       // calculate remainder
    return rho;
}

vector<EI> associates(const EI& alpha)
{
    vector<EI> assoc = { alpha, EI(-1)*alpha, EI(0, 1)*alpha, EI(0, -1)*alpha,
                        EI(1, 1)*alpha, EI(-1, -1)*alpha };
    return assoc;
}

bool operator == (const EI& alpha, const EI& beta)
{
    return (alpha.get_a() == beta.get_a()) && (alpha.get_b() == beta.get_b());
}

bool operator != (const EI& alpha, const EI& beta)
{
    return !(alpha == beta);
}

ostream& operator << (ostream& os, const EI& alpha)
{
    string s = to_string(alpha);
    return os << s;
}

// utilities
// -----------------------------------------------------------------------------------

string to_string(const EI& alpha)
// converts an Eisenstein integer to a string of the form a + b*w
{
    string s;
    if (alpha.get_a() != 0) {
        s += alpha.get_a().get_str();    // we add a to the string
        if (alpha.get_b() != 0) {        // three cases for the b part
            if (alpha.get_b() == 1)
                s += "+\xF9";
            else if (alpha.get_b() == -1)
                s += "-\xF9";
            else if (alpha.get_b() < 0) {
                mpz temp = abs(alpha.get_b());
                s += "-" + temp.get_str() + "\xF9";
            }
            else
                s += "+" + alpha.get_b().get_str() + "\xF9";
        }
    }
    else {
        if (alpha.get_b() == 0)
            s = "0";
        if (alpha.get_b() == 1)
            s = "\xF9";
        if (alpha.get_b() == -1)
            s = "-\xF9";
        if (alpha.get_b() != 1 && alpha.get_b() != -1 && alpha.get_b() != 0)
            s = alpha.get_b().get_str() + "\xF9";
    }
    return s;
}

void printEIproduct(const vector<EI>& v)
// write out a product of EIs with parentheses and powers
{
    if (v.size() == 0) {
        cout << "Error, empty vector" << endl;
        return;
    }
    vector<EI> temp = v;
    sort(temp.begin(), temp.end());

    int n = 1;
    cout << "(" << temp[0] << ")";
    for (int i = 1; i < temp.size(); ++i) {
        if (temp[i] != temp[i - 1]) {
            if (n > 1)
                cout << "^" << n;
            n = 1;
            cout << "(" << temp[i] << ")";
        }
        else
            ++n;
    }
    if (n > 1)
        cout << "^" << n;
    cout << endl;
}

EI EIproduct(const vector<EI>& v)
// iteratively compute a product of Eisenstein integers
{
    EI res = EI(1);
    for (EI alpha : v)
        res = res * alpha;
    return res;
}

EI EIpow(const EI& alpha, mpz n)
// calculates alpha^n
{
    if (alpha == EI(0) && n == 0) return EI(0);     // we define 0^0 = 0
    string bin = decToBin(n);                       // get the power n in binary
    EI res = EI(1);
    for (int i = bin.size() - 1; i >= 0; --i) {     // use repeated squarings
        res = res*res;
        if (bin[i] == '1') res = res*alpha;
    }
    return res;
}

EI EImodularExponentiation(EI alpha, mpz n, EI beta)
// calculates alpha^n mod beta, where n is a positive integer
{
    if (alpha == EI(0) && n == 0) return EI(0); // we define 0^0 = 0
    string bin = decToBin(n);                   // get the power n in binary
    EI res = EI(1);
    for (int i = bin.size() - 1; i >= 0; --i) { // use repeated squarings
        res = res*res % beta;
        if (bin[i] == '1') res = res*alpha % beta;
    }
    return res;
}

EI EIprimary(EI alpha)
// returns the unique primary associate of an odd Eisenstein integer alpha with norm =/= 3
{
    if (alpha.norm() % 3 == 0)
        throw runtime_error("input must not be divisible by 1 - w");
    if (alpha.get_a() % 3 == 0) {
        if (alpha.get_b() % 3 == 1 || alpha.get_b() % 3 == -2) return EI(1,1)*alpha;
        if (alpha.get_b() % 3 == 2 || alpha.get_b() % 3 == -1) return EI(-1,-1)*alpha;
    }
    if (alpha.get_a() % 3 == 1 || alpha.get_a() % 3 == -2) {
        if (alpha.get_b() % 3 == 0) return EI(-1)*alpha;
        if (alpha.get_b() % 3 == 1 || alpha.get_b() % 3 == -2) return EI(0,1)*alpha;
    }
    if (alpha.get_a() % 3 == 2 || alpha.get_a() % 3 == -1) {
        if (alpha.get_b() % 3 == 0) return alpha;
        if (alpha.get_b() % 3 == 2 || alpha.get_b() % 3 == -1) return EI(0,-1)*alpha;
    }
    cout << "FEJL!";

}

// primes and algorithms
// -----------------------------------------------------------------------------------

bool isEIprime(const EI& alpha)
// check if alpha is an Eisenstein prime using the characterization
// note that the function uses a probabilistic primality test, so it may in rare cases
// return false even if alpha is an Eisenstein prime
{
    if (is_prop_prime(alpha.norm()) == 2) return true;
    if (alpha.get_a() != 0 && alpha.get_b() == 0)
        if (abs(alpha.get_a()) % 3 == 2 && is_prop_prime(alpha.get_a()) == 2) return true;
    if (alpha.get_a() == 0 && alpha.get_b() != 0)
        if (abs(alpha.get_b()) % 3 == 2 && is_prop_prime(alpha.get_b()) == 2) return true;
    return false;
}

EI EIgcd(const EI& alpha, const EI& beta)
// Euclidean algorithm for Eisenstein integers
{
    if (alpha % beta == EI(0))
        return beta;
    return EIgcd(beta, alpha % beta);
}

EI EIgcd(const vector<EI>& EIs)
// the gcd of a vector of Eisenstein integers
{
    if (EIs.size() == 0)
        return EI(0);
    if (EIs.size() == 1)
        return EIs[0];

    EI gcd = EIgcd(EIs[0], EIs[1]);
    for (int i = 2; i < EIs.size(); ++i) {
        gcd = EIgcd(gcd, EIs[i]);
    }
    return gcd;
}

bool EIpairwisecoprime(const vector<EI>& EIs)
// returns true, if the GIs in the vector are pairwise coprime
{
    for (int i = 0; i < EIs.size(); ++i) {
        for (int j = i + 1; j < EIs.size(); ++j) {
            if (EIgcd(EIs[i], EIs[j]).norm() != 1) return false;
        }
    }
    return true;
}

vector<EI> EIexgcd(const EI& alpha, const EI& beta)
// Extended Euclidean algorithm for Eisenstein integers
// Returns {gcd, s, t} where s and t are coefficients satisfying alpha*s + beta*t = gcd
{
    if (beta == EI(0)) {
        vector<EI> res = {alpha, EI(1), EI(0)};
        return res;
    }
    else {
        vector<EI> temp = EIexgcd(beta, alpha % beta);
        vector<EI> res = {temp[0], temp[2], temp[1] - quo(alpha, beta)*temp[2]};
        return res;
    }
}

EI EIlcm(const EI& alpha, const EI& beta)
// calculate the lcm of two Eisenstein integers using the Euclidean algorithm
{
    return (alpha * beta) / EIgcd(alpha, beta);
}

EI EIlcm(const vector<EI>& EIs)
// calculate the lcm of a vector of Eisenstein integers
{
    if (EIs.size() == 0)
        return EI(0);
    if (EIs.size() == 1)
        return EIs[0];

    EI lcm = EIlcm(EIs[0], EIs[1]);
    for (int i = 2; i < EIs.size(); ++i) {
        lcm = EIlcm(lcm, EIs[i]);
    }
    return lcm;
}

EI EImodlinearsolve(const EI& alpha, const EI& beta, const EI& gamma)
// returns a solution rho (if one exists) to the equation alpha*rho = beta
// (mod gamma)
{
    vector<EI> s = EIexgcd(alpha, gamma);
    if (beta % s[0] == EI(0)) {
        EI rho = (s[1] * beta)/s[0];
        rho = rho % gamma;
        return rho;
    }
    cout << "No solution to equation: (" << alpha << ")rho" << " = " << beta << " (mod " << gamma << ")" << endl;
    return EI(0);
}

EI EIChRem(const vector<EI> alpha, const vector<EI> moduli)
// The Chinese Remainder Theorem for Z[w], finds rho such that
// alpha[0] = rho mod moduli[0], ..., alpha[n] = rho mod moduli[n] where
// n = alpha.size() = moduli.size()
{
    if (alpha.size() != moduli.size()) {
        cout << "Error, input lengths are not identical" << endl;
        return EI(0);
    }
    if (alpha.size() == 0)
        return EI(0);

    EI rho = alpha[0];      // the solution rho
    EI gamma = moduli[0];   // the product of the moduli
    if (EIpairwisecoprime(moduli)) {
        for (int i = 1; i < alpha.size(); ++i) {
            rho = rho + gamma * EImodlinearsolve(gamma, alpha[i] - rho, moduli[i]);
            gamma = gamma * moduli[i];
        }
        return rho;
    }
    cout << "Error, moduli must be pairwise coprime" << endl;
    return EI(0);
}

vector<EI> EIprimefactor(EI alpha)
// factors alpha as a product of irreducibles and a single unit
{
    vector<mpz> p = primefactor(alpha.norm());
    vector<EI> factors;
    for (int i = 0; i < p.size(); ++i) {
        if (p[i] == 3) {                // if p[i] = 3, 1 - w is a factor
            factors.push_back(EI(1,-1));
            alpha = alpha / EI(1, -1);
        }
        else if (p[i] % 3 == 2) {       // if p[i] % 3 == 2, remove two factors of p[i]
            factors.push_back(p[i]);
            alpha = alpha / p[i];
            ++i;
        }
        else if (p[i] % 3 == 1) {
            mpz k = modularSqrt(-3, p[i]) % p[i];  // find k, such that k^2 = -3 mod p[i] (and reduce mod p[i])
            EI factor = EIgcd(p[i], {1 + k, 2});
            if (alpha % factor != EI(0)) factor = conj(factor); // if the factor does not divide, the conjugate is the desired factor
            factors.push_back(factor);
            alpha = alpha / factor;
        }
    }
    factors.push_back(alpha);       // when the previous loop is done, alpha is a unit
    return factors;
}

EI CubicRes(EI alpha, EI pi)
// uses the definition to compute the residue symbol (alpha, pi)_3 using modular exponentiation
// pi is assumed to be prime
{
    if (pi % EI(1,-1) == EI(0))
        throw runtime_error("second parameter cannot be divisible by 1 - w");
    return EImodularExponentiation(alpha, (pi.norm() - 1)/3, pi);
}

EI CubicResSymbol(EI alpha, EI beta)
// calculates the cubic residue symbol (alpha, beta)_3 where beta is an Eisenstein integer
// not divisible by 1 - w
{
    if (beta % EI(1,-1) == EI(0))
        throw runtime_error("second parameter must not be divisible by 1 - w");
    EI res = EI(1);

    while (true) {
        beta = EIprimary(beta);
        alpha %= beta;

        if (alpha == EI(0)) {
            if (beta.norm() != 1) return EI(0);     // alpha and beta have a common factor
            else return res;
        }

        while (alpha % EI(1, -1) == EI(0)) {        // remove factors of 1 - w and apply
            alpha /= EI(1, -1);                     // the supplementary law for 1 - w
            res *= EIpow(EI(0,1), (2*(beta.get_a() + 1)/3) % 3 + 3);
        }

        EI u = alpha / EIprimary(alpha);            // the inverse of the unit for the primary associate
        if (u == EI(0,1) || u == EI(0, -1)) {
            if (beta.norm() % 9 == 4) res *= EI(0,1);
            if (beta.norm() % 9 == 7) res *= EI(-1,-1);
        }
        if (u == EI(1,1) || u == EI(-1,-1)) {
            if (beta.norm() % 9 == 4) res *= EI(-1,-1);
            if (beta.norm() % 9 == 7) res *= EI(0,1);
        }
        alpha = EIprimary(alpha);
        swap(alpha, beta);
    }
}
