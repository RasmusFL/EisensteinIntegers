#include "GMPinterface.h"

int is_prop_prime(const mpz& a)
// returns 2 if a is definitely prime, 1 if a is probably prime and
// 0 if a is not prime
{
    return (mpz_probab_prime_p(a.get_mpz_t(), 40));
}

void printprimefactors(const vector<primepower>& v)
{
    for (primepower p : v) {
        cout << "(" << p.a << "^" << p.n << ")";
    }
    cout << endl;
}

vector<mpz> primefactor(mpz a)
// returns a vector of prime factors of a (using wheel factorization)
{
    vector<mpz> prime_divisors;
    for (mpz d : {2, 3, 5, 7}) {
        while (a % d == 0) {
            prime_divisors.push_back(d);
            a /= d;
        }
    }
    vector<mpz> increments = {2, 4, 2, 4, 6, 2, 6, 4, 2};
    int i = 0;
    for (mpz d = 11; d*d <= a; d += increments[i++]) {
        while (a % d == 0) {
            prime_divisors.push_back(d);
            a /= d;
        }
        if (i == 8)
            i = 0;
    }
    if (abs(a) > 1)
        prime_divisors.push_back(a);
    return prime_divisors;
}

mpz round(const mpf& a)
// rounds a to the nearest integer
{
    if ((a - floor(a)) < 0.5)
        return mpz(floor(a));
    else
        return mpz(ceil(a));
}

string decToBin(mpz a);

mpz pow(const mpz& a, const mpz& b)
// returns a^b (b must be a non-negative integer)
{
    if (a == 0 && b == 0) return 0;         // we define 0^0 = 0
    string bin = decToBin(b);               // get the power in binary
    mpz res = 1;
    for (int i = bin.size() - 1; i >= 0; --i) {
        res = res*res;
        if (bin[i] == '1') res = res*a;
    }
    return res;
}

vector<primepower> primepowerpairs(const vector<mpz>& v)
// returns a vector of (a, n)s where n is the power of a
{
    vector<mpz> temp = v;
    sort(temp.begin(), temp.end());

    vector<primepower> res;
    unsigned int n = 1;
    for (int i = 0; i < temp.size(); ++i) {
        if (temp[i] != temp[i + 1]) {
            primepower p = primepower(temp[i], n);
            res.push_back(p);
            n = 1;
        }
        else
            ++n;
    }
    return res;
}

mpz gcd(mpz a, mpz b)
// Euclidean algorithm for ordinary integers
{
    if (a % b == 0)
        return b;
    return gcd(b, a % b);
}

mpz Eulerphi(mpz n)
// Euler's phi-function
{
    mpz phi = 1;
    vector<primepower> factors = primepowerpairs(primefactor(n));
    for (primepower p : factors) {
        phi *= (pow(p.a, p.n) - pow(p.a, p.n - 1));
    }
    return phi;
}

string decToBin(mpz a)
// converts a positive integer in base 10 to binary (as a string in reverse!)
{
    string bin;
    while(a != 0) {
        if (a % 2 == 0) {
            a = a/2;
            bin += "0";
        }
        else {
            a = (a - 1)/2;
            bin += "1";
        }
    }
    return bin;
}

mpz modularExponentiation(mpz a, mpz b, mpz n)
// calculates a^b mod n (n must be a positive integer)
{
    string bin = decToBin(b);               // get the power in binary
    mpz res = 1;
    for (int i = bin.size() - 1; i >= 0; --i) {
        res = res*res % n;
        if (bin[i] == '1') res = res*a % n;
    }
    return res;
}

int JacobiSymbol(mpz a, mpz b)
// computes the Jacobi symbol (a/b), where b is an odd positive integer
{
    if (b % 2 == 0 || b < 1)
        throw runtime_error("The second argument has to be an odd positive integer");
    a %= b;
    if (a < 0)
        a += b;  // ensure the first argument is positive (no need to use the first supplementary law)

    int result = 1;
    while (a != 0) {
        while (a % 2 == 0) {            // use the second supplementary law
            a /= 2;
            if (b % 8 == 3 || b % 8 == 5)
                result = -result;
        }
        swap(a,b);
        if (a % 4 == 3 && b % 4 == 3)   // use quadratic reciprocity
            result = -result;
        a %= b;
    }

    if (b == 1)
        return result;
    else
        return 0;
}

mpz modularSqrt(mpz a, mpz p)
// calculates a square root of a modulo p if it exists (p is assumed to be a prime)
// the algorithm is due to Shanks
{
    if (p == 2) return a % p;
    mpz q = p - 1;
    int e = 0;
    while (q % 2 == 0) {
        ++e;
        q /= 2;
    }

    srand(time(0));
    mpz n = rand() % p;
    while (JacobiSymbol(n, p) == 1 || n % p == 0) {
        n = rand() % p;
    }
    mpz z = modularExponentiation(n, q, p);
    mpz y = z;
    mpz r = e;
    mpz t, m;
    mpz x = modularExponentiation(a, (q - 1)/2, p);
    mpz b = a*x*x % p;
    x = a*x % p;

    while (true) {
        if (b % p == 1 || b % p == 1 - p) return x % p;
        m = 1;
        while (modularExponentiation(b, pow(2, m), p) != 1) ++m;
        if (m == r) return 0;   // no square root of a modulo p
        t = modularExponentiation(y, pow(2, r - m - 1), p);
        y = t*t % p;
        r = m % p;
        x = x*t % p;
        b = b*y % p;
    }
}




