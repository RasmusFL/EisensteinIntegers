## EisensteinIntegers

This repository is for all the code, I have written concerning the Eisenstein integers. This manual is practically identical to the one about the GaussianIntegers repository. The code contains algorithms for the following tasks:

* Simple arithmetic including +, -, * and Euclidean division with the Eisenstein integers
* Taking conjugates, determining associates, computing powers, modular exponentiation
* The Euclidean algorithm and the extended Euclidean algorithm
* Solving modular linear equations and the Chinese remainder theorem
* Prime factorization
* Computing the cubic residue symbol

# Prerequisites and setup

Compiling the code requires GMP (The GNU Multiple Precision Arithmetic Library) installed properly. For information on download and setup, consult their [website](https://gmplib.org/).

The files "GMPinterface.h" and "GMPinterface.cpp" handle all the interactions with the main files "eisenstein_integers.h" and "eisenstein_integers.cpp" with the GMP library. Note that the multiple precision integer class `mpz_class` is simply abbreviated to `mpz` and similarly for `mpf_class`. The GMP interface files also contain some basic number theoretic algorithms including the Jacobi symbol, a simple factorization algorithm, Euler's phi function and an algorithm for computing modular square roots (modulo a prime). 

In order to display an Eisenstein integer properly to the console, you must enable ANSI Greek in your `main.cpp` file by having the line SetConsoleOutputCP(1253); in the beginning of the file (requires the windows library). See the `main.cpp` file in the repository for a working example.

# The `EI` class

The main object is the `EI` class which has three constructors. `EI()` is just 0, `EI(mpz a)` is just an Eisenstein integer without any coefficient on ω and `EI(mpz a, mpz b)` is a + bω. Ordinary arithmetic and writing to the console is supported. For example,

```c++
EI alpha {2, 3};
EI beta {-5, 8};
cout << alpha + beta << endl;
cout << alpha - beta << endl;
cout << alpha * beta << endl;
```
will output
```c++
-3+11ω
7-5ω
-34-23ω
```

The division operator / is also supported. This applies Euclidean division and outputs the quotient, so the output of
```c++
cout << beta / alpha;
```
becomes
```c++
4+4ω
```
because 5 - 8ω = (2 + 3ω)(4 + 4ω) - 1. The modulo operator % also works as intended. For example
```c++
cout << beta % alpha; // outputs -1
```

# Some examples of usage

Here are further examples of using the library. What is the greatest common divisor of, say, 36 - 56ω and -144 + 913ω?

```c++
cout << EIgcd({36, -56}, {-144, 913}); // outputs -4-ω
```
The extended Euclidean algorithm takes two Eisenstein integers and outputs a vector with the greatest common divisor and the two coefficients in Bézout's lemma. For example

```c++
vector<EI> v = EIexgcd({36, -56}, {-144, 913});
for (EI a : v)
  cout << a << endl;
```
outputs:

```c++
-4-ω
116+139ω
7+12ω
```
and we can verify that (116 + 139ω)(36 - 56ω) + (7 + 12ω)(-144 + 913ω) = -4 - ω. Let us now solve a modular linear equation. For example, let us solve (36 - 56ω)x = -144 + 913ω mod 587:

```c++
cout << EImodlinearsolve({36, -56}, {-144, 913}, {587}); // outputs -152+179ω
```
And we can verify that (36 - 56ω)(-152 + 179ω) = -144 + 913ω modulo 587. On the other hand, there is no solution to (36 - 56ω)(-152 + 179ω) = -144 + 913ω modulo 8 + 2ω since

```c++
EImodlinearsolve({36, -56}, {-144, 913}, {8,2}); // outputs "No solution to equation: (36ω-56ω)rho = -144+913ω (mod 8+2ω)"
```
The function `EIprimefactor` factors an Eisenstein integer into a product og primes and a single unit (1, -1, ω, -ω, ω^2 or -ω^2). The output is a vector of these Eisenstein integers. To write the factorization in a nice way as a product, use the function `printEIproduct`. As an example,

```c++
printEIproduct(EIprimefactor({894582, -1870956}));
```
will output

```c++
(-1-ω)(1-ω)^5(2)(5+6ω)(95-18ω)(-115+32ω)
```
as the factorization of 894582 - 1870956ω. As a final example, we may compute the cubic residue symbol (alpha/beta)_3, where beta is not divisible by 1 - ω. For example

```c++
cout << CubicResSymbol({36, -26}, {-144, 913}); // outputs -1-ω
```
This means that 36 - 26ω is not a cubic residue modulo -144 + 913ω. 
