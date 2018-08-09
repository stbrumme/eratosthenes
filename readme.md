This is a mirror of my library hosted at https://create.stephan-brumme.com/eratosthenes/

# Parallel Sieve of Eratosthenes

There are 50,847,534 prime numbers between 2 and 1,000,000,000.
For reasons unknown, I wanted to know whether a C++ implementation can find all these numbers on my modest desktop computer (Intel Core i7 860, quad-core, hyperthreading, 2.8GHz) in less than 1 second with a simple algorithm such as the [Sieve of Eratosthenes](https://en.wikipedia.org/wiki/Sieve_of_Eratosthenes).

By only slightly modifying the basic algorithm and unleashing the raw power of [OpenMP](https://www.openmp.org/), I achieved a 20x speed-up.

## More ...
See my website https://create.stephan-brumme.com/eratosthenes/ for code examples and benchmarks.
You'll find the Sieve in PHP and Javascript/Asm.js, too (but not parallelized, though).
