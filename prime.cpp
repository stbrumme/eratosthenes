// //////////////////////////////////////////////////////////
// Copyright (c) 2011 Stephan Brumme. All rights reserved.
// see http://create.stephan-brumme.com/disclaimer.html
//

#ifdef _WIN32
#define USE_WINDOWS_TIMER
#endif

typedef __int32 Number;  // or for the brave ones: typedef __int64 Number;

// number of primes between 2 and x:
//             10 =>           4
//            100 =>          25
//          1,000 =>         168
//         10,000 =>       1,229
//        100,000 =>       9,592
//      1,000,000 =>      78,498
//     10,000,000 =>     664,579
//    100,000,000 =>   5,761,455
//  1,000,000,000 =>  50,847,534
// 10,000,000,000 => 455,052,511
const Number lastNumber = 1*1000*1000*1000LL;

#include <stdio.h>
#include <math.h>

// OpenMP
#include <omp.h>
// block size, 128k is the sweet spot for Core i7 cpus
const int sliceSize = 128*1024;


// simple serial sieve of Eratosthenes
int eratosthenes(Number lastNumber)
{
  // initialize
  char* isPrime = new char[lastNumber+1];
  for (Number i = 0; i <= lastNumber; i++)
    isPrime[i] = 1;

  // find all non-primes
  for (Number i = 2; i*i <= lastNumber; i++)
    if (isPrime[i])
      for (Number j = i*i; j <= lastNumber; j += i)
        isPrime[j] = 0;

  // sieve is complete, count primes
  int found = 0;
  for (Number i = 2; i <= lastNumber; i++)
    found += isPrime[i];

  delete[] isPrime;
  return found;
}


// odd-only sieve
int eratosthenesOdd(Number lastNumber, bool useOpenMP)
{
  // enable/disable OpenMP
  omp_set_num_threads(useOpenMP ? omp_get_num_procs() : 1);

  // instead of i*i <= lastNumber we write i <= lastNumberSquareRoot to help OpenMP
  const Number lastNumberSqrt = (int)sqrt((double)lastNumber);

  Number memorySize = (lastNumber-1)/2;

  // initialize
  char* isPrime = new char[memorySize+1];
#pragma omp parallel for
  for (Number i = 0; i <= memorySize; i++)
    isPrime[i] = 1;

  // find all odd non-primes
#pragma omp parallel for schedule(dynamic)
  for (Number i = 3; i <= lastNumberSqrt; i += 2)
    if (isPrime[i/2])
      for (Number j = i*i; j <= lastNumber; j += 2*i)
        isPrime[j/2] = 0;

  // sieve is complete, count primes
  int found = lastNumber >= 2 ? 1 : 0;
#pragma omp parallel for reduction(+:found)
  for (Number i = 1; i <= memorySize; i++)
    found += isPrime[i];

  delete[] isPrime;
  return found;
}


// process only odd numbers of a specified block
int eratosthenesOddSingleBlock(const Number from, const Number to)
{
  const Number memorySize = (to - from + 1) / 2;

  // initialize
  char* isPrime = new char[memorySize];
  for (Number i = 0; i < memorySize; i++)
    isPrime[i] = 1;

  for (Number i = 3; i*i <= to; i += 2)
  {
    // skip multiples of three: 9, 15, 21, 27, ...
    if (i >= 3*3 && i % 3 == 0)
      continue;
    // skip multiples of five
    if (i >= 5*5 && i % 5 == 0)
      continue;
    // skip multiples of seven
    if (i >= 7*7 && i % 7 == 0)
      continue;
    // skip multiples of eleven
    if (i >= 11*11 && i % 11 == 0)
      continue;
    // skip multiples of thirteen
    if (i >= 13*13 && i % 13 == 0)
      continue;

    // skip numbers before current slice
    Number minJ = ((from+i-1)/i)*i;
    if (minJ < i*i)
      minJ = i*i;

    // start value must be odd
    if ((minJ & 1) == 0)
      minJ += i;

    // find all odd non-primes
    for (Number j = minJ; j <= to; j += 2*i)
    {
      Number index = j - from;
      isPrime[index/2] = 0;
    }
  }

  // count primes in this block
  int found = 0;
  for (Number i = 0; i < memorySize; i++)
    found += isPrime[i];
  // 2 is not odd => include on demand
  if (from <= 2)
    found++;

  delete[] isPrime;
  return found;
}


// process slice-by-slice, odd numbers only
int eratosthenesBlockwise(Number lastNumber, Number sliceSize, bool useOpenMP)
{
  // enable/disable OpenMP
  omp_set_num_threads(useOpenMP ? omp_get_num_procs() : 1);

  int found = 0;
  // each slices covers ["from" ... "to"], incl. "from" and "to"
#pragma omp parallel for reduction(+:found)
  for (Number from = 2; from <= lastNumber; from += sliceSize)
  {
    Number to = from + sliceSize;
    if (to > lastNumber)
      to = lastNumber;

    found += eratosthenesOddSingleBlock(from, to);
  }

  return found;
}


// timing
#ifdef USE_WINDOWS_TIMER
#include <windows.h>
#else
#include <sys/time.h>
#endif

double seconds()
{
#ifdef USE_WINDOWS_TIMER
  LARGE_INTEGER frequency, now;
  QueryPerformanceFrequency(&frequency);
  QueryPerformanceCounter  (&now);
  return now.QuadPart / double(frequency.QuadPart);
#else
  timeval now;
  gettimeofday(&now, NULL);
  return now.tv_sec + now.tv_usec/1000000.0;
#endif
}



int __cdecl main(Number argc, char* argv[])
{
  printf("Primes between 2 and %d\n\n", lastNumber);

  printf("OpenMP uses up to %d threads running on %d processors\n\n",
         omp_get_max_threads(), omp_get_num_procs());

  printf("Simple Sieve\n");
  double startTime = seconds();
  int found = eratosthenes(lastNumber);
  double duration  = seconds() - startTime;
  printf("%d primes found in %.3fs\n\n", found, duration);

  printf("Only odd numbers\n");
  startTime = seconds();
  found = eratosthenesOdd(lastNumber, false);
  duration  = seconds() - startTime;
  printf("%d primes found in %.3fs\n\n", found, duration);

  printf("Only odd numbers, OpenMP\n");
  startTime = seconds();
  found = eratosthenesOdd(lastNumber, true);
  duration  = seconds() - startTime;
  printf("%d primes found in %.3fs\n\n", found, duration);

  printf("Blockwise, only odd numbers\n");
  startTime = seconds();
  found = eratosthenesBlockwise(lastNumber, 2*sliceSize, false);
  duration  = seconds() - startTime;
  printf("%d primes found in %.3fs\n\n", found, duration);

  printf("Blockwise, only odd numbers, OpenMP\n");
  startTime = seconds();
  found = eratosthenesBlockwise(lastNumber, 2*sliceSize, true);
  duration  = seconds() - startTime;
  printf("%d primes found in %.3fs\n\n", found, duration);

  return 0;
}
