#ifndef SHINGLES_HPP_INCLUDED
#define SHINGLES_HPP_INCLUDED

#include <vector>
#include <string>


class Shingles{
    public:
        typedef unsigned int Signature;
        Shingles();
        ~Shingles();
        //Signature compute(const std::vector<int>&) const;
        Signature compute(const std::vector<std::string>&) const;
    private:
        static const unsigned long bigPrimeNumber = 0x7FFFFFFF; // Same as 2^31 - 1
        unsigned int A;     // A, as in a X + b mod bigPrimeNumber
        unsigned int B;     // B, as in a X + b mod bigPrimeNumber
};

#endif  // SHINGLES_HPP_INCLUDED
