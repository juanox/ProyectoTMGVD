#include "Shingles.hpp"
#include <cstddef>      // NULL, std::size_t
#include <cstdlib>      // std::rand, std::srand
#include <ctime>        // std::time


Shingles::Shingles() {
    std::srand(std::time(NULL));
    A = (std::rand() % bigPrimeNumber) + 1;
    B = (std::rand() % bigPrimeNumber) + 1;
}

Shingles::~Shingles(){}

Shingles::Signature
//Shingles::compute(const std::vector<int>& seqs) const {
Shingles::compute(const std::vector<std::string>& seqs) const {

    std::hash<std::string> stringHash;
    Signature minShingleHash = bigPrimeNumber;

    for (auto& val : seqs) {
    //for (unsigned int i=0; i<seqs.size(); i++) {
	//std::string val = seqs[i];
        //std::size_t hashedWord = stringHash(std::to_string(val));
        std::size_t hashedWord = stringHash(val);
        std::size_t shingleID = hashedWord;
        Signature shingleHash = (((unsigned long) A * (unsigned long) shingleID) + B) % bigPrimeNumber;

        if (shingleHash < minShingleHash) {
            minShingleHash = shingleHash;
        }
    }

    return minShingleHash;
}
