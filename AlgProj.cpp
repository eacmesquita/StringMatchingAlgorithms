#include <iostream>
#include<vector>
#include<string>
#include<algorithm>
#include <chrono>
#include <random>
#include <fstream>

using std::chrono::high_resolution_clock;
using std::chrono::duration_cast;
using std::chrono::duration;
using std::chrono::milliseconds;

std::random_device rd;
std::mt19937 mt(rd());
std::uniform_real_distribution<double> dist(0, 1000);

/*
This function will return a randomly generated string, with the alphabet being restricted to A,T,G and C.
The oder rand() and srand() pseudorandom generators do not provide much degree of randomness, and results get repeated.
So, the C++ 11 random library is used here for better randomness and degree of change in the generated strings
This is by far the most time-consuming part of the program.
*/
void GenerateRandomString(char* pArray, size_t size)
{
    static const char alphabet[] = { 'A','T','G','C' };
    static const int alphabetSize = sizeof(alphabet) / sizeof(alphabet[0]);
    for (size_t i = 0; i < size; ++i)
    {
        pArray[i] = alphabet[(int)dist(mt)% alphabetSize];
    }
}

//Algorithms section

//KMP search - compute the prefix
size_t* computePrefix(char* pattern, size_t m)
{
    size_t* prefix = new size_t[m];
    prefix[0] = 0;
    size_t k = 0;
    for (size_t i = 1; i < m; i++)
    {
        while (k > 0 && pattern[k] != pattern[i])
        {
            k = prefix[k - 1];
        }
        if (pattern[k] == pattern[i])
        {
            k = k + 1;
        }
        prefix[i] = k;
    }
    return prefix;
}

//KMP search invoked from main
std::vector<size_t> KMPStringSearch(char* text, size_t n, char* pattern, size_t m)
{
    std::vector<size_t> occurrences;
    size_t* prefix = computePrefix(pattern, m);
    size_t matched_pos = 0;
    for (size_t i = 0; i < n; i++)
    {
        while (matched_pos > 0 && pattern[matched_pos] != text[i])
            matched_pos = prefix[matched_pos - 1];

        if (pattern[matched_pos] == text[i])
            matched_pos = matched_pos + 1;

        if (matched_pos == m)
        {
            occurrences.push_back(i - (m - 1));
            matched_pos = prefix[matched_pos - 1];
        }
    }
    delete[] prefix;
    return occurrences;
}

//Rabin Karp search implementation
std::vector<size_t> RabinKarp(char* text, size_t n, char* pattern, size_t m)
{
    std::vector<size_t> occurrences;
    size_t i, j;
    int p = 0;
    int t = 0;
    int h = 1;
    const int d = 256;
    int q = 101;

    for (i = 0; i < m - 1; i++)
        h = (h * d) % q;

    for (i = 0; i < m; i++)
    {
        p = (d * p + pattern[i]) % q;
        t = (d * t + text[i]) % q;
    }
    for (i = 0; i <= n - m; i++)
    {
        if (p == t)
        {
            for (j = 0; j < m; j++)
            {
                if (text[i + j] != pattern[j])
                    break;
            }
            if (j == m)
            {
                occurrences.push_back(i);
            }
        }
        if (i < n - m)
        {
            t = (d * (t - text[i] * h) + text[i + m]) % q;
            if (t < 0)
                t = (t + q);
        }
    }
    return occurrences;
}

//Algorithms section ends here

void main()
{
    
    std::vector < std::string> vecResults;
    const int ITERATIONS_PER_INPUT_SIZE = 10;
    int NUM_ROUNDS = 10; //how many of the iterations will be executed

    static const size_t nTenTo2 = 100;
    static const size_t nTenTo3 = 10 * nTenTo2;
    static const size_t nTenTo4 = 10 * nTenTo3;
    static const size_t nTenTo5 = 10 * nTenTo4;
    static const size_t nTenTo6 = 10 * nTenTo5;
    static const size_t nTenTo7 = 10 * nTenTo6;
    static const size_t nTenTo8 = 10 * nTenTo7;
    static const size_t nTenTo9 = 10 * nTenTo8;

    char* pattern; char* text;

    size_t mValues[] = {8, 8, 8, 12, 12, 12, 16, 16, 20, 20};
    size_t nValues[] = { nTenTo7, nTenTo8, nTenTo9, nTenTo7, nTenTo8, nTenTo9, nTenTo8, nTenTo9, nTenTo8, nTenTo9 };

    char demo;
    std::cout << "Demo run?(Y/N): ";
    std::cin >> demo;

    if (demo == 'Y') {
        NUM_ROUNDS = 1; //only execute first iteration
    }

    for(int i=0; i < NUM_ROUNDS; i++) {
        vecResults.push_back("m = " + std::to_string(mValues[i]) + " and n = " + std::to_string(nValues[i]));
        double kmpTotalTime = 0, rkTotalTime = 0;
        double kmpAverageTime, rkAverageTime;
        for (int j = 0; j < ITERATIONS_PER_INPUT_SIZE; j++)
        {
            pattern = new char[mValues[i]];
            text = new char[nValues[i]];
            vecResults.push_back("\nIteration " + std::to_string(j));

            std::cout << "\nGenerating pattern of length " + std::to_string(mValues[i]);
            GenerateRandomString(pattern, mValues[i]);
            std::cout << "\nGenerating text of length " + std::to_string(nValues[i]);
            GenerateRandomString(text, nValues[i]);

            std::cout << "\nNow executing KMP search";
            auto t1 = high_resolution_clock::now();
            std::vector<size_t> kmpOccurrences = KMPStringSearch(text, nValues[i], pattern, mValues[i]);
            auto t2 = high_resolution_clock::now();

            std::cout << "\nNow executing Rabin Karp search";
            auto t3 = high_resolution_clock::now();
            std::vector<size_t> RabinKarpOccurrences = RabinKarp(text, nValues[i], pattern, mValues[i]);
            auto t4 = high_resolution_clock::now();

            delete[] pattern;
            delete[] text;

            duration<double, std::milli> kmpTime = t2 - t1;
            duration<double, std::milli> rabinKarpTime = t4 - t3;

            kmpTotalTime += kmpTime.count();
            rkTotalTime += rabinKarpTime.count();

            vecResults.push_back("KMP running time: " + std::to_string(kmpTime.count()) + 
                " ms and Rabin Karp running time: " + std::to_string(rabinKarpTime.count()) + " ms");

            //issue in the implementation
            if (kmpOccurrences.size() != RabinKarpOccurrences.size()) {
                throw "Error - KMP and Rabin Karp did not match in number of occurrences, algorith implementation is incorrect";
            }
            else {
                vecResults.push_back("Number of matches in this iteration - " + std::to_string(kmpOccurrences.size()));
            }          
        }
        kmpAverageTime = kmpTotalTime / (double)ITERATIONS_PER_INPUT_SIZE;
        rkAverageTime = rkTotalTime / (double)ITERATIONS_PER_INPUT_SIZE;

        vecResults.push_back("\nAverage KMP running time r1: " + std::to_string(kmpTotalTime/ITERATIONS_PER_INPUT_SIZE ) 
            + " ms and average Rabin Karp running time r2: " + std::to_string(rkTotalTime/ITERATIONS_PER_INPUT_SIZE) + " ms\n\n");
        }

        std::cout << "\nWriting results to log file";
        std::fstream outputFile;
        outputFile.open("AlgorithmProjectOutput.txt", std::ios::out);
        for (std::string line: vecResults) {
            outputFile << line << "\n";
        }
}
