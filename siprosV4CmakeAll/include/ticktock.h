#pragma once

#include <chrono>
#include <iostream>
#include <fstream>
using namespace std;

#define TICK(x) auto bench_##x = chrono::steady_clock::now();
#define TICKAGAIN(x) bench_##x = chrono::steady_clock::now();
#define TOCK1ST(x)                                                                                                                              \
    {                                                                                                                                           \
        ofstream outputFile;                                                                                                                    \
        outputFile.open("timeUsed.txt", ios::out);                                                                                              \
        outputFile                                                                                                                              \
            << #x ": " << chrono::duration_cast<chrono::duration<double>>(chrono::steady_clock::now() - bench_##x).count() << "s" << std::endl; \
        outputFile.close();                                                                                                                     \
    }
#define TOCK(x)                                                                                                                                 \
    {                                                                                                                                           \
        ofstream outputFile;                                                                                                                    \
        outputFile.open("timeUsed.txt", ios::app);                                                                                              \
        outputFile                                                                                                                              \
            << #x ": " << chrono::duration_cast<chrono::duration<double>>(chrono::steady_clock::now() - bench_##x).count() << "s" << std::endl; \
        outputFile.close();                                                                                                                     \
    }
