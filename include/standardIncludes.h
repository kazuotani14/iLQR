#pragma once
#include <vector>
#include <eigen/Eigen/Core>
#include <eigen/Eigen/Eigenvalues>
#include <iostream>
#ifndef DEBUG
#define ASSERT(x) {}
#else
#define ASSERT(x) assert(x)
#endif

using namespace std;
using namespace Eigen;
template<typename T>
inline T sqr(const T &val){ return val*val; }
const double pi = M_PI;
const double timeDelta = 1.0/60.0;
