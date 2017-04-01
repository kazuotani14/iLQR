#ifndef _STANDARD_INCLUDES_H_
#define _STANDARD_INCLUDES_H_

#include <vector>
#include "eigen/Eigen/Core"
#include "eigen/Eigen/Eigenvalues"
#include "eigen/Eigen/StdVector"
#include <iostream>
#include <cstdio>
#include <math.h>
#include <time.h>
#include <ctime>
#define EIGEN_USE_NEW_STDVECTOR

#define Eye2 Eigen::Matrix2d::Identity(2,2);

using Eigen::Vector2d;
using Eigen::VectorXd;
using Eigen::MatrixXd;

typedef Eigen::Vector2d Vec2d;
typedef Eigen::VectorXd VecXd;
typedef Eigen::MatrixXd MatXd;
typedef std::vector<VecXd, Eigen::aligned_allocator<VecXd> > VecOfVecXd;
typedef std::vector<MatXd, Eigen::aligned_allocator<MatXd> >  VecOfMatXd;

//---------------------------------
// Math helper functions
//---------------------------------
const double pi = M_PI;

template<typename T>
inline T sqr(const T &val){ return val*val; }

template<typename T>
inline T cube(const T &val){ return val*val*val; }

template <typename T>
inline int sgn(T &val) {return (T(0) < val) - (val < T(0)); }

inline double sabs(double x, double y)
{
  //Differentiable "soft" absolute value function
  return sqrt(sqr(x)+sqr(y))-y;
}

//---------------------------------
// Floating-point modulo
// The result (the remainder) has same sign as the divisor.
// Similar to matlab's mod(); Not similar to fmod() -   Mod(-3,4)= 1   fmod(-3,4)= -3
// From http://stackoverflow.com/questions/4633177/c-how-to-wrap-a-float-to-the-interval-pi-pi
template<typename T>
inline T Mod(T x, T y)
{
    //static_assert(!std::numeric_limits<T>::is_exact , "Mod: floating-point type expected");

    if (0. == y)
        return x;
    double m= x - y * floor(x/y);
    // handle boundary cases resulted from floating-point cut off:
    if (y > 0)              // modulo range: [0..y)
    {
        if (m>=y)           // Mod(-1e-16             , 360.    ): m= 360.
            return 0;

        if (m<0 )
        {
            if (y+m == y)
                return 0  ; // just in case...
            else
                return y+m; // Mod(106.81415022205296 , _TWO_PI ): m= -1.421e-14
        }
    }
    else                    // modulo range: (y..0]
    {
        if (m<=y)           // Mod(1e-16              , -360.   ): m= -360.
            return 0;

        if (m>0 )
        {
            if (y+m == y)
                return 0  ; // just in case...
            else
                return y+m; // Mod(-106.81415022205296, -_TWO_PI): m= 1.421e-14
        }
    }

    return m;
}

inline double wrap_to_pi(double angle)
{
  return Mod(angle+pi, 2*pi) - pi;
}

//---------------------------------
// Eigen-specific helper functions
//---------------------------------

inline VectorXd elem_square(const VectorXd &vec)
{
  return vec.array().square().matrix();
}

inline VectorXd elem_sqrt(const VectorXd &vec)
{
  return vec.array().sqrt().matrix();
}

inline VectorXd sabs(const VectorXd &vec, const VectorXd &p)
{
  //Differentiable "soft" absolute value function
  VecXd sum = elem_sqrt(elem_square(vec)+elem_square(p));
  return sum - p;
}

#endif
