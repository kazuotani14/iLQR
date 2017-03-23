#ifndef _LOCO_CAR_H_
#define _LOCO_CAR_H_

#include "standardIncludes.h"
#include "iLQR.h"

class LocoCar : public iLQR
{
  // Model parameters
  const double m;         // mass (kg)
  const double g ;
  const double L;        // wheelbase (m)
  const double b;      // CoG to rear axle
  const double a;          // CoG to front axle
  const double G_front;   // calculated load or specify front rear load directly
  const double G_rear;
  const double C_x;       // longitudinal stiffness
  const double C_alpha;   //lateral stiffness
  const double Iz;     // rotational inertia
  const double mu;      //5.2/G_rear
  const double mu_spin;  //4.3/G_rear

  Vec2d tire_dyn(double Ux, double Ux_cmd, double mu, double mu_slide,
                    double Fz, double C_x, double C_alpha, double alpha);

  // void finite_difference( )

public:

  Vec2d obs; //position of obstacle in map frame. set obs(0) to 9999 when it doesn't exist

  LocoCar(): m(2.35), g(9.81), L(0.257), b(0.14328), C_x(50), C_alpha(45),
             Iz(0.045), mu(0.75), mu_spin(0.2), a(0.11372), G_front(12.852550506),
             G_rear(10.200949494) { }

  VecXd dynamics(const VecXd &x, const VecXd &u) ; //dynamics

  double cost(const VecXd &x, const VecXd &u) ;

  double final_cost(const VecXd &x) ;

};

#endif
