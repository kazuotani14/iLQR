#ifndef _LOCO_CAR_H_
#define _LOCO_CAR_H_

#include "standardIncludes.h"
#include "iLQR.h"

class LocoCar : public iLQR
{
  // Model parameters
  const double m;         // mass (kg)
  const double g ;
  const double L;         // wheelbase (m)
  const double b;         // CoG to rear axle
  const double a;         // CoG to front axle
  const double G_front;   // calculate load or specify front rear load directly
  const double G_rear;
  const double C_x;       // longitudinal stiffness
  const double C_alpha;   //lateral stiffness
  const double Iz;        // rotational inertia
  const double mu;        //5.2/G_rear
  const double mu_spin;   //4.3/G_rear

  //Cost weights
  Vec2d cu;             //control cost
  Vec2d cdu;            //change in control cost
  VecXd cx;             //running cost for velocities
  VecXd px;             //smoothness scales for running cost
  double c_drift;       //small reward for drifing
  double kp_obs;        //obstacle avoidance costs
  double kv_obs;
  double dist_thres;

  Vec2d tire_dyn(double Ux, double Ux_cmd, double mu, double mu_slide,
                    double Fz, double C_x, double C_alpha, double alpha);

  VecXd dynamics(const VecXd &x, const VecXd &u);
  double cost(const VecXd &x, const VecXd &u);
  double final_cost(const VecXd &x);
  VecXd integrate_dynamics(const VecXd &x, const VecXd u);

public:
  Vec2d obs; //position of obstacle in map frame. set obs to (0,0)

  LocoCar(): m(2.35), g(9.81), L(0.257), a(0.11372), b(0.14328),
             C_x(65), C_alpha(55), Iz(0.025), mu(0.45), mu_spin(0.2),
             G_front(12.852550506), G_rear(10.200949494),
             c_drift(-0.001), kp_obs(0.5), kv_obs(0.1), dist_thres(0.5)
 {
   cu << 1e-3, 1e-3;
   cdu << 5, 5e-1;
   cx.resize(6);
   px.resize(6);
   cx << 0.5, 0.1, 0.4, 0.05, 0.005, 0.002;
   px << 0.01, 0.01, 0.1, 0.01, 0.01, 0.1;
 }
};

#endif
