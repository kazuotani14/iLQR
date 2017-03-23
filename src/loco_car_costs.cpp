#include "loco_car.h"

// TODO vectorize so we can do parallel line search. for now, do sequential

double LocoCar::cost(const VecXd &x, const VecXd &u)
{

  // Input: n=8 state vector(s). 6 states, 2 precalculated du
  // columns of x and u will be each hypothesis

  // TODO make these member variables
  //Coefficients/weights
  Vec2d cu(1e-3, 1e-3); //control cost
  Vec2d cdu(1e-1, 1e-2); cdu *= 50;  //change in control cost

  Eigen::Vector3d cx(0.5, 0.1, 0.4); //running cost for pose
  Eigen::Vector3d cdx(0.05, 0.005, 0.002); // running cost for velocities
  Eigen::Vector3d px(0.01, 0.01, 0.1); //smoothness scales for running cost
  double c_drift = -0.001;

  // obstacle costs
  double kp_obs = 0.5;
  double kv_obs = 0.1;
  double dist_thres = 0.5;

  //control cost
  double lu = cu.dot(elem_square(u-x.segment<2>(3)));
  // double ldu = cdu.dot(elem_square(x.tail(2)));
  double ldu = 0; //TODO implement this

  //running cost
  Eigen::Vector3d dist = x.head(3) - x_d.head(3);
  Eigen::Vector3d vdist = x.segment<3>(3) - x_d.segment<3>(3); //4:6 in matlab
  double lx = cx.dot(sabs(dist, px)) + cdx.dot(sabs(vdist, px));

  //small reward for drifing
  double ld = c_drift * (sabs(x(5),1) - 0.2);

  //obstacle avoidance cost
  double lobs = 0;
  if (obs.sum()>1e-7)
  {
    Vec2d pos = x.segment<2>(0);
    Vec2d vel = x.segment<2>(3);
    Vec2d vec2obs = obs - pos;
    double dist2obs = vec2obs.norm();
    double velnorm = vel.norm();

    double Ustatic = sqr((1/dist2obs - 1/dist_thres));
    if (Ustatic >= dist_thres)
      Ustatic = 0;

    double dotprod = vec2obs.dot(vel); // vec . vel = |vec||vel|cos(theta)
    double Udynamic = dotprod / (dist2obs*velnorm);
    if (Udynamic < 0)
      Udynamic = 0;

    lobs = kp_obs*Ustatic + kv_obs*Udynamic;
  }

  // std::cout << "costs: " << lu << ' ' << lx << ' ' << ldu << ' ' << ld << ' ' << lobs <<'\n';
  double total_cost = lu + lx + ldu + ld + lobs;
  return total_cost;
}; //cost

double LocoCar::final_cost(const VecXd &x)
{
  // Weights
  VecXd cf(6); cf << 10, 10, 1, 0.1, 0.1, 0.1; // final cost
  VecXd pf(6); pf << 0.01, 0.01, 0.1, 0.1, 0.1, 0.1; //smoothness scales

  VecXd dist = x.segment<6>(0) - x_d;

  double final_cost = cf.dot((sabs(dist,pf)) + elem_square(sabs(dist,pf)));

  return final_cost;
};  //final_cost
