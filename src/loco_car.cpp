#include "loco_car.h"

Eigen::Vector2d LocoCar::tire_dyn(double Ux, double Ux_cmd, double mu, double mu_slide,
                  double Fz, double C_x, double C_alpha, double alpha)
{
  double Fx, Fy;
  Eigen::Vector2d F_tire;

  //longitude wheel slip
  double K;
  if (Ux_cmd == Ux){
    K = 0;
  }
  else if (Ux == 0){
    Fx = sgn(Ux_cmd)*mu*Fz;
    Fy = 0;
    F_tire << Fx, Fy;
    return F_tire;
  }
  else{
    K = (Ux_cmd-Ux)/std::abs(Ux);
  }

  //instead of avoiding -1, now look for positive equivalent
  int reverse = 1;
  if (K < 0){
    reverse = -1;
    K = std::abs(K);
  }

  // alpha > pi/2 cannot be adapted to this formula
  // because of the use of tan(). Use the equivalent angle instead.
  if (std::abs(alpha) > pi/2){
    alpha = (std::abs(alpha)-pi/2)*sgn(alpha);
  }

  double gamm = sqrt(sqr(C_x)*sqr(K/(1+K))+ sqr(C_alpha)*sqr(tan(alpha)/(1+K)));

  double F;
  if (gamm <= 3*mu*Fz){
    F = gamm - 1/(3*mu*Fz)*sqr(gamm) + 1/(27*sqr(mu)*sqr(Fz))*(sqr(gamm)*gamm);
  }
  else{
    // more accurate modeling with peak friction value
    F = (mu_slide + (mu-mu_slide)/(1 + sqr((gamm-3*mu*Fz)/27)))*Fz;
  }

  if (gamm == 0){
    Fx = 0;
    Fy = 0;
  }
  else{
    Fx = C_x/gamm * (K/(1+K)) * F * reverse;
    Fy = -C_alpha/gamm * (tan(alpha)/(1+K)) * F;
  }

  F_tire << Fx, Fy;
  return F_tire;
}

Eigen::VectorXd LocoCar::dynamics(const Eigen::VectorXd &x, const Eigen::Vector2d &u)
{
  double pos_x = x(0); double pos_y = x(1); double pos_phi = x(2);
  double Ux = x(3); double Uy = x(4); double r = x(5);

  double Ux_cmd = u(0); double delta = u(1);

  // Lateral slip angles
  double alpha_F, alpha_R;
  if (Ux==0 && Uy==0){ // no slip
    alpha_F = 0;
    alpha_R = 0;
  }
  else if (Ux==0){ // perfect side slide
    alpha_F = pi/2*sgn(Uy) - delta;
    alpha_R = pi/2*sgn(Uy);
  }
  else if (Ux < 0){ // rare ken block situations
    alpha_F = (sgn(Uy)*pi)-atan((Uy+a*r)/std::abs(Ux))-delta;
    alpha_R = (sgn(Uy)*pi)-atan((Uy-b*r)/std::abs(Ux));
  }
  else{ //normal situation
    alpha_F = atan((Uy+a*r)/std::abs(Ux))-delta;
    alpha_R = atan((Uy-b*r)/std::abs(Ux));
  }

  // safety that keeps alphas in valid range
  alpha_F = wrap_to_pi(alpha_F);
  alpha_R = wrap_to_pi(alpha_R);
  // std::cout << alpha_F << ' ' << alpha_R << '\n';

  Eigen::Vector2d Ff = tire_dyn(Ux, Ux, mu, mu_spin, G_front, C_x, C_alpha, alpha_F);
  Eigen::Vector2d Fr = tire_dyn(Ux, Ux_cmd, mu, mu_spin, G_rear, C_x, C_alpha, alpha_R);
  double Fxf, Fyf, Fxr, Fyr;
  Fxf = Ff(0); Fyf = Ff(1); Fxr = Fr(0); Fyr = Fr(1);

  // Vehicle body dynamics
  // std::cout << Fxf << ' ' << Fyf << ' ' << Fxr << ' ' << Fyr << '\n';
  double r_dot = (a*Fyf*cos(delta)-b*Fyr)/Iz;
  double Ux_dot = (Fxr-Fyf*sin(delta))/m+r*Uy;
  double Uy_dot = (Fyf*cos(delta)+Fyr)/m-r*Ux;

  // translate dx to terrain frame
  double U = sqrt(sqr(Ux)+sqr(Uy));
  double beta;
  if (Ux == 0 && Uy == 0){
    beta = 0;
  }
  else if (Ux == 0){
    beta = pi/2*sgn(Uy);
  }
  else if (Ux < 0 && Uy == 0){
      beta = pi;
  }
  else if (Ux < 0){
    beta = sgn(Uy)*pi-atan(Uy/std::abs(Ux));
  }
  else{
    beta = atan(Uy/std::abs(Ux));
  }
  beta = wrap_to_pi(beta);

  double Ux_terrain = U*cos(beta+pos_phi);
  double Uy_terrain = U*sin(beta+pos_phi);

  Eigen::VectorXd dx(6);
  dx << Ux_terrain, Uy_terrain, r, Ux_dot, Uy_dot, r_dot;
  return dx;
}

double LocoCar::cost(const Eigen::VectorXd &x, const Eigen::Vector2d &u){
  double cost;

  return cost;
}; // cost function

double LocoCar::final_cost(const Eigen::VectorXd &x){
  double final_cost;

  return final_cost;
};  // final cost

int main() {
  Eigen::VectorXd x(6);
  x << 1,2,3,1,2,3;
  Eigen::Vector2d u(1,1);

  LocoCar loco;
  Eigen::VectorXd dx(6);
  dx = loco.dynamics(x,u);
  print_vec(dx);
  return 0;
}
