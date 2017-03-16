#include "standardIncludes.h"
#include "loco_car.h"
#include <fstream>

using namespace std;
using namespace Eigen;

int main()
{
  int K = 50; // length of trajectory
  LocoCar loco;

  Eigen::VectorXd x0(6);
  x0 << 0, 0, 0, 2, 0, 0;
  Eigen::VectorXd u0(2);
  u0 << 1, 0.5;
  Eigen::VectorXd x_d(6);
  x_d << 3, 3, pi/2, 0, 0, 0;
  Eigen::VectorXd obs(2);
  obs << 1.5, 0;

  loco.obs = obs;

  loco.generateFeedbackController(x0, x_d, K, 2, u0); // generates trajectory xs

  std::cout<< "saving stuff..\n";
  ofstream ofs("locotrajectory.txt", ios::out);
  for (int k = 0; k<K-1; k++)
  {
    ofs << loco.us[k][0] << ", " << loco.us[k][1] << endl;
  }
  return 1;
}
