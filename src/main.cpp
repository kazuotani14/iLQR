#include "standardIncludes.h"
#include "acrobot.h"
#include <fstream>

using namespace std;
using namespace Eigen;

int main()
{
  int K = 500; // length of trajectory
  Acrobot acrobot;
  acrobot.init();
  Vector4d startState(0,0,0,0); // angle 0,1, angular vel 0,1
  Vector4d endState(pi,0,0,0);  // raised vertically above the pivot
  acrobot.generateFeedbackController(startState, endState, K, 1); // generates trajectory xs

  ofstream ofs("acrobottrajectory.txt", ios::out);
  for (int k = 0; k<K; k++)
  {
    ofs << acrobot.xs[k][0] << ", " << acrobot.xs[k][1] << endl;
  }
  return 1;
}
