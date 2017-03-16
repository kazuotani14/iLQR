#include "iLQR2.h"

void iLQR::boxQP(H,g,lower,upper,x0)
{
// Do this last. Fill with simpler solution first (clamp? see Tassa)
/* H: 2x2 matrix
 * g: 2x1 vector
 * lower: 2x1 vector (from control_limits)
 * upper: 2x1 vector (from control_limits)
 * x0:
 */

// TODO simpler solution (ex. clamp)

} //boxQP
