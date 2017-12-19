#ifndef MPC_H
#define MPC_H

#include <vector>
#include "Eigen-3.3/Eigen/Core"

using namespace std;

class MPC {
  size_t N;//number of timesteps to evaluate
  double dt;//timestep in seconds
  uint n_latency;//control latency in timesteps
  vector<double> sangle_hist;//store previously sent values of sangle, to apply as constraints in future steps
  vector<double> acc_hist;//store previously sent values of a, to apply as constraints in future steps
  double v_target;
  long step_count; //keep track of number of steps, for debugging statements
  bool debug;//turn on debugging switches, including a glitch of a=-1 in order to measure actual latency
  
  //moving these here so they can depend on values received in init
  size_t x_start;
  size_t y_start;
  size_t psi_start;
  size_t v_start;
  size_t cte_start;
  size_t epsi_start;
  size_t sangle_start;
  size_t acc_start;
  double Lf;

  //store target points, and use getters to access
  vector<double> x_target;
  vector<double> y_target;
  
 public:  
  MPC();

  virtual ~MPC();

  //initialize constants, indicies, sangle and throttle (a) history (for latency compensation)
  void init(size_t N_in, double dt_in, int latency_in, double v_target_in);
  // Solve the model given an initial state and polynomial coefficients.
  // Return the first actuatotions.
  vector<double> Solve(Eigen::VectorXd state, Eigen::VectorXd coeffs);

  //separate getter functions to return predicted points' x and y ordinate arrays
  //only for displaying green line
  vector<double> get_x_target();
  vector<double> get_y_target();
};

#endif /* MPC_H */
