#include "MPC.h"
#include <cppad/cppad.hpp>
#include <cppad/ipopt/solve.hpp>
#include "Eigen-3.3/Eigen/Core"
//#include "matplotlibcpp.h"

using CppAD::AD;

// TODO: Set the timestep length and duration: done in MPC.h

// This value assumes the model presented in the classroom is used.
//
// It was obtained by measuring the radius formed by running the vehicle in the
// simulator around in a circle with a constant steering angle and velocity on a
// flat terrain.
//
// Lf was tuned until the the radius formed by the simulating the model
// presented in the classroom matched the previous radius.
//
// This is the length from front to CoG that has a similar radius.
//

// MPC class definition implementation.
//
MPC::MPC() {}
MPC::~MPC() {}

class FG_eval {
 public:
  // Fitted polynomial coefficients
  Eigen::VectorXd coeffs;
  size_t N;
  size_t x_start;
  size_t y_start;
  size_t psi_start;
  size_t v_start;
  size_t cte_start;
  size_t epsi_start;
  size_t sangle_start;
  size_t acc_start;
  double Lf;
  double v_target;
  double dt;
  
  FG_eval(Eigen::VectorXd coeffs, size_t N, double Lf,double v_target,double dt) {
    this->coeffs = coeffs;
    this->N = N;
    this->Lf = Lf;
    this->v_target = v_target;
    this->dt = dt;

    x_start     = 0;
    y_start     = x_start + N;
    psi_start   = y_start + N;
    v_start     = psi_start + N;
    cte_start   = v_start + N;
    epsi_start  = cte_start + N;
    sangle_start = epsi_start + N;
    acc_start     = sangle_start + N - 1;

  }

  typedef CPPAD_TESTVECTOR(AD<double>) ADvector;
  void operator()(ADvector& fg, const ADvector& vars) {
    // TODO: implement MPC
    // `fg` a vector of the cost constraints, `vars` is a vector of variable values (state & actuators)
    // NOTE: You'll probably go back and forth between this function and
    // the Solver function below.

    // The cost is stored is the first element of `fg`.
    // Any additions to the cost should be added to `fg[0]`.
    fg[0] = 0;
    // Reference State Cost
    // TODO: Define the cost related the reference state and
    // any anything you think may be beneficial.


    for (uint t = 0; t < N ; t++) {
      AD<double> cte = vars[cte_start + t];
      AD<double> epsi = vars[epsi_start + t];
      AD<double> x   = vars[x_start + t];
      AD<double> y   = vars[y_start + t];
      AD<double> v   = vars[v_start + t];

      //the final destination point
      AD<double> xf  = vars[y_start - 1];
      AD<double> yf  = vars[psi_start - 1];
      //distance from destination point (squared)
      AD<double> dist2 = CppAD::pow(x - xf,2) + CppAD::pow(y - yf,2);
      //velocity error
      AD<double> verror = v - v_target;

      AD<double> fg_next = 0;
      fg_next +=  10 * CppAD::pow(cte,2);//cte cross track error: range is 0 - 3 meters
      fg_next += 200 * CppAD::pow(epsi,2);//epsi angle error; realistic range is 0 - 0.1 radians
      ////fg_next += 0.1 * dist2; //distance from farthest point; range is ~100 meters.  Excluded for now- speed target more appropriate
      fg_next +=   2 * CppAD::pow(verror,2);//velocity error; normal range is 0 - 5 kph (need to add integral term to hit target)
      if(t < N-1){
        AD<double> sangle = vars[sangle_start + t]; //steering angle
        AD<double> acc = vars[acc_start + t]; //acceleration

        fg_next += 5  * CppAD::pow(sangle,2); //minimize steering angle; range is 0-0.2 radians- penalize large steering angles much more than small ones
        fg_next += 5  * CppAD::pow(acc,2); //minimize acceleration; range is 0-0.3
      
        if (t < N-2){
          AD<double> sangle1 = vars[sangle_start + t + 1];
          AD<double> acc1 = vars[acc_start + t + 1];
          AD<double> d_sangle = sangle1 - sangle;
          AD<double> d_acc = acc1 - acc;

          fg_next += 100 * CppAD::pow(d_sangle,2);
          fg_next +=  10 * CppAD::pow(d_acc,2);
        }//if t < N-2
      }//if t < N-1

      //fg[0] += fg_next / log(t+2); //weigh far points less than near ones, to reduce jumping when new points are added
      fg[0] += fg_next; //bad idea- never mind
    }// for t over N
    
    //
    // Setup Constraints
    //
    // NOTE: In this section you'll setup the model constraints.

    // Initial constraints
    //
    // We add 1 to each of the starting indices due to cost being located at
    // index 0 of `fg`.
    // This bumps up the position of all the other values.
    fg[1 + x_start] = vars[x_start];
    fg[1 + y_start] = vars[y_start];
    fg[1 + psi_start] = vars[psi_start];
    fg[1 + v_start] = vars[v_start];
    fg[1 + cte_start] = vars[cte_start];
    fg[1 + epsi_start] = vars[epsi_start];

    // The rest of the constraints
    for (uint t = 1; t < N ; t++) {
      AD<double> x0 = vars[x_start + t - 1];
      AD<double> x1 = vars[x_start + t];
      AD<double> y0 = vars[y_start + t - 1];
      AD<double> y1 = vars[y_start + t];
      AD<double> psi0 = vars[psi_start + t - 1];
      AD<double> psi1 = vars[psi_start + t];
      AD<double> v0 = vars[v_start + t - 1];
      AD<double> v1 = vars[v_start + t];
      AD<double> cte0 = vars[cte_start + t - 1];
      AD<double> cte1 = vars[cte_start + t];
      AD<double> epsi0 = vars[epsi_start + t - 1];
      AD<double> epsi1 = vars[epsi_start + t];

      AD<double> sangle0 = vars[sangle_start + t - 1];
      AD<double> acc0 = vars[acc_start + t - 1];

      //calcualte  y coordinate at x0
      AD<double> f0 = 0;
      AD<double> ecks = 1;
      for (uint i = 0; i < coeffs.size(); i++){
        f0 += coeffs[i] * ecks;
        ecks *= x0;
      }
      
      //calculate  derivative of track at x0: eg 'ax^3 + bx^2 + cx + d' becomes '3ax^2 + 2bx + c'
      AD<double> deriv = 0;
      AD<double> mult = 1;
      ecks = 1;
      for (uint i = 1; i <coeffs.size(); i++){ //start at 1 to ignore constant term ('d', above)
        deriv += coeffs[i] * ecks * mult;
        ecks *= x0;
        mult += 1;
      }
      //convert derivative into  angle
      AD<double> psides0 = CppAD::atan(deriv);

      // NOTE: The use of `AD<double>` and use of `CppAD`!
      // This is also CppAD can compute derivatives and pass
      // these to the solver.

      // TODO: Setup the rest of the model constraints
      fg[1 + x_start + t]    = x1 - (x0 + v0 * CppAD::cos(psi0) * dt);
      fg[1 + y_start + t]    = y1 - (y0 + v0 * CppAD::sin(psi0) * dt);
      //update for different convention of simulator
      //fg[1 + psi_start + t]  = psi1 - (psi0 + (v0/Lf) *  sangle0 * dt);
      fg[1 + psi_start + t]  = psi1 - (psi0 - (v0/Lf) *  sangle0 * dt);
      fg[1 + v_start + t]    = v1 - (v0 + acc0 * dt);
      //fg[1 + cte_start + t]  = cte1 - (cte0 + v0 * sin(epsi0) * dt);
      fg[1 + cte_start + t] = cte1 - ((f0 - y0) + (v0 * CppAD::sin(epsi0) * dt));
      //fg[1 + epsi_start + t] = epsi1 - (epsi0 + v0/Lf * sangle0 * dt);
      fg[1 + epsi_start + t] = epsi1 - ((psi0 - psides0) + (v0/Lf) * sangle0 * dt);
    }
  }
};


void MPC::init(size_t N_in, double dt_in, int n_latency_in,double v_target_in){
  N = N_in;
  dt = dt_in;
  n_latency = n_latency_in;
  v_target = v_target_in;

  debug = false;//DONE: set to false for real runs
  
  x_start     = 0;
  y_start     = x_start + N;
  psi_start   = y_start + N;
  v_start     = psi_start + N;
  cte_start   = v_start + N;
  epsi_start  = cte_start + N;
  sangle_start = epsi_start + N;
  acc_start     = sangle_start + N - 1;

  x_target.resize(N);
  y_target.resize(N);

  sangle_hist.resize(n_latency);
  acc_hist.resize(n_latency);
  for(uint i = 0; i < n_latency; i++){
    sangle_hist[i] = 0;
    acc_hist[i] = 0;
  }

  Lf = 2.67;
  step_count = 0;
  cout << "N: "<<N<<"\n";
  cout << "dt: "<<dt<<"\n";
  cout << "n_latency: "<<n_latency<<"\n";
  cout << "v_target: "<<v_target<<"\n";
  cout << "x_start: "<<x_start<<"\n";
  cout << "y_start: "<<y_start<<"\n";
  cout << "psi_start: "<<psi_start<<"\n";
  cout << "v_start: "<<v_start<<"\n";
  cout << "cte_start: "<<cte_start<<"\n";
  cout << "epsi_start: "<<epsi_start<<"\n";
  cout << "sangle_start: "<<sangle_start<<"\n";
  cout << "acc_start: "<<acc_start<<"\n";
}


vector<double> MPC::Solve(Eigen::VectorXd state, Eigen::VectorXd coeffs) {
  //cout <<"starting Solve\n";
  bool ok = true;
  size_t i;
  typedef CPPAD_TESTVECTOR(double) Dvector;
  //cout <<"setting state names\n";
  double x    = state[0];
  double y    = state[1];
  double psi  = state[2];
  double v    = state[3];
  double cte  = state[4];
  double epsi = state[5];
  
  // TODO: Set the number of model variables (includes both states and inputs).
  // For example: If the state is a 4 element vector, the actuators is a 2
  // element vector and there are 10 timesteps. The number of variables is:
  //
  // 4 * 10 + 2 * 9
  //ajg: 6 state variables, 2 actuators
  size_t n_vars = (6 * N) + (2 * (N - 1));
  // TODO: Set the number of constraints
  size_t n_constraints = N * 6;
  //cout <<"n_vars,n_constraints are "<<n_vars<<","<<n_constraints<<"\n";
  // Initial value of the independent variables.
  // SHOULD BE 0 besides initial state.
  Dvector vars(n_vars);
  //cout <<"initializing vars\n";
  for (i = 0; i < n_vars; i++) {
    vars[i] = 0.0;
  }
  //cout <<"initialized vars\n";
  

  vars[x_start] = x;
  vars[y_start] = y;
  vars[psi_start] = psi;
  vars[v_start] = v;
  vars[cte_start] = cte;
  vars[epsi_start] = epsi;

  Dvector vars_lowerbound(n_vars);
  Dvector vars_upperbound(n_vars);
  // TODO: Set lower and upper limits for variables.
  for (i = 0; i < sangle_start; i++) {
    //cout <<"setting large bound "<<i<<"\n";
    vars_lowerbound[i] = -1.0e19;
    vars_upperbound[i] = 1.0e19;
  }
  //constrain steering angles that haven't been applied yet (due to latency) to the values that were already sent
  uint j = 0;
  for (i = sangle_start; i < sangle_start + n_latency; i++){
    vars_lowerbound[i] = sangle_hist[j];
    vars_upperbound[i] = sangle_hist[j];
    j++;
  }
  // The upper and lower limits of sangle are set to -25 and 25
  // degrees (values in radians).
  // Consider changing to something smaller, if desired
  //    (would mean I'd rather have a larger cte than a large turn value
  //     I'd rather attack that through cost function fudge factors)
  for (i = sangle_start + n_latency; i < acc_start; i++) {
    //cout <<"setting sangle bound "<<i<<"\n";
    vars_lowerbound[i] = -0.436332;
    vars_upperbound[i] = 0.436332;
  }
  j = 0;
  //constrain acceleration values that haven't been applied yet (due to latency) to the values that were already sent
  for (i = acc_start; i < acc_start + n_latency; i++){
    vars_lowerbound[i] = acc_hist[j];
    vars_upperbound[i] = acc_hist[j];
    j++;
  }
  //constrain other accelerations to the allowed range
  for (i = acc_start + n_latency; i < n_vars; i++) {
    //cout <<"setting a bound "<<i<<"\n";
    vars_lowerbound[i] = -1.0;
    vars_upperbound[i] = 1.0;
  }
  //cout <<"finished setting bounds\n";
  // Lower and upper limits for the constraints
  // Should be 0 besides initial state.
  Dvector constraints_lowerbound(n_constraints);
  Dvector constraints_upperbound(n_constraints);
  for (i = 0; i < n_constraints; i++) {
    constraints_lowerbound[i] = 0;
    constraints_upperbound[i] = 0;
  }
  constraints_lowerbound[x_start] = x;
  constraints_lowerbound[y_start] = y;
  constraints_lowerbound[psi_start] = psi;
  constraints_lowerbound[v_start] = v;
  constraints_lowerbound[cte_start] = cte;
  constraints_lowerbound[epsi_start] = epsi;

  constraints_upperbound[x_start] = x;
  constraints_upperbound[y_start] = y;
  constraints_upperbound[psi_start] = psi;
  constraints_upperbound[v_start] = v;
  constraints_upperbound[cte_start] = cte;
  constraints_upperbound[epsi_start] = epsi;
  
  // object that computes objective and constraints
  //cout <<"creating fg_eval\n";
  FG_eval fg_eval(coeffs,N,Lf,v_target,dt);

  //
  // NOTE: You don't have to worry about these options
  //
  // options for IPOPT solver
  std::string options;
  // Uncomment this if you'd like more print information
  options += "Integer print_level  0\n";
  // NOTE: Setting sparse to true allows the solver to take advantage
  // of sparse routines, this makes the computation MUCH FASTER. If you
  // can uncomment 1 of these and see if it makes a difference or not but
  // if you uncomment both the computation time should go up in orders of
  // magnitude.
  options += "Sparse  true        forward\n";
  options += "Sparse  true        reverse\n";
  // NOTE: Currently the solver has a maximum time limit of 0.5 seconds.
  // Change this as you see fit.
  options += "Numeric max_cpu_time          0.5\n";

  // place to return solution
  CppAD::ipopt::solve_result<Dvector> solution;

  // solve the problem
  //cout <<"calling solve\n";
  CppAD::ipopt::solve<Dvector, FG_eval>(
      options, vars, vars_lowerbound, vars_upperbound, constraints_lowerbound,
      constraints_upperbound, fg_eval, solution);

  // Check some of the solution values
  ok &= solution.status == CppAD::ipopt::solve_result<Dvector>::success;

  // Cost
  //auto cost = solution.obj_value;
  //std::cout << "Cost " << cost << std::endl;

  //c++ isn't as easy as python
  //x_target = solution.x[x_start:y_start - 1];
  //y_target = solution.x[y_start:psi_start - 1];
  for (i = 0; i < N; i++){
    x_target[i] = solution.x[x_start+i];
    y_target[i] = solution.x[y_start+i];
  }
  
  // TODO: Return the first actuator values. The variables can be accessed with
  // `solution.x[i]`.
  //
  // {...} is shorthand for creating a vector, so auto x1 = {1.0,2.0}
  // creates a 2 element double vector.
  //To account for latency, return the step far enough in to be applicable in <latency> amount of time

  double found_sangle = solution.x[sangle_start+n_latency];
  double found_acc = solution.x[acc_start+n_latency];
  if (debug && (step_count % 100 == 0)){
    found_acc = -1; //create a glitch that I can trace through the logs to measure actual latency
  }
  
  for (i=0; i < n_latency - 1; i++){
    sangle_hist[i] = sangle_hist[i+1];
    acc_hist[i] = acc_hist[i+1];
  }
  sangle_hist[n_latency - 1] = found_sangle;
  acc_hist[n_latency - 1] = found_acc;
  step_count++;
  return {found_sangle, found_acc};
}

vector<double> MPC::get_x_target(){
  return x_target;
}

vector<double> MPC::get_y_target(){
  return y_target;
}
