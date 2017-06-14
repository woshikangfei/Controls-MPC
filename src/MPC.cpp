#include "MPC.h"
#include <cppad/cppad.hpp>
#include <cppad/ipopt/solve.hpp>
#include "Eigen-3.3/Eigen/Core"

using CppAD::AD;

// TODO: Set the timestep length and duration
size_t N = 10;
double dt = 0.06;

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
const double Lf = 2.67;

double ref_cte = 0;
double ref_epsi = 0;
double ref_v = 70;

size_t 
_x = 0;
size_t start_y = start_x + N;
size_t start_psi = start_y + N;
size_t start_v = start_psi + N;
size_t start_cte = start_v + N;
size_t start_epsi = start_cte + N;
size_t start_delta = start_epsi + N;
size_t start_a = start_delta + N -1;


class FG_eval {
 public:
  // Fitted polynomial coefficients
  Eigen::VectorXd coeffs;
  FG_eval(Eigen::VectorXd coeffs) { this->coeffs = coeffs; }

  typedef CPPAD_TESTVECTOR(AD<double>) ADvector;
  void operator()(ADvector& fg, const ADvector& vars) {

   
    fg[0] = 0;

    for(int i = 0; i < N; i++) {
      fg[0] += CppAD::pow(vars[start_cte + i] - ref_cte, 2);
      fg[0] += CppAD::pow(vars[start_epsi + i] - ref_epsi, 2);
      fg[0] += CppAD::pow(vars[start_v + i] - ref_v, 2);
    }

    
    for(int i = 0; i < N-1; i++) {
      fg[0] += CppAD::pow(vars[start_delta + i], 2);
      fg[0] += CppAD::pow(vars[start_a + i], 2);
      fg[0] += CppAD::pow(vars[start_cte+i+1] - vars[start_cte+i], 2);
    }

    
    for(int i = 0; i < N-2; i++) {
      fg[0] += 500*CppAD::pow(vars[start_delta + i+ 1] - vars[start_delta + i], 2);
      fg[0] += CppAD::pow(vars[start_a + i + 1] - vars[start_a + i], 2);
    }

    // initial contraints
    fg[start_x + 1]    = vars[start_x];
    fg[start_y + 1]    = vars[start_y];
    fg[start_psi + 1]  = vars[start_psi];
    fg[start_v + 1]    = vars[start_v];
    fg[start_cte + 1]  = vars[start_cte];
    fg[start_epsi + 1] = vars[start_epsi];

    // rest of contraints to zeros
    for(int i = 0; i < N -1; i++) {
      
      AD<double> x1 = vars[start_x + i +1];
      AD<double> y1 = vars[start_y + i +1];
      AD<double> psi1 = vars[start_psi + i +1];
      AD<double> v1 = vars[start_v + i +1];
      AD<double> cte1 = vars[start_cte + i +1];
      AD<double> epsi1 = vars[start_epsi + i +1];

      
      AD<double> x0 = vars[start_x + i];
      AD<double> y0 = vars[start_y + i];
      AD<double> psi0 = vars[start_psi + i];
      AD<double> v0 = vars[start_v + i];
      AD<double> cte0 = vars[start_cte + i];
      AD<double> epsi0 = vars[start_epsi + i];

      
      AD<double> delta0 = vars[start_delta +i];
      AD<double> a0 = vars[start_a +i];

      AD<double> f0 = coeffs[0]+coeffs[1]*x0+coeffs[2]*CppAD::pow(x0, 2)+coeffs[3]*CppAD::pow(x0, 3);
      AD<double> psides0 = CppAD::atan(coeffs[1]+2*coeffs[2]*x0+3*coeffs[3]*CppAD::pow(x0, 2));

      
      fg[2 + start_x + i] = x1 - (x0 + v0*CppAD::cos(psi0)*dt);
      fg[2 + start_y + i] = y1 - (y0 + v0*CppAD::sin(psi0)*dt);
      fg[2 + start_psi + i] = psi1 - (psi0 + (v0*delta0/Lf)*dt);
      fg[2 + start_v + i] = v1 - (v0 + a0*dt);
      fg[2 + start_cte + i] = cte1 - (f0 - y0 + v0*CppAD::sin(epsi0)*dt);
      fg[2 + start_epsi + i] = epsi1 - ((psi0 - psides0) + (v0*delta0/Lf)*dt);
    }
   
  }
};

//
// MPC class definition implementation.
//
MPC::MPC() {
  solx_.resize(N);
  soly_.resize(N);
}
MPC::~MPC() {}

vector<double> MPC::Solve(Eigen::VectorXd state, Eigen::VectorXd coeffs) {
  bool ok = true;
  size_t i;
  typedef CPPAD_TESTVECTOR(double) Dvector;

  // TODO: Set the number of model variables (includes both states and inputs).
  // For example: If the state is a 4 element vector, the actuators is a 2
  // element vector and there are 10 timesteps. The number of variables is:
  //
  // 4 * 10 + 2 * 9
  size_t n_vars = 6*N + 2*(N-1);
  // TODO: Set the number of constraints
  size_t n_constraints = 6*N;
 

  double x    = state[0];
  double y    = state[1];
  double psi  = state[2];
  double v    = state[3];
  double cte  = state[4];
  double epsi = state[5];
 

  // Initial value of the independent variables.
  // SHOULD BE 0 besides initial state.
  Dvector vars(n_vars);
  for (int i = 0; i < n_vars; i++) {
    vars[i] = 0;
  }

  vars[start_x]    = x;
  vars[start_y]    = y;
  vars[start_psi]  = psi;
  vars[start_v]    = v;
  vars[start_cte]  = cte;
  vars[start_epsi] = epsi;
   
  Dvector vars_lowerbound(n_vars);
  Dvector vars_upperbound(n_vars);
  // TODO: Set lower and upper limits for variables.
  // Set lower and upper limits for variables.
  for(int i = 0; i < start_delta; i++) {
    vars_lowerbound[i] = -1.0e19;
    vars_upperbound[i] = 1.0e19;
  }

  
  for(int i = start_delta; i < start_a; i++) {
    vars_lowerbound[i] = -0.436332; 
    vars_upperbound[i] = 0.436332;
  }

  //   upper and lower limits
  for(int i = start_a; i < n_vars; i++) {
    vars_lowerbound[i] = -1.0;
    vars_upperbound[i] = 1.0;
  }
 
 
  // Lower and upper limits for the constraints
  // Should be 0 besides initial state.
  Dvector constraints_lowerbound(n_constraints);
  Dvector constraints_upperbound(n_constraints);
  for (int i = 0; i < n_constraints; i++) {
    constraints_lowerbound[i] = 0;
    constraints_upperbound[i] = 0;
  }
 
  constraints_lowerbound[start_x] = x;
  constraints_lowerbound[start_y] = y;
  constraints_lowerbound[start_psi] = psi;
  constraints_lowerbound[start_v] = v;
  constraints_lowerbound[start_cte] = cte;
  constraints_lowerbound[start_epsi] = epsi;

  constraints_upperbound[start_x] = x;
  constraints_upperbound[start_y] = y;
  constraints_upperbound[start_psi] = psi;
  constraints_upperbound[start_v] = v;
  constraints_upperbound[start_cte] = cte;
  constraints_upperbound[start_epsi] = epsi;

  // object that computes objective and constraints
  FG_eval fg_eval(coeffs);

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
  CppAD::ipopt::solve<Dvector, FG_eval>(
      options, vars, vars_lowerbound, vars_upperbound, constraints_lowerbound,
      constraints_upperbound, fg_eval, solution);

  // Check some of the solution values
  ok &= solution.status == CppAD::ipopt::solve_result<Dvector>::success;

  // Cost
  auto cost = solution.obj_value;
  std::cout << "Cost " << cost << std::endl;

  // TODO: Return the first actuator values. The variables can be accessed with
  // `solution.x[i]`.
  //
  // {...} is shorthand for creating a vector, so auto x1 = {1.0,2.0}
  // creates a 2 element double vector.
  return { solution.x[start_delta+1], solution.x[start_a+1] };
}
