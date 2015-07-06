#ifndef grid_h
#define grid_h
class proc;
class params;
class gridsize;
class communicator;

#include "tensor0.h"
#include "tensor1.h"
#include "particle.h"
#include <mpi.h>
#include "poisson.h"
#include "poisson_IO.h"
#include <fstream>
#include <iostream>
class grid
{
  params *param_;
  proc *pc_;
  gridsize *size_;
  communicator *com_;
  double RK4_preCoeff[4];
  double RK4_postCoeff[4];
  double mean_energy_transferred;
  double C_BC; //the convective velocity used in outflow BC
  std::fstream touch_check;
  std::ofstream stat_TKE;
  std::ofstream stat_CMax;
  std::ofstream stat_CMin;
  std::ofstream stat_CMean;
  std::ofstream stat_RhoMax;
  std::ofstream stat_RhoMin;
  std::ofstream stat_RhoMean;
  std::ofstream stat_ParticleMaxCFL;
  std::ofstream stat_GasMaxCFL;
  std::ofstream stat_GasMaxDiffCFL;
  std::ofstream stat_NumIteration;
  std::ofstream stat_BalanceIndex;
  std::ofstream stat_Uin;
  std::ofstream stat_Uout;
  std::ofstream stat_dx_o_eta;
  std::ofstream stat_ST;
  std::ofstream stat_time;
  bool Is_touch_; //=1 if touch=1 =0 if not
  int sign_fnc(double); //=1 if double>=0 otherwise is -1
  double ABS(double); //=|double|
  void Write_info(); //Wrtie information to info.txt
  void open_stat_file(char*,std::ofstream&); //open a stat file
  void store_stat_binary(char *,tensor0&); //get the name of file and store the reduced value using communicator

  poisson PS_;
  poisson_IO PS_IO_;
  //SOME STATISTICS
  double U_avg_in; 
  double U_avg_out;
 public:
  particle part; //particle tracking class
  tensor1 RU; //rho*u stored at cell faces
  tensor1 RU_int; //previous RK4 substep value
  tensor1 RU_new; //next RK4 substep value
  tensor1 RU_np1; //next timestep value
  tensor1 RU_WP; //momentum with pressure gradient effect (should match the divergence condition)
  tensor1 RHS_RU;
  tensor1 U; //u stored at cell center
  tensor0 P; //pressure stored at cell center
  tensor0 dP; //pressure delta form (used in Poisson equation) stored at cell center
  tensor0 RHS_Pois; //one part of RHS of Poisson equation
  tensor0 C; //particle concentration stored at cell center
  tensor0 Rho; //Gas density stored at cell center
  tensor0 Rho_int;
  tensor0 Rho_new;
  tensor0 Rho_np1;
  tensor0 RHS_Rho;
  tensor1 Rho_face; //Gas density stored at cell faces
  tensor0 T; //Gas temperature stored at cell cneter
  tensor1 dummy; //dummy array for computations
  tensor1 dummy2; //dummy array for computations
  tensor1 dummy3; //dummy array for computations
  tensor0 divergence; //to store divergence of momentum/velocity
  tensor0 RHS_Part_Temp; //Interpolated RHS of particle energy equation (due to the algorithm it has to be saved)
  double P0; // mean thermodynamic pressure
  double P0_int;
  double P0_new;
  double P0_np1;
  double dP0_dt;    // P0 rate of change
  double sigma_RHS;   // sum of RHS of Poisson equation over all the domain exept for the P0 rate of change term
  double T_cur;
  int num_timestep;
  int RK4_count;
  grid(){}
  grid(gridsize*,params*,proc*,communicator*);
  ~grid();
  void Initialize(grid&); //the argument is usefull only for IO case!
  void Store();
  void TimeAdvance(); //advance one time step
  void Update_Rho();
  void Update_Rho(grid&); //For Inflow Outflow simulation
  void Update_RU_WOP();
  void Update_RU_WOP(grid&); //For InFlow OutFlow simulation
  void Update_P0();
  void Compute_Div_U_new();
  void Compute_RHS_Pois();
  void Solve_Poisson();
  void Solve_Poisson_IO();
  void Update_RU_WP();
  void TimeAdvance_RK4();
  void Test_Poisson();
  void Test_Poisson_IO();
  void Update_Particle();
  void Update_Particle(grid&); //For InFlow OutFlow Simulation
  void Statistics();
  bool Touch(grid&); //=0/1 in touch.check (1: exit code now, 0: contunue), argument is reference to IO simulation
  void Correct_mean_RU(); //make mean zero in y and z directions, and U_INFTY in x direction
  void Correct_OutFlow(); //make sure sum(div(u_new_WOP))=sum(div(u_new)) by changing outflow!
};
#endif               
