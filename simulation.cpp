#include <iostream>
#include <mpi.h>
#include "params.h"
#include "proc.h"
#include "gridsize.h"
#include "tensor0.h"
#include "tensor1.h"
#include "communicator.h"
#include "poisson.h"
#include "poisson_IO.h"
#include "grid.h"
int main (int argc,char *argv[] ) 
{
  if (argc != 3) 
    {
      std::cout << "Please supply parameter files!" << std::endl;
      exit(1);
    }
  MPI_Init(&argc, &argv);
  params PARAM(argv[1]);
  params PARAM_IO(argv[2]);
  proc PROC;
  gridsize GSIZE(&PARAM,&PROC);
  gridsize GSIZE_IO(&PARAM_IO,&PROC);
  communicator COMM(&GSIZE,&PARAM,&PROC);
  communicator COMM_IO(&GSIZE_IO,&PARAM_IO,&PROC);
  grid GRID(&GSIZE,&PARAM,&PROC,&COMM);
  grid GRID_IO(&GSIZE_IO,&PARAM_IO,&PROC,&COMM_IO);
  GRID.Initialize(GRID);
  GRID_IO.Initialize(GRID);
  
  //GRID_IO.Test_Poisson_IO();

  do {
    MPI_Barrier(MPI_COMM_WORLD);
      PARAM.update(GRID.T_cur);
      PARAM_IO.update(GRID_IO.T_cur); //FOR IO
      for (GRID.RK4_count=0;GRID.RK4_count<4;GRID.RK4_count++)
	{
	  //First time advance box
	  GRID.Update_Rho();
	  GRID.Update_P0();
	  GRID.Update_Particle();
	  GRID.Update_RU_WOP(); 
	  GRID.Compute_Div_U_new();
	  GRID.Compute_RHS_Pois();
	  GRID.Solve_Poisson();
	  GRID.Update_RU_WP();
	  //Now time advance InFlow-OutFlow
	  GRID_IO.RK4_count=GRID.RK4_count;
	  GRID_IO.Update_Rho(GRID); //RU_int, Rho, Rho_int, and Rho_np1  ===>>>>   Rho_np1 , Rho_new
	  GRID_IO.Update_P0();
	  GRID_IO.Update_Particle(GRID); //RU_int, Rho_face
	  GRID_IO.Update_RU_WOP(GRID); //RHS_RU (from particle update or zero), U@face(from particle), RU_int, RU_np1, Rho_face ====>>>> RU_new, RU_np1 (expect inflow ghost cells for u) 
	  GRID_IO.Compute_Div_U_new(); //Rho_new , RHS_Part_Temp ====>>>> divergence (only for interior cells)
	  GRID_IO.Correct_OutFlow(); // Rho_new ====>>>>> Rho_face   Update RU_new and RU_np1 at the outlet
	  GRID_IO.Compute_RHS_Pois(); //Rho_face, RU_new , divergence ====>>> RHS_POIS = div(grad(RU_new/Rho_face)) - divergence
	  GRID_IO.Solve_Poisson_IO();
	  GRID.TimeAdvance_RK4();
	  GRID_IO.Update_RU_WP();
	  GRID_IO.TimeAdvance_RK4();
	}
      GRID.Correct_mean_RU(); //make all components zero mean
      GRID_IO.Correct_mean_RU(); //only make the y and z components zero mean
      GRID.TimeAdvance();
      GRID_IO.TimeAdvance();
      GRID.Statistics();
      GRID_IO.Statistics();
      GRID.Store();
      GRID_IO.Store();
  }  while ((GRID.T_cur<PARAM.T_final())&&(!GRID.Touch(GRID_IO)));
    MPI_Finalize();
}
  
