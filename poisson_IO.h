#ifndef poissonIO_h
#define poissonIO_h
class params;
class proc;
class gridsize;
class tensor0;
class tensor1;
class communicator;
#include <mpi.h>
extern "C" 
{
#include "/home/ggeraci/TOOLS/fft/fft_2d.h" //required for SANDIA parallel fft package (using FFTW)
/* required for Hypre package */
#include <math.h>
#include "_hypre_utilities.h"
#include "HYPRE_krylov.h"
#include "HYPRE.h"
#include "HYPRE_parcsr_ls.h"
}
class poisson_IO
{
  params *param_;
  proc *pc_;
  gridsize *size_;
  double epsilon_; //convergence criterion
  int length_; //length of data each process has which is asumed to be the same for all processes
  /* fft package required variables */
  fft_plan_2d *plan;
  FFT_DATA *fft_data;
  int nbuff;
  int in_ilo,in_ihi,in_jlo,in_jhi,in_klo,in_khi;
  int out_ilo,out_ihi,out_jlo,out_jhi,out_klo,out_khi;
  tensor0 fft_real; //to store fft2 of RHS
  tensor0 fft_imag;

  /* assuming constant mesh size */
  double inv_dx2;
  double inv_dy2;
  double inv_dz2;

  //Thomas algorithm: 
  double dy;
  double dz;
  double jj,kk;
  double dy2;
  double dz2;
  tensor1 thomas_c;
  tensor1 thomas_d;
  int Maximum_comm;
  int COMMUT_COUNT; //count the number of communication in the Thomas algorithm
  MPI_Request *request_send;
  MPI_Request *request_recv;
  MPI_Status *stat;
  double **send_buf;
  double **recv_buf;


  int bs_; //bordersize
  int Nx_tot_,Ny_tot_,Nz_tot_,NxNy_tot_; //global grid size
  int Nx_,Ny_,Nz_; //local grid size (including the ghost cells)
  int ilower,iupper; //range of rows "owned" by this process in the Big matrix


  tensor0 center_coeff,RHS_new; //extra memory required for Iterative Poisson solve
  tensor1 dummy;//extra memory required for Iterative Poisson solve

  void CCP_FFT(tensor0&,tensor0&); //Solve Constant Coefficient Possion using FFT
  void FFT(tensor0&);
  void IFFT(tensor0&);

  void thomas_forward_real(tensor0&,int,int);
  void thomas_backward_real(tensor0&,int,int);
  void thomas_forward_imag(tensor0&,int,int);
  void thomas_backward_imag(tensor0&,int,int);

  int num_iteration1_; //actual number of iteration required for convergence in outer loop
  double num_iteration2_; //actual number of iteration required for convergence in linear solver
  int num_iter_; //number of iteration at each lin. sys. solve
 public:
  poisson_IO(){}
  ~poisson_IO();
  poisson_IO(params*,proc*,gridsize*,communicator*);
  void Solve(tensor1&,tensor0&,tensor0&); //(C,P,RHS) div(C*grad(P))=RHS   P should includes your initial guess, Note: C is a vector quantity corresponding to coefficeints evaluated on cel faces, VCP:Variable Coefficient Poisson
  int num_iteration1() const {return num_iteration1_;}//This function returns the actual number of iteration it takes to converge in outer loop
  double num_iteration2() const {return num_iteration2_;}//This function returns the average number of iteration it takes to converge in lin. solver
  double FFT_TIME;
  double SOLVER_TIME;
};

#endif