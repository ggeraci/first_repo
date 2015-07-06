#include "params.h"
#include "proc.h"
#include "gridsize.h"
#include "tensor0.h"
#include "tensor1.h"
#include <mpi.h>
#include <complex>
#include "communicator.h"
#include "poisson_IO.h"
#include <sys/time.h>
#include <ctime>
poisson_IO::poisson_IO(params* par,proc* pro,gridsize* gri,communicator* com):center_coeff(gri,com),dummy(gri,com),RHS_new(gri,com), fft_real(gri,com), fft_imag(gri,com), thomas_c(gri,com), thomas_d(gri,com)
{
  param_=par;
  pc_=pro;
  size_=gri;
  in_ilo=gri->il();
  in_ihi=gri->ih();
  in_jlo=gri->jl();
  in_jhi=gri->jh();
  in_klo=gri->kl();
  in_khi=gri->kh();
  out_ilo=in_ilo;
  out_ihi=in_ihi;
  out_jlo=in_jlo;
  out_jhi=in_jhi;
  out_klo=in_klo;
  out_khi=in_khi;
  bs_=size_->bs();
  Nx_tot_=gri->Nx_tot();
  Ny_tot_=gri->Ny_tot();
  Nz_tot_=gri->Nz_tot();
  NxNy_tot_=Nx_tot_*Ny_tot_;
  Nx_=gri->Nx(); Ny_=gri->Ny(); Nz_=gri->Nz();
  epsilon_=par->epsilon();
  //timing solvers:
  FFT_TIME=0;
  SOLVER_TIME=0;
  /* FFT INITIALIZE */
  nbuff=(in_jhi-in_jlo)*(in_khi-in_klo);
  plan=fft_2d_create_plan(pc_->XMates(),gri->Ny_tot(),gri->Nz_tot(),in_jlo,in_jhi,in_klo,in_khi,out_jlo,out_jhi,out_klo,out_khi,0,0,&nbuff);
  fft_data=new FFT_DATA[nbuff];
  inv_dx2=1./(size_->dx()*size_->dx());
  inv_dy2=1./(size_->dy()*size_->dy());
  inv_dz2=1./(size_->dz()*size_->dz());
  // Thomas initialize
  Maximum_comm=2*Ny_*Nz_;
  request_send=new MPI_Request[Maximum_comm];
  request_recv=new MPI_Request[Maximum_comm];
  stat=new MPI_Status[Maximum_comm];
  send_buf=new double*[Maximum_comm];
  recv_buf=new double*[Maximum_comm];
  for (int i(0);i<Maximum_comm;i++)
    {
      send_buf[i]=new double[2];
      recv_buf[i]=new double[2];
    }
  dy=size_->dy();
  dz=size_->dz();
  dy2=2./(dy*dy);
  dz2=2./(dz*dz);
}

poisson_IO::~poisson_IO()
{
  //  fft_3d_destroy_plan(plan);	
  delete[] fft_data;
  //thomas:
  for (int i(0);i<Maximum_comm;i++)
    {
      delete[] send_buf[i];
      delete[] recv_buf[i];
    }
  delete[] send_buf;
  delete[] recv_buf;
  delete[] request_send;
  delete[] request_recv;
  delete[] stat;
}

void poisson_IO::Solve(tensor1 &C,tensor0 &P,tensor0 &RHS)
{
  center_coeff.Equal_I_F2C(C); //compute average valuie of coefficient at cell center
  center_coeff.Equal_Divide(-1.,center_coeff); //compute minus inverse of cell center coefficients
  int i(0);
  num_iteration2_=0;
  double error(100);
  while((i<param_->Iteration())&&(error>epsilon_)) //iteration loop
    {
      dummy.Equal_Grad_C2F(P);
      dummy*=C;
      RHS_new.Equal_Div_F2C(dummy); //RHS_new=div(C grad(P))
      RHS_new-=RHS;
      RHS_new*=center_coeff; //RHS_new is the new RHS
      CCP_FFT(RHS_new,dummy.x); // dp is stored in dummy.x ,assumed P_final=P_guess+dp
      error=dummy.x.max_abs();
      P+=dummy.x; //correct P_guess
      P.Update_Ghosts();
      //correct P for inflow
      if (pc_->IsIn())
	for (int j(0);j<Ny_;j++)
	  for (int k(0);k<Nz_;k++)
	    P(0,j,k)=P(1,j,k);
      //correct P for outflow
      if (pc_->IsOut())
	for (int j(0);j<Ny_;j++)
	  for (int k(0);k<Nz_;k++)
	    P(Nx_-1,j,k)=P(Nx_-2,j,k);  
      i++;
      num_iteration2_+=num_iter_;
    }
  num_iteration1_=i;
  num_iteration2_/=i;
}

void poisson_IO::CCP_FFT(tensor0& RHS,tensor0 &P)
{
  clock_t startTime=clock();
  FFT(RHS);
  FFT_TIME+=double( clock() - startTime ) / (double)CLOCKS_PER_SEC;
  //here we have RHS(x,k2,k3) computed, and need to solve the tridiagonal linear system
  startTime=clock();

  // Thomas forward
  COMMUT_COUNT=0;
  for (int k=in_klo;k<=in_khi;k++)
    for (int j=in_jlo;j<=in_jhi;j++)
      {
	thomas_forward_real(fft_real,j,k);
	thomas_forward_imag(fft_imag,j,k);
     }
  //make sure all sends are done
  if (!pc_->IsOut())
    MPI_Waitall(COMMUT_COUNT,&request_send[0],&stat[0]);
  // Thomas backward
  COMMUT_COUNT=0;
  for (int k=in_klo;k<=in_khi;k++)
    for (int j=in_jlo;j<=in_jhi;j++)
      {
	thomas_backward_real(fft_real,j,k);
	thomas_backward_imag(fft_imag,j,k);
     }
  //make sure all sends are done
  if (!pc_->IsIn())
    MPI_Waitall(COMMUT_COUNT,&request_send[0],&stat[0]);

  SOLVER_TIME+=double( clock() - startTime ) / (double)CLOCKS_PER_SEC;
  // Now we have P(i,k2,k3) computed in fft_real and fft_imag 
  startTime=clock();
  IFFT(P);
  FFT_TIME+=double( clock() - startTime ) / (double)CLOCKS_PER_SEC;
}

void poisson_IO::FFT(tensor0 &RHS)
{
  //load data
  int count=0;
  for (int i=bs_;i<Nx_-bs_;i++)
    {
      count=0;
      for (int k=in_klo;k<=in_khi;k++)
	for (int j=in_jlo;j<=in_jhi;j++)
	  {
	    fft_data[count].im=0;
	    fft_data[count++].re=RHS(i,j-in_jlo+bs_,k-in_klo+bs_);
	  }
      //Take FT of RHS
      fft_2d(fft_data,fft_data,1,plan);
      //store in the 3D tensor
      count=0;
      for (int k=in_klo;k<=in_khi;k++)
	for (int j=in_jlo;j<=in_jhi;j++)
	  {
	    //store the FFTed data
	    fft_real(i,j-in_jlo+bs_,k-in_klo+bs_)=fft_data[count].re;
	    fft_imag(i,j-in_jlo+bs_,k-in_klo+bs_)=fft_data[count].im;
	    count++;
	  } 
    }
}

void poisson_IO::IFFT(tensor0 &P)
{
  int tot=size_->Ny_tot()*size_->Nz_tot();
  int count=0;
  for (int i=bs_;i<Nx_-bs_;i++)
    {
      count=0;
      for (int k=in_klo;k<=in_khi;k++)
	for (int j=in_jlo;j<=in_jhi;j++)
	  {
	    fft_data[count].im=fft_imag(i,j-in_jlo+bs_,k-in_klo+bs_);
	    fft_data[count++].re=fft_real(i,j-in_jlo+bs_,k-in_klo+bs_);
	  }
      //Take IFT for this slice
      fft_2d(fft_data,fft_data,-1,plan);
      //scale and store solution in P
      count=0;
      for (int k=in_klo;k<=in_khi;k++)
	for (int j=in_jlo;j<=in_jhi;j++)
	  {
	    //store the FFTed data
	    P(i,j-in_jlo+bs_,k-in_klo+bs_)=fft_data[count++].re/(tot);
	  }
    }
}

void poisson_IO::thomas_forward_real(tensor0 &val,int j,int k)
{
  //compute wavenumbers
  double TWO_PI(2*3.141592653589793);
  double Two_PI_Over_Ly = TWO_PI/size_->Ly();
  double Two_PI_Over_Lz = TWO_PI/size_->Lz();
  jj = ( ( j<Ny_tot_/2 ) ? j : j-Ny_tot_ ) * Two_PI_Over_Ly;
  kk = ( ( k<Nz_tot_/2 ) ? k : k-Nz_tot_ ) * Two_PI_Over_Lz;

  //local j and k indecis
  int jl=j-in_jlo+bs_;
  int kl=k-in_klo+bs_;
  
  double a=inv_dx2;
  double b0=-((dy2)*(1-cos(jj*dy))+(dz2)*(1-cos(kk*dz)))-inv_dx2;
  double b=-((dy2)*(1-cos(jj*dy))+(dz2)*(1-cos(kk*dz)))-2*inv_dx2;
  double c=inv_dx2;

  if (pc_->IsIn()) //first slice
    {
      thomas_c.x(bs_,jl,kl)=c/b0;
      thomas_d.x(bs_,jl,kl)=val(bs_,jl,kl)/b0;
      if (!pc_->IsOut()) //ends in middle
	{
	  for (int i(bs_+1);i<Nx_-bs_;i++)
	    {
	      thomas_c.x(i,jl,kl)=c/(b-thomas_c.x(i-1,jl,kl)*a);
	      thomas_d.x(i,jl,kl)=(val(i,jl,kl)-thomas_d.x(i-1,jl,kl)*a)/(b-thomas_c.x(i-1,jl,kl)*a);
	    }
	  //Send the result to the right process
	  send_buf[COMMUT_COUNT][0]=thomas_c.x(Nx_-bs_-1,jl,kl);
	  send_buf[COMMUT_COUNT][1]=thomas_d.x(Nx_-bs_-1,jl,kl);
	  MPI_Isend(send_buf[COMMUT_COUNT],2,MPI_DOUBLE,pc_->RIGHT(),COMMUT_COUNT,MPI_COMM_WORLD,&request_send[COMMUT_COUNT]);
	  COMMUT_COUNT++;
	}
      else //ends int the last slice
	{
	  for (int i(bs_+1);i<Nx_-bs_-1;i++)
	    {
	      thomas_c.x(i,jl,kl)=c/(b-thomas_c.x(i-1,jl,kl)*a);
	      thomas_d.x(i,jl,kl)=(val(i,jl,kl)-thomas_d.x(i-1,jl,kl)*a)/(b-thomas_c.x(i-1,jl,kl)*a);
	    }
	  if (j+k==0)
	    thomas_d.x(Nx_-bs_-1,jl,kl)=0;
	  else
	    thomas_d.x(Nx_-bs_-1,jl,kl)=(val(Nx_-bs_-1,jl,kl)-thomas_d.x(Nx_-bs_-2,jl,kl)*a)/(b0-thomas_c.x(Nx_-bs_-2,jl,kl)*a);
	}
    }
  else // middle slice
    {
      MPI_Recv(recv_buf[COMMUT_COUNT],2,MPI_DOUBLE,pc_->LEFT(),COMMUT_COUNT,MPI_COMM_WORLD,&stat[COMMUT_COUNT]);
      thomas_c.x(0,jl,kl)=recv_buf[COMMUT_COUNT][0];
      thomas_d.x(0,jl,kl)=recv_buf[COMMUT_COUNT][1];
      if (!pc_->IsOut()) //ends in the middle
	{
	  for (int i(bs_);i<Nx_-bs_;i++)
	    {
	      thomas_c.x(i,jl,kl)=c/(b-thomas_c.x(i-1,jl,kl)*a);
	      thomas_d.x(i,jl,kl)=(val(i,jl,kl)-thomas_d.x(i-1,jl,kl)*a)/(b-thomas_c.x(i-1,jl,kl)*a);
	    }
	  //Send the result to the right process
	  send_buf[COMMUT_COUNT][0]=thomas_c.x(Nx_-bs_-1,jl,kl);
	  send_buf[COMMUT_COUNT][1]=thomas_d.x(Nx_-bs_-1,jl,kl);
	  MPI_Isend(send_buf[COMMUT_COUNT],2,MPI_DOUBLE,pc_->RIGHT(),COMMUT_COUNT,MPI_COMM_WORLD,&request_send[COMMUT_COUNT]);
	}
      else //ends in the last slice
	{
	  for (int i(bs_);i<Nx_-bs_-1;i++)
	    {
	      thomas_c.x(i,jl,kl)=c/(b-thomas_c.x(i-1,jl,kl)*a);
	      thomas_d.x(i,jl,kl)=(val(i,jl,kl)-thomas_d.x(i-1,jl,kl)*a)/(b-thomas_c.x(i-1,jl,kl)*a);
	    }
	  if (j+k==0)
	    thomas_d.x(Nx_-bs_-1,jl,kl)=0;
	  else
	    thomas_d.x(Nx_-bs_-1,jl,kl)=(val(Nx_-bs_-1,jl,kl)-thomas_d.x(Nx_-bs_-2,jl,kl)*a)/(b0-thomas_c.x(Nx_-bs_-2,jl,kl)*a);
	}
      COMMUT_COUNT++;
    }
}

void poisson_IO::thomas_backward_real(tensor0 &val,int j,int k)
{
  //compute wavenumbers
  double TWO_PI(2*3.141592653589793);
  double Two_PI_Over_Ly = TWO_PI/size_->Ly();
  double Two_PI_Over_Lz = TWO_PI/size_->Lz();
  jj = ( ( j<Ny_tot_/2 ) ? j : j-Ny_tot_ ) * Two_PI_Over_Ly;
  kk = ( ( k<Nz_tot_/2 ) ? k : k-Nz_tot_ ) * Two_PI_Over_Lz;

  //local j and k indecis
  int jl=j-in_jlo+bs_;
  int kl=k-in_klo+bs_;
  
  if (pc_->IsOut()) //the last process does not recv data for forward step
    {
      val(Nx_-bs_-1,jl,kl)=thomas_d.x(Nx_-bs_-1,jl,kl);
      for (int i(Nx_-bs_-2);i>=bs_;i--)
	val(i,jl,kl)=thomas_d.x(i,jl,kl)-thomas_c.x(i,jl,kl)*val(i+1,jl,kl);	 
      if (!pc_->IsIn())
	{
	  //Send the result to the right process
	  send_buf[COMMUT_COUNT][0]=val(bs_,jl,kl);
	  MPI_Isend(send_buf[COMMUT_COUNT],1,MPI_DOUBLE,pc_->LEFT(),COMMUT_COUNT,MPI_COMM_WORLD,&request_send[COMMUT_COUNT]);
	  COMMUT_COUNT++;
	}
    }
  else
    {
      MPI_Recv(recv_buf[COMMUT_COUNT],1,MPI_DOUBLE,pc_->RIGHT(),COMMUT_COUNT,MPI_COMM_WORLD,&stat[COMMUT_COUNT]);
      val(Nx_-1,jl,kl)=recv_buf[COMMUT_COUNT][0];
      for (int i(Nx_-bs_-1);i>=bs_;i--)
	val(i,jl,kl)=thomas_d.x(i,jl,kl)-thomas_c.x(i,jl,kl)*val(i+1,jl,kl);	 
      if (!pc_->IsIn())
	{
	  //Send the result to the right process
	  send_buf[COMMUT_COUNT][0]=val(bs_,jl,kl);
	  MPI_Isend(send_buf[COMMUT_COUNT],1,MPI_DOUBLE,pc_->LEFT(),COMMUT_COUNT,MPI_COMM_WORLD,&request_send[COMMUT_COUNT]);
	}
      COMMUT_COUNT++;
    }
}

void poisson_IO::thomas_forward_imag(tensor0 &val,int j,int k)
{
  //compute wavenumbers
  double TWO_PI(2*3.141592653589793);
  double Two_PI_Over_Ly = TWO_PI/size_->Ly();
  double Two_PI_Over_Lz = TWO_PI/size_->Lz();
  jj = ( ( j<Ny_tot_/2 ) ? j : j-Ny_tot_ ) * Two_PI_Over_Ly;
  kk = ( ( k<Nz_tot_/2 ) ? k : k-Nz_tot_ ) * Two_PI_Over_Lz;

  //local j and k indecis
  int jl=j-in_jlo+bs_;
  int kl=k-in_klo+bs_;
  
  double a=inv_dx2;
  double b0=-((dy2)*(1-cos(jj*dy))+(dz2)*(1-cos(kk*dz)))-inv_dx2;
  double b=-((dy2)*(1-cos(jj*dy))+(dz2)*(1-cos(kk*dz)))-2*inv_dx2;
  double c=inv_dx2;

  if (pc_->IsIn()) //first slice
    {
      thomas_c.y(bs_,jl,kl)=c/b0;
      thomas_d.y(bs_,jl,kl)=val(bs_,jl,kl)/b0;
      if (!pc_->IsOut()) //ends in middle
	{
	  for (int i(bs_+1);i<Nx_-bs_;i++)
	    {
	      thomas_c.y(i,jl,kl)=c/(b-thomas_c.y(i-1,jl,kl)*a);
	      thomas_d.y(i,jl,kl)=(val(i,jl,kl)-thomas_d.y(i-1,jl,kl)*a)/(b-thomas_c.y(i-1,jl,kl)*a);
	    }
	  //Send the result to the right process
	  send_buf[COMMUT_COUNT][0]=thomas_c.y(Nx_-bs_-1,jl,kl);
	  send_buf[COMMUT_COUNT][1]=thomas_d.y(Nx_-bs_-1,jl,kl);
	  MPI_Isend(send_buf[COMMUT_COUNT],2,MPI_DOUBLE,pc_->RIGHT(),COMMUT_COUNT,MPI_COMM_WORLD,&request_send[COMMUT_COUNT]);
	  COMMUT_COUNT++;
	}
      else //ends int the last slice
	{
	  for (int i(bs_+1);i<Nx_-bs_-1;i++)
	    {
	      thomas_c.y(i,jl,kl)=c/(b-thomas_c.y(i-1,jl,kl)*a);
	      thomas_d.y(i,jl,kl)=(val(i,jl,kl)-thomas_d.y(i-1,jl,kl)*a)/(b-thomas_c.y(i-1,jl,kl)*a);
	    }
	  if (j+k==0)
	    thomas_d.y(Nx_-bs_-1,jl,kl)=0;
	  else
	    thomas_d.y(Nx_-bs_-1,jl,kl)=(val(Nx_-bs_-1,jl,kl)-thomas_d.y(Nx_-bs_-2,jl,kl)*a)/(b0-thomas_c.y(Nx_-bs_-2,jl,kl)*a);
	}
    }
  else // middle slice
    {
      MPI_Recv(recv_buf[COMMUT_COUNT],2,MPI_DOUBLE,pc_->LEFT(),COMMUT_COUNT,MPI_COMM_WORLD,&stat[COMMUT_COUNT]);
      thomas_c.y(0,jl,kl)=recv_buf[COMMUT_COUNT][0];
      thomas_d.y(0,jl,kl)=recv_buf[COMMUT_COUNT][1];
      if (!pc_->IsOut()) //ends in the middle
	{
	  for (int i(bs_);i<Nx_-bs_;i++)
	    {
	      thomas_c.y(i,jl,kl)=c/(b-thomas_c.y(i-1,jl,kl)*a);
	      thomas_d.y(i,jl,kl)=(val(i,jl,kl)-thomas_d.y(i-1,jl,kl)*a)/(b-thomas_c.y(i-1,jl,kl)*a);
	    }
	  //Send the result to the right process
	  send_buf[COMMUT_COUNT][0]=thomas_c.y(Nx_-bs_-1,jl,kl);
	  send_buf[COMMUT_COUNT][1]=thomas_d.y(Nx_-bs_-1,jl,kl);
	  MPI_Isend(send_buf[COMMUT_COUNT],2,MPI_DOUBLE,pc_->RIGHT(),COMMUT_COUNT,MPI_COMM_WORLD,&request_send[COMMUT_COUNT]);
	}
      else //ends in the last slice
	{
	  for (int i(bs_);i<Nx_-bs_-1;i++)
	    {
	      thomas_c.y(i,jl,kl)=c/(b-thomas_c.y(i-1,jl,kl)*a);
	      thomas_d.y(i,jl,kl)=(val(i,jl,kl)-thomas_d.y(i-1,jl,kl)*a)/(b-thomas_c.y(i-1,jl,kl)*a);
	    }
	  if (j+k==0)
	    thomas_d.y(Nx_-bs_-1,jl,kl)=0;
	  else
	    thomas_d.y(Nx_-bs_-1,jl,kl)=(val(Nx_-bs_-1,jl,kl)-thomas_d.y(Nx_-bs_-2,jl,kl)*a)/(b0-thomas_c.y(Nx_-bs_-2,jl,kl)*a);
	}
      COMMUT_COUNT++;
    }
}

void poisson_IO::thomas_backward_imag(tensor0 &val,int j,int k)
{
  //compute wavenumbers
  double TWO_PI(2*3.141592653589793);
  double Two_PI_Over_Ly = TWO_PI/size_->Ly();
  double Two_PI_Over_Lz = TWO_PI/size_->Lz();
  jj = ( ( j<Ny_tot_/2 ) ? j : j-Ny_tot_ ) * Two_PI_Over_Ly;
  kk = ( ( k<Nz_tot_/2 ) ? k : k-Nz_tot_ ) * Two_PI_Over_Lz;

  //local j and k indecis
  int jl=j-in_jlo+bs_;
  int kl=k-in_klo+bs_;
  
  if (pc_->IsOut()) //the last process does not recv data for forward step
    {
      val(Nx_-bs_-1,jl,kl)=thomas_d.y(Nx_-bs_-1,jl,kl);
      for (int i(Nx_-bs_-2);i>=bs_;i--)
	val(i,jl,kl)=thomas_d.y(i,jl,kl)-thomas_c.y(i,jl,kl)*val(i+1,jl,kl);	 
      if (!pc_->IsIn())
	{
	  //Send the result to the right process
	  send_buf[COMMUT_COUNT][0]=val(bs_,jl,kl);
	  MPI_Isend(send_buf[COMMUT_COUNT],1,MPI_DOUBLE,pc_->LEFT(),COMMUT_COUNT,MPI_COMM_WORLD,&request_send[COMMUT_COUNT]);
	  COMMUT_COUNT++;
	}
    }
  else
    {
      MPI_Recv(recv_buf[COMMUT_COUNT],1,MPI_DOUBLE,pc_->RIGHT(),COMMUT_COUNT,MPI_COMM_WORLD,&stat[COMMUT_COUNT]);
      val(Nx_-1,jl,kl)=recv_buf[COMMUT_COUNT][0];
      for (int i(Nx_-bs_-1);i>=bs_;i--)
	val(i,jl,kl)=thomas_d.y(i,jl,kl)-thomas_c.y(i,jl,kl)*val(i+1,jl,kl);	 
      if (!pc_->IsIn())
	{
	  //Send the result to the right process
	  send_buf[COMMUT_COUNT][0]=val(bs_,jl,kl);
	  MPI_Isend(send_buf[COMMUT_COUNT],1,MPI_DOUBLE,pc_->LEFT(),COMMUT_COUNT,MPI_COMM_WORLD,&request_send[COMMUT_COUNT]);
	}
      COMMUT_COUNT++;
    }
}
