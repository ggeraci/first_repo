#include "tensor0.h"
#include "tensor1.h"
#include "communicator.h"
#include <mpi.h>
#include <math.h>
///////////////////////////////////////////////////////////////////////////////Class of first order tensor (vector field)////////////////////////////////////////////////////////

void tensor1::Update_Ghosts()
{
  if (p_->IsBox()) comm_->SEND_RECV_BLOCKING(*this); else comm_->SEND_RECV_BLOCKING_INOUT(*this);
}

void tensor1::Update_Ghosts_CUM()
{
  if (p_->IsBox()) comm_->SEND_RECV_BLOCKING_CUM(*this); else comm_->SEND_RECV_BLOCKING_CUM_INOUT(*this);
}

tensor1::tensor1(gridsize* p,communicator *c):x(p,c),y(p,c),z(p,c)
{
  p_=p;
  comm_=c;
}  

tensor1& tensor1::operator=(const double &a)
{
  x=a;
  y=a;
  z=a;
  return *this;
}
tensor1& tensor1::operator+=(const double &a)
{
  x+=a;
  y+=a;
  z+=a;
  return *this;
}
tensor1& tensor1::operator*=(const double &a)
{
  x*=a;
  y*=a;
  z*=a;
  return *this;
}
tensor1& tensor1::operator-=(const double &a)
{
  x-=a;
  y-=a;
  z-=a;
  return *this;
}
tensor1& tensor1::operator/=(const double &a)
{
  x/=a;
  y/=a;
  z/=a;
  return *this;
}

tensor1& tensor1::operator=(const tensor1 &a)
{
  x=a.x;
  y=a.y;
  z=a.z;
  return *this;
}
tensor1& tensor1::operator+=(const tensor1 &a)
{
  x+=a.x;
  y+=a.y;
  z+=a.z;
  return *this;
}
tensor1& tensor1::operator*=(const tensor1 &a)
{
  x*=a.x;
  y*=a.y;
  z*=a.z;
  return *this;
}
tensor1& tensor1::operator-=(const tensor1 &a)
{
  x-=a.x;
  y-=a.y;
  z-=a.z;
  return *this;
}
tensor1& tensor1::operator/=(const tensor1 &a)
{
  x/=a.x;
  y/=a.y;
  z/=a.z;
  return *this;
}

tensor1& tensor1::operator=(const tensor0 &a)
{
  x=a;
  y=a;
  z=a;
  return *this;
}
tensor1& tensor1::operator+=(const tensor0 &a)
{
  x+=a;
  y+=a;
  z+=a;
  return *this;
}
tensor1& tensor1::operator*=(const tensor0 &a)
{
  x*=a;
  y*=a;
  z*=a;
  return *this;
}
tensor1& tensor1::operator-=(const tensor0 &a)
{
  x-=a;
  y-=a;
  z-=a;
  return *this;
}
tensor1& tensor1::operator/=(const tensor0 &a)
{
  x/=a;
  y/=a;
  z/=a;
  return *this;
}
//////////////////////////////////////////////////////////////////////Arithmatic//////////////////////////////////////////////////
void tensor1::PlusEqual_Mult(tensor1 &a,tensor1 &b)
{
  x.PlusEqual_Mult(a.x,b.x);
  y.PlusEqual_Mult(a.y,b.y);
  z.PlusEqual_Mult(a.z,b.z);
}

void tensor1::PlusEqual_Mult(tensor1 &a,double b)
{
  x.PlusEqual_Mult(a.x,b);
  y.PlusEqual_Mult(a.y,b);
  z.PlusEqual_Mult(a.z,b);
}

void tensor1::PlusEqual_Mult(double a,tensor1 &b)
{
  x.PlusEqual_Mult(a,b.x);
  y.PlusEqual_Mult(a,b.y);
  z.PlusEqual_Mult(a,b.z);
}
/////////////
void tensor1::Equal_Divide(tensor1 &a,tensor1 &b)
{
  x.Equal_Divide(a.x,b.x);
  y.Equal_Divide(a.y,b.y);
  z.Equal_Divide(a.z,b.z);
}

void tensor1::Equal_Divide(tensor1 &a,double b)
{
  x.Equal_Divide(a.x,b);
  y.Equal_Divide(a.y,b);
  z.Equal_Divide(a.z,b);
}

void tensor1::Equal_Divide(double a,tensor1 &b)
{
  x.Equal_Divide(a,b.x);
  y.Equal_Divide(a,b.y);
  z.Equal_Divide(a,b.z);
}

void tensor1::Equal_Divide(tensor1 &a,tensor0 &b)
{
  x.Equal_Divide(a.x,b);
  y.Equal_Divide(a.y,b);
  z.Equal_Divide(a.z,b);
}

void tensor1::Equal_Divide(tensor0 &a,tensor1 &b)
{
  x.Equal_Divide(a,b.x);
  y.Equal_Divide(a,b.y);
  z.Equal_Divide(a,b.z);
}
//////////
void tensor1::Equal_Mult(tensor1 &a,tensor1 &b)
{
  x.Equal_Mult(a.x,b.x);
  y.Equal_Mult(a.y,b.y);
  z.Equal_Mult(a.z,b.z);
}

void tensor1::Equal_Mult(tensor1 &a,double b)
{
  x.Equal_Mult(a.x,b);
  y.Equal_Mult(a.y,b);
  z.Equal_Mult(a.z,b);
}

void tensor1::Equal_Mult(double a,tensor1 &b)
{
  x.Equal_Mult(a,b.x);
  y.Equal_Mult(a,b.y);
  z.Equal_Mult(a,b.z);
}
/////////
void tensor1::Equal_LinComb(double a, tensor1 &T,double b, tensor1 &R)
{
  x.Equal_LinComb(a,T.x,b,R.x);
  y.Equal_LinComb(a,T.y,b,R.y);
  z.Equal_LinComb(a,T.z,b,R.z);
}

void tensor1::Equal_LinComb(tensor1 &A, tensor1 &T,double b, tensor1 &R)
{
  x.Equal_LinComb(A.x,T.x,b,R.x);
  y.Equal_LinComb(A.y,T.y,b,R.y);
  z.Equal_LinComb(A.z,T.z,b,R.z);
}
void tensor1::Equal_LinComb(tensor1 &A, tensor1 &T,tensor1 &B, tensor1 &R)
{
  x.Equal_LinComb(A.x,T.x,B.x,R.x);
  y.Equal_LinComb(A.y,T.y,B.y,R.y);
  z.Equal_LinComb(A.z,T.z,B.z,R.z);
}
void tensor1::Equal_LinComb(double a, tensor1 &T,tensor1 &B, tensor1 &R)
{
  x.Equal_LinComb(a,T.x,B.x,R.x);
  y.Equal_LinComb(a,T.y,B.y,R.y);
  z.Equal_LinComb(a,T.z,B.z,R.z);
}

void tensor1::Equal_LinComb(double a, tensor1 &T,double b, tensor1 &R,double c, tensor1 &S)
{
  x.Equal_LinComb(a,T.x,b,R.x,c,S.x);
  y.Equal_LinComb(a,T.y,b,R.y,c,S.y);
  z.Equal_LinComb(a,T.z,b,R.z,c,S.z);
}
///////////////////////////////Differentiation (reuse basics)///////////////////////
void tensor1::Equal_Grad_C2F(tensor0 &T)
{
  for (int i(bs()); i<Nx(); i++)
    for (int j(0); j<Ny(); j++)
      for (int k(0); k<Nz(); k++)
	x(i,j,k)=T.ddx_C2F(i,j,k);

  for (int i(0); i<Nx(); i++)
    for (int j(bs()); j<Ny(); j++)
      for (int k(0); k<Nz(); k++)
	y(i,j,k)=T.ddy_C2F(i,j,k);

  for (int i(0); i<Nx(); i++)
    for (int j(0); j<Ny(); j++)
      for (int k(bs()); k<Nz(); k++)
	z(i,j,k)=T.ddz_C2F(i,j,k);

  Update_Ghosts();
}

void tensor1::Equal_Del2(tensor1 &T)
{
  x.Equal_Del2(T.x);
  y.Equal_Del2(T.y);
  z.Equal_Del2(T.z);
}
////////////////////////////////////Interpolation (reuse basics)//////////////////////////
void tensor1::Equal_I_C2F(tensor0& T)
{
  x.Equal_Ix_C2F(T);
  y.Equal_Iy_C2F(T);
  z.Equal_Iz_C2F(T);
}

void tensor1::Equal_I_F2C(tensor1& T)
{
  x.Equal_Ix_F2C(T.x);
  y.Equal_Iy_F2C(T.y);
  z.Equal_Iz_F2C(T.z);
}

void tensor1::Equal_Ix_C2F(tensor1& T)
{
  x.Equal_Ix_C2F(T.x);
  y.Equal_Ix_C2F(T.y);
  z.Equal_Ix_C2F(T.z);
}

void tensor1::Equal_Iy_C2F(tensor1& T)
{
  x.Equal_Iy_C2F(T.x);
  y.Equal_Iy_C2F(T.y);
  z.Equal_Iy_C2F(T.z);
}

void tensor1::Equal_Iz_C2F(tensor1& T)
{
  x.Equal_Iz_C2F(T.x);
  y.Equal_Iz_C2F(T.y);
  z.Equal_Iz_C2F(T.z);
}
/////////////////////////////////////////////////////////////////////Statistics/////////////////////////////
void tensor1::make_mean_zero()
{
  x.make_mean_zero();
  y.make_mean_zero();
  z.make_mean_zero();
}

void tensor1::make_reduced_mean_zero()
{
  x.make_reduced_mean_zero();
  y.make_reduced_mean_zero();
  z.make_reduced_mean_zero();
}

double tensor1::max()
{
  double maxX=x.max();
  double maxY=y.max();
  double maxZ=z.max();
  maxX=(maxY>maxX)?maxY:maxX;
  maxX=(maxZ>maxX)?maxZ:maxX;
  return maxX;
}

double tensor1::min()
{
  double minX=x.min();
  double minY=y.min();
  double minZ=z.min();
  minX=(minY<minX)?minY:minX;
  minX=(minZ<minX)?minZ:minX;
  return minX;
}
double tensor1::sum()
{
  return x.sum()+y.sum()+z.sum();
}

double tensor1::sum_squares()
{
  return x.sum_squares()+y.sum_squares()+z.sum_squares();
}

double tensor1::mean_squares()
{
  return sum_squares()/(x.parameter()->size_tot());
}

double tensor1::rms()
{
  return sqrt(mean_squares());
}
double tensor1::max_cfl(double dt)
{
  double maxX=x.max()*dt/p_->dx();
  double maxY=y.max()*dt/p_->dy();
  double maxZ=z.max()*dt/p_->dz();
  maxX=(maxY>maxX)?maxY:maxX;
  maxX=(maxZ>maxX)?maxZ:maxX;
  return maxX;
}

void tensor1::reduceX()
{
  comm_->REDUCEX(x);
  comm_->REDUCEX(y);
  comm_->REDUCEX(z);
}


//Post process
void tensor1::Equal_TKE2(tensor1& RU,tensor1& U, tensor0& R) //assume all three are already reduced
{
  for (int i(0); i<Nx();i++)
    for (int j(0);j<Ny();j++)
      for (int k(0);k<Nz();k++)
	{
	  x(i,j,k)=0.5*R(i,j,k)*(U.x(i,j,k)-RU.x.ptr_reduce_[i]/R.ptr_reduce_[i])*(U.x(i,j,k)-RU.x.ptr_reduce_[i]/R.ptr_reduce_[i]);
	  y(i,j,k)=0.5*R(i,j,k)*(U.y(i,j,k)-RU.y.ptr_reduce_[i]/R.ptr_reduce_[i])*(U.y(i,j,k)-RU.y.ptr_reduce_[i]/R.ptr_reduce_[i]);
	  z(i,j,k)=0.5*R(i,j,k)*(U.z(i,j,k)-RU.z.ptr_reduce_[i]/R.ptr_reduce_[i])*(U.z(i,j,k)-RU.z.ptr_reduce_[i]/R.ptr_reduce_[i]);
	}
  reduceX();
}

void tensor1::Equal_UiUi2(tensor1& RU,tensor1& U, tensor0& R) //assume all three are already reduced
{
  for (int i(0); i<Nx();i++)
    for (int j(0);j<Ny();j++)
      for (int k(0);k<Nz();k++)
	{
	  x(i,j,k)=(U.x(i,j,k)-RU.x.ptr_reduce_[i]/R.ptr_reduce_[i])*(U.x(i,j,k)-RU.x.ptr_reduce_[i]/R.ptr_reduce_[i]);
	  y(i,j,k)=(U.y(i,j,k)-RU.y.ptr_reduce_[i]/R.ptr_reduce_[i])*(U.y(i,j,k)-RU.y.ptr_reduce_[i]/R.ptr_reduce_[i]);
	  z(i,j,k)=(U.z(i,j,k)-RU.z.ptr_reduce_[i]/R.ptr_reduce_[i])*(U.z(i,j,k)-RU.z.ptr_reduce_[i]/R.ptr_reduce_[i]);
	}
  reduceX();
}

void tensor1::Equal_Cor_Rho_U(tensor1& U, tensor0& R) //assume all three are already reduced
{
  for (int i(0); i<Nx();i++)
    for (int j(0);j<Ny();j++)
      for (int k(0);k<Nz();k++)
	{
	  x(i,j,k)=(U.x(i,j,k)-U.x.ptr_reduce_[i])*(R(i,j,k)-R.ptr_reduce_[i]);
	  y(i,j,k)=(U.y(i,j,k)-U.y.ptr_reduce_[i])*(R(i,j,k)-R.ptr_reduce_[i]);
	  z(i,j,k)=(U.z(i,j,k)-U.z.ptr_reduce_[i])*(R(i,j,k)-R.ptr_reduce_[i]);
	}
  reduceX();
}


void tensor1::Equal_TKE(tensor1& U, tensor0& R) //assume all three are already reduced
{
  for (int i(0); i<Nx();i++)
    for (int j(0);j<Ny();j++)
      for (int k(0);k<Nz();k++)
	{
	  x(i,j,k)=0.5*R(i,j,k)*(U.x(i,j,k)-U.x.ptr_reduce_[i])*(U.x(i,j,k)-U.x.ptr_reduce_[i]);
	  y(i,j,k)=0.5*R(i,j,k)*(U.y(i,j,k)-U.y.ptr_reduce_[i])*(U.y(i,j,k)-U.y.ptr_reduce_[i]);
	  z(i,j,k)=0.5*R(i,j,k)*(U.z(i,j,k)-U.z.ptr_reduce_[i])*(U.z(i,j,k)-U.z.ptr_reduce_[i]);
	}
  reduceX();
}

void tensor1::Equal_UiUi(tensor1& U, tensor0& R) //assume all three are already reduced
{
  for (int i(0); i<Nx();i++)
    for (int j(0);j<Ny();j++)
      for (int k(0);k<Nz();k++)
	{
	  x(i,j,k)=(U.x(i,j,k)-U.x.ptr_reduce_[i])*(U.x(i,j,k)-U.x.ptr_reduce_[i]);
	  y(i,j,k)=(U.y(i,j,k)-U.y.ptr_reduce_[i])*(U.y(i,j,k)-U.y.ptr_reduce_[i]);
	  z(i,j,k)=(U.z(i,j,k)-U.z.ptr_reduce_[i])*(U.z(i,j,k)-U.z.ptr_reduce_[i]);
	}
  reduceX();
}
