#include "tensor0.h"
#include "tensor1.h"
#include "communicator.h"
#include <mpi.h>
#include <math.h>
#include <iostream>
double tensor0::ABS(double a)
{
  return (a>=0)?a:-a;
}
///////////////////////////////////////////////////////////////////////////////Class of first order tensor (vector field)////////////////////////////////////////////////////////
tensor0::tensor0(gridsize* p,communicator *c)
{
  p_=p;
  comm_=c;
  ptr=new double[p_->size()];
  ptr_reduce_=new double[p_->Nx()];
  inv_dx=1./dx();
  inv_dy=1./dy();
  inv_dz=1./dz();
  inv_2dx=1./(2*dx());
  inv_2dy=1./(2*dy());
  inv_2dz=1./(2*dz());
  inv_dx2=1./(dx()*dx());
  inv_dy2=1./(dy()*dy());
  inv_dz2=1./(dz()*dz());
  Nx_=Nx();
  Ny_=Ny();
  Nz_=Nz();
  NxNy_=Nx_*Ny_;
}

tensor0::~tensor0()
{
  delete[] ptr;
}

void tensor0::Show()
{
  for (int i(0);i<Nx_;i++)
    {
    for (int j(0);j<Ny_;j++)
      {
	std::cout<<i<<","<<j<<std::endl;
	for (int k(0);k<Nz_;k++)
	  {
	    std::cout<<(*this)(i,j,k)<<"  ";
	  }
	std::cout<<std::endl;
      }
    std::cout<<std::endl<<std::endl;
    }
}

void tensor0::ShowIn()
{
  for (int i(bs());i<Nx_-bs();i++)
    {
      for (int j(bs());j<Ny_-bs();j++)
      {
	std::cout<<i<<","<<j<<std::endl;
	for (int k(bs());k<Nz_-bs();k++)
	  {
	    std::cout<<(*this)(i,j,k)<<"  ";
	  }
	std::cout<<std::endl;
      }
    std::cout<<std::endl<<std::endl;
    }
}

double tensor0::operator() (int i,int j,int k) const
{
  return ptr[k*NxNy_+j*Nx_+i];
}

double& tensor0::operator() (int i,int j,int k)
{
  return ptr[k*NxNy_+j*Nx_+i];
}


void tensor0::Update_Ghosts()
{
  if (p_->IsBox()) comm_->SEND_RECV_BLOCKING(*this); else comm_->SEND_RECV_BLOCKING_INOUT(*this);
}

void tensor0::Update_Ghosts_CUM()
{
 if (p_->IsBox()) comm_->SEND_RECV_BLOCKING_CUM(*this); else comm_->SEND_RECV_BLOCKING_CUM_INOUT(*this);
}

tensor0& tensor0::operator=(const double &a)
{
  for (int i(0);i<p_->size();i++)
    ptr[i]=a;
  return *this;
}
tensor0& tensor0::operator+=(const double &a)
{
  for (int i(0);i<p_->size();i++)
    ptr[i]+=a;
  return *this;
}
tensor0& tensor0::operator*=(const double &a)
{
  for (int i(0);i<p_->size();i++)
    ptr[i]*=a;
  return *this;
}
tensor0& tensor0::operator-=(const double &a)
{
  for (int i(0);i<p_->size();i++)
    ptr[i]-=a;
  return *this;
}
tensor0& tensor0::operator/=(const double &a)
{
  for (int i(0);i<p_->size();i++)
    ptr[i]/=a;
  return *this;
}

tensor0& tensor0::operator=(const tensor0 &a)
{
  for (int i(0);i<p_->size();i++)
    ptr[i]=a.pointer()[i];
  return *this;
}
tensor0& tensor0::operator+=(const tensor0 &a)
{
  for (int i(0);i<p_->size();i++)
    ptr[i]+=a.pointer()[i];
  return *this;
}
tensor0& tensor0::operator*=(const tensor0 &a)
{
  for (int i(0);i<p_->size();i++)
    ptr[i]*=a.pointer()[i];
  return *this;
}
tensor0& tensor0::operator-=(const tensor0 &a)
{
  for (int i(0);i<p_->size();i++)
    ptr[i]-=a.pointer()[i];
  return *this;
}
tensor0& tensor0::operator/=(const tensor0 &a)
{
  for (int i(0);i<p_->size();i++)
    ptr[i]/=a.pointer()[i];
  return *this;
}

//////////////////////////////////////////////////////////////////////Arithmatic//////////////////////////////////////////////////
void tensor0::PlusEqual_Mult(tensor0 &a,tensor0 &b) //T+=a*b
{
  for (int i(0);i<p_->size();i++)
    ptr[i]+=a.pointer()[i]*b.pointer()[i];
}

void tensor0::PlusEqual_Mult(tensor0 &a,double b)
{
  for (int i(0);i<p_->size();i++)
    ptr[i]+=a.pointer()[i]*b;
}

void tensor0::PlusEqual_Mult(double a,tensor0 &b)
{
  for (int i(0);i<p_->size();i++)
    ptr[i]+=a*b.pointer()[i];
}

/////////////

void tensor0::Equal_Divide(tensor0 &a,tensor0 &b) //T=a/b
{
  for (int i(0);i<p_->size();i++)
    ptr[i]=a.pointer()[i]/b.pointer()[i];
}

void tensor0::Equal_Divide(tensor0 &a,double b)
{
  for (int i(0);i<p_->size();i++)
    ptr[i]=a.pointer()[i]/b;
}

void tensor0::Equal_Divide(double a,tensor0 &b)
{
  for (int i(0);i<p_->size();i++)
    ptr[i]=a/b.pointer()[i];
}
//////////
void tensor0::Equal_Mult(tensor0 &a,tensor0 &b) //T=a*b
{
  for (int i(0);i<p_->size();i++)
    ptr[i]=a.pointer()[i]*b.pointer()[i];
}

void tensor0::Equal_Mult(tensor0 &a,double b)
{
  for (int i(0);i<p_->size();i++)
    ptr[i]=a.pointer()[i]*b;
}

void tensor0::Equal_Mult(double a,tensor0 &b)
{
  for (int i(0);i<p_->size();i++)
    ptr[i]=a*b.pointer()[i];
}
/////////

void tensor0::Equal_LinComb(double a, tensor0 &T,double b, tensor0 &R)
{
  for (int i(0);i<p_->size();i++)
    ptr[i]=a*T.pointer()[i]+b*R.pointer()[i];
}

void tensor0::Equal_LinComb(tensor0 &A, tensor0 &T,double b, tensor0 &R)
{
  for (int i(0);i<p_->size();i++)
    ptr[i]=A.pointer()[i]*T.pointer()[i]+b*R.pointer()[i];
}

void tensor0::Equal_LinComb(tensor0& A, tensor0 &T,tensor0& B, tensor0 &R)
{
  for (int i(0);i<p_->size();i++)
    ptr[i]=A.pointer()[i]*T.pointer()[i]+B.pointer()[i]*R.pointer()[i];
}

void tensor0::Equal_LinComb(double a, tensor0 &T,tensor0& B, tensor0 &R)
{
  for (int i(0);i<p_->size();i++)
    ptr[i]=a*T.pointer()[i]+B.pointer()[i]*R.pointer()[i];
}

void tensor0::Equal_LinComb(double a, tensor0 &T,double b, tensor0 &R,double c,tensor0 &S)
{
  for (int i(0);i<p_->size();i++)
    ptr[i]=a*T.pointer()[i]+b*R.pointer()[i]+c*S.pointer()[i];
}
////////////////////////////////////////////////////////////////Differentiation(Basics)/////////////////////////////////////////////////////////////////
double tensor0::ddx_F2C(int i,int j,int k)
{
  return ((*this)(i+1,j,k)-(*this)(i,j,k))*inv_dx;
}

double tensor0::ddy_F2C(int i,int j,int k)
{
  return ((*this)(i,j+1,k)-(*this)(i,j,k))*inv_dy;
}

double tensor0::ddz_F2C(int i,int j,int k)
{
  return ((*this)(i,j,k+1)-(*this)(i,j,k))*inv_dz;
}

double tensor0::ddx_C2C(int i,int j,int k)
{
  return ((*this)(i+1,j,k)-(*this)(i-1,j,k))*inv_2dx;
}

double tensor0::ddy_C2C(int i,int j,int k)
{
  return ((*this)(i,j+1,k)-(*this)(i,j-1,k))*inv_2dy;
}

double tensor0::ddz_C2C(int i,int j,int k)
{
  return ((*this)(i,j,k+1)-(*this)(i,j,k-1))*inv_2dz;
}

double tensor0::ddx_C2F(int i,int j,int k)
{
  return ((*this)(i,j,k)-(*this)(i-1,j,k))*inv_dx;
}

double tensor0::ddy_C2F(int i,int j,int k)
{
  return ((*this)(i,j,k)-(*this)(i,j-1,k))*inv_dy;
}

double tensor0::ddz_C2F(int i,int j,int k)
{
  return ((*this)(i,j,k)-(*this)(i,j,k-1))*inv_dz;
}
//second derivative is infact: d2dx2= ddx_F2C(ddx_C2F) = ddx_C2F(ddx_F2C)
double tensor0::d2dx2(int i,int j,int k)
{
  return ((*this)(i+1,j,k)-2*(*this)(i,j,k)+(*this)(i-1,j,k))*inv_dx2;
}

double tensor0::d2dy2(int i,int j,int k)
{
  return ((*this)(i,j+1,k)-2*(*this)(i,j,k)+(*this)(i,j-1,k))*inv_dy2;
}

double tensor0::d2dz2(int i,int j,int k)
{
  return ((*this)(i,j,k+1)-2*(*this)(i,j,k)+(*this)(i,j,k-1))*inv_dz2;
}
///////////////////////////////Differentiation (reuse basics)///////////////////////
void tensor0::Equal_Div_F2C(tensor1 &T)
{
  for (int i(0); i<Nx()-bs(); i++)
    for (int j(0); j<Ny()-bs(); j++)
      for (int k(0); k<Nz()-bs(); k++)
	(*this)(i,j,k)=T.x.ddx_F2C(i,j,k)+T.y.ddy_F2C(i,j,k)+T.z.ddz_F2C(i,j,k);
  Update_Ghosts();
}

void tensor0::Equal_Del2(tensor0 &T)
{
  for (int i(bs()); i<Nx()-bs(); i++)
    for (int j(bs()); j<Ny()-bs(); j++)
      for (int k(bs()); k<Nz()-bs(); k++)
	(*this)(i,j,k)=T.d2dx2(i,j,k)+T.d2dy2(i,j,k)+T.d2dz2(i,j,k);
  Update_Ghosts();
}

////////////////////////////////////////////////////Interpolation (Basics)//////////////////////////////////////////////////////////////////////
double tensor0::Ix_F2C(int i,int j,int k)
{
  return ((*this)(i,j,k)+(*this)(i+1,j,k))*0.5;
}

double tensor0::Iy_F2C(int i,int j,int k)
{
  return ((*this)(i,j,k)+(*this)(i,j+1,k))*0.5;
}

double tensor0::Iz_F2C(int i,int j,int k)
{
  return ((*this)(i,j,k)+(*this)(i,j,k+1))*0.5;
}

double tensor0::Ix_C2F(int i,int j,int k)
{
  return ((*this)(i,j,k)+(*this)(i-1,j,k))*0.5;
}

double tensor0::Iy_C2F(int i,int j,int k)
{
  return ((*this)(i,j,k)+(*this)(i,j-1,k))*0.5;
}

double tensor0::Iz_C2F(int i,int j,int k)
{
  return ((*this)(i,j,k)+(*this)(i,j,k-1))*0.5;
}

////////////////////////////////////Interpolation (reuse basics)//////////////////////////
void tensor0::Equal_Ix_C2F(tensor0& T)
{
  for (int i(bs()); i<Nx(); i++)
    for (int j(0); j<Ny(); j++)
      for (int k(0); k<Nz(); k++)
	(*this)(i,j,k)=T.Ix_C2F(i,j,k);
  Update_Ghosts();
}

void tensor0::Equal_Iy_C2F(tensor0& T)
{
  for (int i(0); i<Nx(); i++)
    for (int j(bs()); j<Ny(); j++)
      for (int k(0); k<Nz(); k++)
	(*this)(i,j,k)=T.Iy_C2F(i,j,k);
  Update_Ghosts();
}

void tensor0::Equal_Iz_C2F(tensor0& T)
{
  for (int i(0); i<Nx(); i++)
    for (int j(0); j<Ny(); j++)
      for (int k(bs()); k<Nz(); k++)
	(*this)(i,j,k)=T.Iz_C2F(i,j,k);
  Update_Ghosts();
}


void tensor0::Equal_Ix_F2C(tensor0& T)
{
  for (int i(0); i<Nx()-bs(); i++)
    for (int j(0); j<Ny(); j++)
      for (int k(0); k<Nz(); k++)
	(*this)(i,j,k)=T.Ix_F2C(i,j,k);
  Update_Ghosts();
}

void tensor0::Equal_Iy_F2C(tensor0& T)
{
  for (int i(0); i<Nx(); i++)
    for (int j(0); j<Ny()-bs(); j++)
      for (int k(0); k<Nz(); k++)
	(*this)(i,j,k)=T.Iy_F2C(i,j,k);
  Update_Ghosts();
}

void tensor0::Equal_Iz_F2C(tensor0& T)
{
  for (int i(0); i<Nx(); i++)
    for (int j(0); j<Ny(); j++)
      for (int k(0); k<Nz()-bs(); k++)
	(*this)(i,j,k)=T.Iz_F2C(i,j,k);
  Update_Ghosts();
}

void tensor0::Equal_I_F2C(tensor1& T)
{
  double OneThird=1./3.;
  for (int i(0); i<Nx()-bs(); i++)
    for (int j(0); j<Ny()-bs(); j++)
      for (int k(0); k<Nz()-bs(); k++)
	(*this)(i,j,k)=OneThird*(T.x.Ix_F2C(i,j,k)+T.y.Iy_F2C(i,j,k)+T.z.Iz_F2C(i,j,k));
  Update_Ghosts();
}
/////////////////////////////////////////////////////////////////////Statistics/////////////////////////////
double tensor0::sum()
{
  double LocalSum(0);
  for (int i(bs()); i<Nx()-bs(); i++)
    for (int j(bs()); j<Ny()-bs(); j++)
      for (int k(bs()); k<Nz()-bs(); k++)
	LocalSum+=(*this)(i,j,k);
  double GlobalSum;
  MPI_Reduce(&LocalSum,&GlobalSum,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
  MPI_Bcast(&GlobalSum,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
  return GlobalSum; 
}

double tensor0::mean()
{
  double GlobalSum=sum();
  return GlobalSum/p_->size_tot();
}

void tensor0::make_mean_zero()
{
  double Mean=mean();
  (*this)-=Mean;
}

void tensor0::make_reduced_mean_zero()
{
  reduceX();
  for (int i(0); i<Nx();i++)
    for (int j(0);j<Ny();j++)
      for (int k(0);k<Nz();k++)
	(*this)(i,j,k)-=ptr_reduce_[i];
}

double tensor0::max()
{
  double LocalMax=(*this)(bs(),bs(),bs());
  for (int i(bs()); i<Nx()-bs();i++)
    for (int j(bs());j<Ny()-bs();j++)
      for (int k(bs());k<Nz()-bs();k++)
	if ((*this)(i,j,k)>LocalMax) LocalMax=(*this)(i,j,k);
  double GlobalMax;
  MPI_Reduce(&LocalMax,&GlobalMax,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);
  MPI_Bcast(&GlobalMax,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
  return GlobalMax;
}

double tensor0::max_tot()
{
  double LocalMax=(*this)(bs(),bs(),bs());
  for (int i(0); i<Nx();i++)
    for (int j(0);j<Ny();j++)
      for (int k(0);k<Nz();k++)
	if ((*this)(i,j,k)>LocalMax) LocalMax=(*this)(i,j,k);
  double GlobalMax;
  MPI_Reduce(&LocalMax,&GlobalMax,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);
  MPI_Bcast(&GlobalMax,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
  return GlobalMax;
}

double tensor0::max_abs()
{
  double LocalMax=ABS((*this)(bs(),bs(),bs()));
  for (int i(bs()); i<Nx()-bs();i++)
    for (int j(bs());j<Ny()-bs();j++)
      for (int k(bs());k<Nz()-bs();k++)
	if (ABS((*this)(i,j,k))>LocalMax) LocalMax=ABS((*this)(i,j,k));
  double GlobalMax;
  MPI_Reduce(&LocalMax,&GlobalMax,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);
  MPI_Bcast(&GlobalMax,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
  return GlobalMax;
}


double tensor0::min()
{
  double LocalMin=(*this)(bs(),bs(),bs());
  for (int i(bs()); i<Nx()-bs();i++)
    for (int j(bs());j<Ny()-bs();j++)
      for (int k(bs());k<Nz()-bs();k++)
	if ((*this)(i,j,k)<LocalMin) LocalMin=(*this)(i,j,k);
  double GlobalMin;
  MPI_Reduce(&LocalMin,&GlobalMin,1,MPI_DOUBLE,MPI_MIN,0,MPI_COMM_WORLD);
  MPI_Bcast(&GlobalMin,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
  return GlobalMin;
}

double tensor0::min_tot()
{
  double LocalMin=(*this)(bs(),bs(),bs());
  for (int i(0); i<Nx();i++)
    for (int j(0);j<Ny();j++)
      for (int k(0);k<Nz();k++)
	if ((*this)(i,j,k)<LocalMin) LocalMin=(*this)(i,j,k);
  double GlobalMin;
  MPI_Reduce(&LocalMin,&GlobalMin,1,MPI_DOUBLE,MPI_MIN,0,MPI_COMM_WORLD);
  MPI_Bcast(&GlobalMin,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
  return GlobalMin;
}

double tensor0::sum_squares()
{
  double LocalSumSquares(0);
  for (int i(bs()); i<Nx()-bs(); i++)
    for (int j(bs()); j<Ny()-bs(); j++)
      for (int k(bs()); k<Nz()-bs(); k++)
	LocalSumSquares+=(*this)(i,j,k)*(*this)(i,j,k);
  double GlobalSumSquares;
  MPI_Reduce(&LocalSumSquares,&GlobalSumSquares,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
  MPI_Bcast(&GlobalSumSquares,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
  return GlobalSumSquares; 
}

double tensor0::mean_squares()
{
  return sum_squares()/p_->size_tot();
}

double tensor0::rms()
{
  return sqrt(mean_squares());
}

double tensor0::var()
{
  return mean_squares()-mean()*mean();
}

double tensor0::sd()
{
  return sqrt(var());
}

double tensor0::corr(tensor0& T)
{
  double mu1=mean();
  double mu2=T.mean();
  double sd1=sd();
  double sd2=sd();
  double LocalCov=0;
  for (int i(bs()); i<Nx()-bs(); i++)
    for (int j(bs()); j<Ny()-bs(); j++)
      for (int k(bs()); k<Nz()-bs(); k++)
	LocalCov+=((*this)(i,j,k)-mu1)*(T(i,j,k)-mu2);
  double GlobalCov;
  MPI_Reduce(&LocalCov,&GlobalCov,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
  MPI_Bcast(&GlobalCov,1,MPI_DOUBLE,0,MPI_COMM_WORLD);  
  return GlobalCov/(p_->size_tot()*sd1*sd2);
}

void tensor0::reduceX()
{
  comm_->REDUCEX(*this);
}


///////////////////////////////////////POST PROCESS///////////////////////////////////
void tensor0::Equal_UiUi(tensor1& U)
{
  for (int i(0); i<Nx();i++)
    for (int j(0);j<Ny();j++)
      for (int k(0);k<Nz();k++)
	(*this)(i,j,k)=(U.x(i,j,k)*U.x(i,j,k)+U.y(i,j,k)*U.y(i,j,k)+U.z(i,j,k)*U.z(i,j,k));
  reduceX();
}

void tensor0::Equal_SijSij(tensor1& U)
{
  for (int i=bs();i<Nx()-bs();i++)
    for (int j=bs();j<Ny()-bs();j++)
      for (int k=bs();k<Nz()-bs();k++) 
	{
	  (*this)(i,j,k)=pow(U.x.ddx_F2C(i,j,k),2)+pow(U.y.ddy_F2C(i,j,k),2)+pow(U.z.ddz_F2C(i,j,k),2);
	  (*this)(i,j,k)+=0.125*pow(U.x.ddy_C2C(i,j,k)+U.x.ddy_C2C(i+1,j,k)+U.y.ddx_C2C(i,j,k)+U.y.ddx_C2C(i,j+1,k),2);
	  (*this)(i,j,k)+=0.125*pow(U.x.ddz_C2C(i,j,k)+U.x.ddz_C2C(i+1,j,k)+U.z.ddx_C2C(i,j,k)+U.z.ddx_C2C(i,j,k+1),2);
	  (*this)(i,j,k)+=0.125*pow(U.y.ddz_C2C(i,j,k)+U.y.ddz_C2C(i,j+1,k)+U.z.ddy_C2C(i,j,k)+U.z.ddy_C2C(i,j,k+1),2);
	}
}

void tensor0::Equal_WijWij(tensor1& U)
{
  for (int i=bs();i<Nx()-bs();i++)
    for (int j=bs();j<Ny()-bs();j++)
      for (int k=bs();k<Nz()-bs();k++) 
	{
	  (*this)(i,j,k)+=0.125*pow(U.x.ddy_C2C(i,j,k)+U.x.ddy_C2C(i+1,j,k)-U.y.ddx_C2C(i,j,k)-U.y.ddx_C2C(i,j+1,k),2);
	  (*this)(i,j,k)+=0.125*pow(U.x.ddz_C2C(i,j,k)+U.x.ddz_C2C(i+1,j,k)-U.z.ddx_C2C(i,j,k)-U.z.ddx_C2C(i,j,k+1),2);
	  (*this)(i,j,k)+=0.125*pow(U.y.ddz_C2C(i,j,k)+U.y.ddz_C2C(i,j+1,k)-U.z.ddy_C2C(i,j,k)-U.z.ddy_C2C(i,j,k+1),2);
	}
}

void tensor0::Equal_TKE(tensor1& U,tensor1& RU)
{
  for (int i(0); i<Nx();i++)
    for (int j(0);j<Ny();j++)
      for (int k(0);k<Nz();k++)
	(*this)(i,j,k)=.5*(U.x(i,j,k)*RU.x(i,j,k)+U.y(i,j,k)*RU.y(i,j,k)+U.z(i,j,k)*RU.z(i,j,k));
  reduceX();
}

void tensor0::Equal_dissipation(tensor0& nu,tensor0& sijsij)
{
   for (int i(0); i<Nx();i++)
    for (int j(0);j<Ny();j++)
      for (int k(0);k<Nz();k++)
	(*this)(i,j,k)=2*nu(i,j,k)*sijsij(i,j,k);
  reduceX();
}

void tensor0::Equal_nu(double mu,tensor0& Rho)
{
  Equal_Divide(mu,Rho);
  reduceX();
}

void tensor0::Equal_eta(tensor0 &nu,tensor0 &epsilon)
{
  for (int i(0); i<Nx();i++)
    for (int j(0);j<Ny();j++)
      for (int k(0);k<Nz();k++)
	(*this)(i,j,k)=pow(nu(i,j,k)*nu(i,j,k)*nu(i,j,k)/epsilon(i,j,k),0.25);
  reduceX();
}

void tensor0::Equal_ST(double Tp, tensor0 &sijsij)
{
  double Tp2=2*Tp*Tp;
  for (int i(0); i<Nx();i++)
    for (int j(0);j<Ny();j++)
      for (int k(0);k<Nz();k++)
	(*this)(i,j,k)=sqrt(Tp2*sijsij(i,j,k));
  reduceX();
}

void tensor0::Equal_Reduced_sd(tensor0& T)
{
  T.reduceX();
  for (int i(0); i<Nx();i++)
    for (int j(0);j<Ny();j++)
      for (int k(0);k<Nz();k++)
	(*this)(i,j,k)=(T(i,j,k)-T.ptr_reduce_[i])*(T(i,j,k)-T.ptr_reduce_[i]);
  reduceX();
  for (int i(0); i<Nx();i++)
    ptr_reduce_[i]=sqrt(ptr_reduce_[i]);
}
