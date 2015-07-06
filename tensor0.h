#ifndef tensor0_h
#define tensor0_h
#include "gridsize.h"
class communicator;
class tensor1;
class tensor0
{
  communicator* comm_;
  gridsize* p_;
  double *ptr; //pointer to the first element, notE: we are using 1D array here and overloading () in order to have fancy 3 indecies.
  double inv_dx,inv_dy,inv_dz;
  double inv_2dx,inv_2dy,inv_2dz;
  double inv_dx2,inv_dy2,inv_dz2;
  int Nx_,Ny_,Nz_,NxNy_;
  double ABS(double);
  public:  
  double *ptr_reduce_; //pointer to a 1D array used for reduce functions (receive buffer) //communicator use this!
  int Nx() const{return p_->Nx();}
  int Ny() const{return p_->Ny();}
  int Nz() const{return p_->Nz();}
  int bs() const{return p_->bs();}
  double Lx() const{return p_->Lx();}
  double Ly() const{return p_->Ly();}
  double Lz() const{return p_->Lz();}
  double dx() const{return p_->dx();}
  double dy() const{return p_->dy();}
  double dz() const{return p_->dz();}
  double *pointer() const{return ptr;}
  tensor0* my_ptr() {return this;}
  int size() const{return p_->size();}
  gridsize* parameter() const {return p_;}
  void Update_Ghosts();
  void Update_Ghosts_CUM(); //using cumulative send_recv (useful for particle projection)
  void Show(); //print out values to the screen
  void ShowIn(); //print out interior values to the screen
  tensor0(gridsize*,communicator*);
  tensor0(){}
  ~tensor0();
  /////////////////////////////////////////operator overloading//////////////////////////
  tensor0& operator=(const double&);
  tensor0& operator+=(const double&);
  tensor0& operator-=(const double&);
  tensor0& operator*=(const double&);
  tensor0& operator/=(const double&);
  tensor0& operator=(const tensor0&);
  tensor0& operator+=(const tensor0&);
  tensor0& operator-=(const tensor0&);
  tensor0& operator*=(const tensor0&);
  tensor0& operator/=(const tensor0&);
  double operator() (int,int,int) const;
  double& operator() (int,int,int);
  /////////////////////////////////////////statistics//////////////////////////////////MPI IMPLEMENTATION/////////
  double sum(); //return sum of all elements in T
  double mean(); //return mean value of T
  void make_mean_zero(); //make mean zero
  void make_reduced_mean_zero(); //make mean of each slice zero
  double max(); //return maximum value
  double max_tot(); //return maximum value including ghost cells
  double max_abs(); //return maximum of absolute value of elements
  double min(); //return minimum value
  double min_tot(); //return minimum value including ghost cells
  double sum_squares(); //return sum of squares of elements
  double rms(); //return root mean squred of elements
  double mean_squares(); //return mean of sum of squares of elements
  double var(); //return variance
  double sd();//return standard deviation
  double corr(tensor0&); //return correlation between this tensor and tensor T
  void reduceX(); //reduce 3D tensor to a 1D tensor in x direction, by averaging all elements
  /////////////////////////////////////POST PROCESS///////////////////////////////////////////////////////////////////////////////
  //Input of dissipation: Nu and SijSij
  void Equal_dissipation(tensor0&,tensor0&); //epsilon=2*nu*SijSij
  //Input TKE: fluctuation of RU and U
  void Equal_TKE(tensor1&,tensor1&); //compute TKE of the given fluctuations velocity field and density
  //Input of uiui: fluctuation of U
  void Equal_UiUi(tensor1&); //calculate UiUi
  //Input of Nu: Re, Rho
  void Equal_nu(double,tensor0&); //nu=1/(Re*Rho)
  //Input of eta: nu and epsilon
  void Equal_eta(tensor0&,tensor0&); //eta=(nu^3/epsilon)^(1/4)
  //Input of ST: Tp, SijSij
  void Equal_ST(double,tensor0&); //ST=sqrt(Tp*Tp*2*SijSij)
  //Input of SijSij: U
  void Equal_SijSij(tensor1&); //compute SijSij of the given fluctuation velocity field
  //Input of WijWij: U
  void Equal_WijWij(tensor1&); //compute WijWij of the given fluctuation velocity field
  //Input tensor T
  void Equal_Reduced_sd(tensor0&); //element=x^2 - E(X)  reduced array=sqrt(E(X^2 - E(X)^2))
  /////////////////////////////////////////arithmatic//////////////////////////////////name format: operator_operator_operator_.. and in prenthesis corresponding data in order
  void PlusEqual_Mult(tensor0&,double); //T=T+A*B
  void PlusEqual_Mult(tensor0&,tensor0&);
  void PlusEqual_Mult(double,tensor0&);
  void Equal_Divide(tensor0&,double); //T=A/B
  void Equal_Divide(tensor0&,tensor0&);
  void Equal_Divide(double,tensor0&);
  void Equal_Mult(tensor0&,double); //T=A*B
  void Equal_Mult(tensor0&,tensor0&);
  void Equal_Mult(double,tensor0&);
  void Equal_LinComb(double,tensor0&,double,tensor0&); //T=A*S+B*R
  void Equal_LinComb(tensor0&,tensor0&,tensor0&,tensor0&);
  void Equal_LinComb(double,tensor0&,tensor0&,tensor0&);
  void Equal_LinComb(tensor0&,tensor0&,double,tensor0&);
  void Equal_LinComb(double,tensor0&,double,tensor0&,double,tensor0&); //T=A*S+B*R+C*T  (specially for lagrangian TWC)
  ///////////////////////////////////////differentiation/////////////////////////////
  double ddx_F2C(int,int,int); //take derivative in x direction from data stored at cell face to cell center
  double ddy_F2C(int,int,int);
  double ddz_F2C(int,int,int);
  double ddx_C2C(int,int,int);
  double ddy_C2C(int,int,int);
  double ddz_C2C(int,int,int);
  double ddx_C2F(int,int,int);
  double ddy_C2F(int,int,int);
  double ddz_C2F(int,int,int);
  double d2dx2(int i,int j,int k); //2nd derivative in x direction C2C and F2F
  double d2dy2(int i,int j,int k);
  double d2dz2(int i,int j,int k);
  void Equal_Div_F2C(tensor1&);
  void Equal_Del2(tensor0&);
  /////////////////////////////////////interpolation/////////////////////////////////
  double Ix_F2C(int,int,int); //interpolate in x direction from faces to center
  double Iy_F2C(int,int,int);
  double Iz_F2C(int,int,int);
  double Ix_C2F(int,int,int);
  double Iy_C2F(int,int,int);
  double Iz_C2F(int,int,int);
  void Equal_Ix_C2F(tensor0&);
  void Equal_Iy_C2F(tensor0&);
  void Equal_Iz_C2F(tensor0&);
  void Equal_Ix_F2C(tensor0&);
  void Equal_Iy_F2C(tensor0&);
  void Equal_Iz_F2C(tensor0&);
  void Equal_I_F2C(tensor1&); //cell center value equals to mean of 6 cell faces
};
#endif
