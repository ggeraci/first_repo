#include "proc.h"
#include <math.h>
proc::proc()
{
  MPI_Comm_rank(MPI_COMM_WORLD, &myRank_);
  MPI_Comm_size(MPI_COMM_WORLD, &totalNumProcessors_);
  NumProcCalculator();
  proc_i_=myRank_%NumProcX_;
  proc_j_=(myRank_/NumProcX_)%NumProcY_;
  proc_k_=myRank_/(NumProcX_*NumProcY_);

  procRIGHT_=((proc_i_+1)%NumProcX_)+proc_j_*NumProcX_+proc_k_*NumProcX_*NumProcY_;
  procLEFT_=((proc_i_-1+NumProcX_)%NumProcX_)+proc_j_*NumProcX_+proc_k_*NumProcX_*NumProcY_;
  procTOP_=proc_i_+((proc_j_+1)%NumProcY_)*NumProcX_+proc_k_*NumProcX_*NumProcY_;
  procBOT_=proc_i_+((proc_j_-1+NumProcY_)%NumProcY_)*NumProcX_+proc_k_*NumProcX_*NumProcY_;
  procFRONT_=proc_i_+proc_j_*NumProcX_+((proc_k_+1)%NumProcZ_)*NumProcX_*NumProcY_;
  procREAR_=proc_i_+proc_j_*NumProcX_+((proc_k_-1+NumProcZ_)%NumProcZ_)*NumProcX_*NumProcY_;
  procROOT_=0;

  //update list of processes with same i in he logical grid
  XMates_ranks=new int[NumProcY_*NumProcZ_];
  //update list of processes with same j&k in he logical grid
  YZMates_ranks=new int[NumProcX_];
  int count=0;
  for (int i(0);i<totalNumProcessors_;i++)
    if ((i%NumProcX_)==proc_i_) XMates_ranks[count++]=i;
  count=0;
  for (int i(0);i<totalNumProcessors_;i++)
    if (((i/NumProcX_)%NumProcY_==proc_j_)&&(i/(NumProcX_*NumProcY_)==proc_k_)) YZMates_ranks[count++]=i;
  //define a group process with same i in the logical grid
  MPI_Comm_group(MPI_COMM_WORLD,&all_group_);
  MPI_Group_incl(all_group_,NumProcY_*NumProcZ_,XMates_ranks,&XMates_group_);
  MPI_Comm_create(MPI_COMM_WORLD,XMates_group_,&XMates_comm_);
  MPI_Group_incl(all_group_,NumProcX_,YZMates_ranks,&YZMates_group_);
  MPI_Comm_create(MPI_COMM_WORLD,YZMates_group_,&YZMates_comm_);

  MPI_Comm_rank(YZMates_comm_, &myRank_YZ_);
  MPI_Comm_rank(XMates_comm_, &myRank_X_);
}


void proc::NumProcCalculator()
{

  int n=totalNumProcessors_;
  int m=n;
  if (m%4==0)
    {
      m/=4;
      NumProcX_ = 4;
    }
  else
    {
      NumProcX_ = 1;
    }
  int a=int(pow(m,1./3.))+1;

  for (int i(a);i>0;i--)
    if (m%i==0)
      {
        NumProcX_ *= i;
        break;
      }
  //if (n%(NumProcX_*4)==0) NumProcX_*=4;                                                                   
  n = n / NumProcX_;
  a = int(sqrt(n)) + 1;
  for (int i(a);i>0;i--)
    if (n%i == 0)
      {
        NumProcY_ = i;
        break;
      }
  NumProcZ_ = n / NumProcY_;
}

proc::~proc()
{
  //delete[] XMates_ranks;
  //delete[] YZMates_ranks;
  /*
  MPI_Comm_free(&XMates_comm_);
  MPI_Comm_free(&YZMates_comm_);
  MPI_Group_free(&all_group_);
  MPI_Group_free(&XMates_group_);
  MPI_Group_free(&YZMates_group_);
  */
}


