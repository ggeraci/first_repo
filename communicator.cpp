#include <fstream>
#include <sstream>
#include <math.h>
#include "proc.h"
#include "params.h"
#include "gridsize.h"
#include "tensor0.h"
#include "tensor1.h"
#include "communicator.h"
communicator::communicator(gridsize* s,params* p,proc *pc)
{
  size_=s;
  param_=p;
  pc_=pc;
  Nx=s->Nx();
  Ny=s->Ny();
  Nz=s->Nz();
  bs=s->bs();
  proc_i=pc->I();
  proc_j=pc->J();
  proc_k=pc->K();
  procRIGHT=pc->RIGHT();
  procLEFT=pc->LEFT();
  procTOP=pc->TOP();
  procBOT=pc->BOT();
  procFRONT=pc->FRONT();
  procREAR=pc->REAR();
  //////////////////////////////data storage initializing (full tensor)////////////////////////////////////
  TNOD=s->Nx_tot()*s->Ny_tot()*s->Nz_tot()*sizeof(double);
  gsizes[0]=s->Nz_tot();  gsizes[1]=s->Ny_tot();  gsizes[2]=s->Nx_tot();
  psizes[0]=pc->NZ();  psizes[1]=pc->NY();  psizes[2]=pc->NX();
  lsizes[0]=s->Nz()-2*s->bs();  lsizes[1]=s->Ny()-2*s->bs();  lsizes[2]=s->Nx()-2*s->bs();
  start_indices_g[0]=s->kl();  start_indices_g[1]=s->jl();  start_indices_g[2]=s->il();
  memsizes[0]=s->Nz();  memsizes[1]=s->Ny();  memsizes[2]=s->Nx();
  start_indices_l[0]=s->bs();  start_indices_l[1]=s->bs();  start_indices_l[2]=s->bs();
  MPI_Type_create_subarray(3,gsizes,lsizes,start_indices_g,MPI_ORDER_C,MPI_DOUBLE,&filetype);
  MPI_Type_commit(&filetype);
  MPI_Type_create_subarray(3,memsizes,lsizes,start_indices_l,MPI_ORDER_C,MPI_DOUBLE,&memtype);
  MPI_Type_commit(&memtype);
  //////////////////////////////data storage initializing (one slice)////////////////////////////////////
  TNOD_S=s->Ny_tot()*s->Nz_tot()*sizeof(double);
  gsizes_S[0]=s->Nz_tot();  gsizes_S[1]=s->Ny_tot();  gsizes_S[2]=1;
  psizes_S[0]=pc->NZ();  psizes_S[1]=pc->NY();  psizes_S[2]=1;
  lsizes_S[0]=s->Nz()-2*s->bs();  lsizes_S[1]=s->Ny()-2*s->bs();  lsizes_S[2]=1;
  start_indices_g_S[0]=s->kl();  start_indices_g_S[1]=s->jl();  start_indices_g_S[2]=0;
  memsizes_S[0]=s->Nz();  memsizes_S[1]=s->Ny();  memsizes_S[2]=s->Nx();
  start_indices_l_S[0]=s->bs();  start_indices_l_S[1]=s->bs();  start_indices_l_S[2]=s->Nx()-1;
  MPI_Type_create_subarray(3,gsizes_S,lsizes_S,start_indices_g_S,MPI_ORDER_C,MPI_DOUBLE,&filetype_S);
  MPI_Type_commit(&filetype_S);
  MPI_Type_create_subarray(3,memsizes_S,lsizes_S,start_indices_l_S,MPI_ORDER_C,MPI_DOUBLE,&memtype_S);
  MPI_Type_commit(&memtype_S);
  /////////////////////////////communication buffer data storage//////////////////
  N_max=((s->Nx()>s->Ny())?s->Nx():s->Ny());
  N_max=((N_max>s->Nz())?N_max:s->Nz());
  N_max++;
  Sbuf=new double*[2];
  Rbuf=new double*[2];
  Sbuf[0]=new double[N_max*N_max*s->bs()];
  Sbuf[1]=new double[N_max*N_max*s->bs()];
  Rbuf[0]=new double[N_max*N_max*s->bs()];
  Rbuf[1]=new double[N_max*N_max*s->bs()];
}

communicator::~communicator()
{
  //MPI_Type_free(&filetype);
  //MPI_Type_free(&memtype);
  for (int i(0);i<2;i++)
    {
      delete[] Sbuf[i];
      delete[] Rbuf[i];
    }
  delete [] Rbuf;
  delete [] Sbuf;
}

void communicator::read_S(tensor0 &T,char *adrs)
{
  MPI_File_open(pc_->XMates(),adrs,MPI_MODE_RDONLY,MPI_INFO_NULL,&Myfile);
  MPI_File_set_view(Myfile,0,MPI_DOUBLE,filetype_S,"native",MPI_INFO_NULL);
  MPI_File_read_all(Myfile,T.pointer(),1,memtype_S,&status);
  MPI_File_close(&Myfile); 
}


void communicator::read(tensor0 &T,char *adrs)
{
  MPI_File_open(MPI_COMM_WORLD,adrs,MPI_MODE_RDONLY,MPI_INFO_NULL,&Myfile);
  MPI_File_set_view(Myfile,0,MPI_DOUBLE,filetype,"native",MPI_INFO_NULL);
  MPI_File_read_all(Myfile,T.pointer(),1,memtype,&status);
  MPI_File_close(&Myfile); 
  T.Update_Ghosts();
}

void communicator::read(tensor1 &T,char *adrs)
{
  MPI_File_open(MPI_COMM_WORLD,adrs,MPI_MODE_RDONLY,MPI_INFO_NULL,&Myfile);
  MPI_File_set_view(Myfile,0,MPI_DOUBLE,filetype,"native",MPI_INFO_NULL);
  MPI_File_read_all(Myfile,T.x.pointer(),1,memtype,&status);
  MPI_File_set_view(Myfile,TNOD,MPI_DOUBLE,filetype,"native",MPI_INFO_NULL);
  MPI_File_read_all(Myfile,T.y.pointer(),1,memtype,&status);
  MPI_File_set_view(Myfile,2*TNOD,MPI_DOUBLE,filetype,"native",MPI_INFO_NULL);
  MPI_File_read_all(Myfile,T.z.pointer(),1,memtype,&status);
  MPI_File_close(&Myfile); 
  T.Update_Ghosts();
}

void communicator::read_S(tensor1 &T,char *adrs)
{
  MPI_File_open(pc_->XMates(),adrs,MPI_MODE_RDONLY,MPI_INFO_NULL,&Myfile);
  MPI_File_set_view(Myfile,0,MPI_DOUBLE,filetype_S,"native",MPI_INFO_NULL);
  MPI_File_read_all(Myfile,T.x.pointer(),1,memtype_S,&status);
  MPI_File_set_view(Myfile,TNOD_S,MPI_DOUBLE,filetype_S,"native",MPI_INFO_NULL);
  MPI_File_read_all(Myfile,T.y.pointer(),1,memtype_S,&status);
  MPI_File_set_view(Myfile,2*TNOD_S,MPI_DOUBLE,filetype_S,"native",MPI_INFO_NULL);
  MPI_File_read_all(Myfile,T.z.pointer(),1,memtype_S,&status);
  MPI_File_close(&Myfile); 
}

void communicator::write_S(tensor0 &T,char *adrs)
{
  MPI_File_open(pc_->XMates(),adrs,MPI_MODE_CREATE|MPI_MODE_WRONLY,MPI_INFO_NULL,&Myfile);
  MPI_File_set_view(Myfile,0,MPI_DOUBLE,filetype_S,"native",MPI_INFO_NULL);
  MPI_File_write_all(Myfile,T.pointer(),1,memtype_S,&status);
  MPI_File_close(&Myfile);
}


void communicator::write(tensor0 &T,char *adrs)
{
  MPI_File_open(MPI_COMM_WORLD,adrs,MPI_MODE_CREATE|MPI_MODE_WRONLY,MPI_INFO_NULL,&Myfile);
  MPI_File_set_view(Myfile,0,MPI_DOUBLE,filetype,"native",MPI_INFO_NULL);
  MPI_File_write_all(Myfile,T.pointer(),1,memtype,&status);
  MPI_File_close(&Myfile);
}

void communicator::write(tensor1 &T,char *adrs)
{
  MPI_File_open(MPI_COMM_WORLD,adrs,MPI_MODE_CREATE|MPI_MODE_WRONLY,MPI_INFO_NULL,&Myfile);
  MPI_File_set_view(Myfile,0,MPI_DOUBLE,filetype,"native",MPI_INFO_NULL);
  MPI_File_write_all(Myfile,T.x.pointer(),1,memtype,&status);
  MPI_File_set_view(Myfile,TNOD,MPI_DOUBLE,filetype,"native",MPI_INFO_NULL);
  MPI_File_write_all(Myfile,T.y.pointer(),1,memtype,&status);
  MPI_File_set_view(Myfile,2*TNOD,MPI_DOUBLE,filetype,"native",MPI_INFO_NULL);
  MPI_File_write_all(Myfile,T.z.pointer(),1,memtype,&status);
  MPI_File_close(&Myfile);
}

void communicator::write_S(tensor1 &T,char *adrs)
{
  MPI_File_open(pc_->XMates(),adrs,MPI_MODE_CREATE|MPI_MODE_WRONLY,MPI_INFO_NULL,&Myfile);
  MPI_File_set_view(Myfile,0,MPI_DOUBLE,filetype_S,"native",MPI_INFO_NULL);
  MPI_File_write_all(Myfile,T.x.pointer(),1,memtype_S,&status);
  MPI_File_set_view(Myfile,TNOD_S,MPI_DOUBLE,filetype_S,"native",MPI_INFO_NULL);
  MPI_File_write_all(Myfile,T.y.pointer(),1,memtype_S,&status);
  MPI_File_set_view(Myfile,2*TNOD_S,MPI_DOUBLE,filetype_S,"native",MPI_INFO_NULL);
  MPI_File_write_all(Myfile,T.z.pointer(),1,memtype_S,&status);
  MPI_File_close(&Myfile);
}

void communicator::write_reduced(tensor0 &T,char *adrs)
{
  if (pc_->IsRootX())
    {
      MPI_File_open(pc_->YZMates(),adrs,MPI_MODE_CREATE|MPI_MODE_APPEND|MPI_MODE_WRONLY,MPI_INFO_NULL,&Myfile);
      MPI_Offset offset;
      MPI_File_get_size(Myfile,&offset);
      offset+=size_->il()*sizeof(double);
      MPI_File_write_at_all(Myfile,offset,&T.ptr_reduce_[bs],Nx-2*bs,MPI_DOUBLE,MPI_STATUS_IGNORE);
      MPI_File_close(&Myfile);
    }
}




void communicator::SEND_RECV_BLOCKING(tensor0 &T)
{
  MPI_Request request_send[2];
  MPI_Request request_recv[2];
  MPI_Status stat[2];
  //////////////FIRST LEFT AND RIGHT/////////////////
  MPI_Irecv(Rbuf[0],bs*(Ny-2*bs)*(Nz-2*bs),MPI_DOUBLE,procRIGHT,1,MPI_COMM_WORLD,&request_recv[0]);
  MPI_Irecv(Rbuf[1],bs*(Ny-2*bs)*(Nz-2*bs),MPI_DOUBLE,procLEFT,0,MPI_COMM_WORLD,&request_recv[1]);
  int c=0;
  for (int i(Nx-2*bs);i<Nx-bs;i++)
    for (int j(bs);j<Ny-bs;j++)
      for (int k(bs);k<Nz-bs;k++)
	Sbuf[0][c++]=T(i,j,k);
  MPI_Isend(Sbuf[0],bs*(Ny-2*bs)*(Nz-2*bs),MPI_DOUBLE,procRIGHT,0,MPI_COMM_WORLD,&request_send[0]);

  c=0;
  for (int i(bs);i<2*bs;i++)
    for (int j(bs);j<Ny-bs;j++)
      for (int k(bs);k<Nz-bs;k++)
	Sbuf[1][c++]=T(i,j,k);
  MPI_Isend(Sbuf[1],bs*(Ny-2*bs)*(Nz-2*bs),MPI_DOUBLE,procLEFT,1,MPI_COMM_WORLD,&request_send[1]);

  MPI_Waitall(2,&request_recv[0],&stat[0]);

  c=0;
  for (int i(Nx-bs);i<Nx;i++)
    for (int j(bs);j<Ny-bs;j++)
      for (int k(bs);k<Nz-bs;k++)
	T(i,j,k)=Rbuf[0][c++];
  c=0;
  for (int i(0);i<bs;i++)
    for (int j(bs);j<Ny-bs;j++)
      for (int k(bs);k<Nz-bs;k++)
	T(i,j,k)= Rbuf[1][c++];
  MPI_Waitall(2,&request_send[0],&stat[0]);


  /////////////////Second for TOP AND BOTTOM////////
  MPI_Irecv(Rbuf[0],(Nx)*bs*(Nz-2*bs),MPI_DOUBLE,procTOP,1,MPI_COMM_WORLD,&request_recv[0]);
  MPI_Irecv(Rbuf[1],(Nx)*bs*(Nz-2*bs),MPI_DOUBLE,procBOT,0,MPI_COMM_WORLD,&request_recv[1]);
  c=0;
  for (int i(0);i<Nx;i++)
    for (int j(Ny-2*bs);j<Ny-bs;j++)
      for (int k(bs);k<Nz-bs;k++)
	Sbuf[0][c++]=T(i,j,k);
  MPI_Isend(Sbuf[0],(Nx)*bs*(Nz-2*bs),MPI_DOUBLE,procTOP,0,MPI_COMM_WORLD,&request_send[0]);
  
  c=0;
  for (int i(0);i<Nx;i++)
    for (int j(bs);j<2*bs;j++)
      for (int k(bs);k<Nz-bs;k++)
	Sbuf[1][c++]=T(i,j,k);
  MPI_Isend(Sbuf[1],(Nx)*bs*(Nz-2*bs),MPI_DOUBLE,procBOT,1,MPI_COMM_WORLD,&request_send[1]);
  
  MPI_Waitall(2,&request_recv[0],&stat[0]);
  c=0;
  for (int i(0);i<Nx;i++)
    for (int j(Ny-bs);j<Ny;j++)
      for (int k(bs);k<Nz-bs;k++)
	T(i,j,k)= Rbuf[0][c++];

  c=0;
  for (int i(0);i<Nx;i++)
    for (int j(0);j<bs;j++)
      for (int k(bs);k<Nz-bs;k++)
	T(i,j,k)= Rbuf[1][c++];
  
  MPI_Waitall(2,&request_send[0],&stat[0]);
  
  
  ///////////////////////////THIRD FOR FRONT AND REAR/////////////////////
   MPI_Irecv(Rbuf[0],(Nx)*(Ny)*bs,MPI_DOUBLE,procFRONT,1,MPI_COMM_WORLD,&request_recv[0]);
   MPI_Irecv(Rbuf[1],(Nx)*(Ny)*bs,MPI_DOUBLE,procREAR,0,MPI_COMM_WORLD,&request_recv[1]);
   c=0;;
   for (int i(0);i<Nx;i++)
     for (int j(0);j<Ny;j++)
       for (int k(Nz-2*bs);k<Nz-bs;k++)
	 Sbuf[0][c++]=T(i,j,k);

   MPI_Isend(Sbuf[0],(Nx)*(Ny)*bs,MPI_DOUBLE,procFRONT,0,MPI_COMM_WORLD,&request_send[0]);
   
   c=0;
   for (int i(0);i<Nx;i++)
     for (int j(0);j<Ny;j++)
       for (int k(bs);k<2*bs;k++)
	 Sbuf[1][c++]=T(i,j,k);

   MPI_Isend(Sbuf[1],(Nx)*(Ny)*bs,MPI_DOUBLE,procREAR,1,MPI_COMM_WORLD,&request_send[1]);
   MPI_Waitall(2,&request_recv[0],&stat[0]);
   c=0;;
   for (int i(0);i<Nx;i++)
     for (int j(0);j<Ny;j++)
       for (int k(Nz-bs);k<Nz;k++)
	 T(i,j,k)=Rbuf[0][c++];

   c=0;
   for (int i(0);i<Nx;i++)
     for (int j(0);j<Ny;j++)
       for (int k(0);k<bs;k++)
	 T(i,j,k)= Rbuf[1][c++];
   MPI_Waitall(2,&request_send[0],&stat[0]);
}

void communicator::SEND_RECV_BLOCKING_INOUT(tensor0 &T)
{
  MPI_Request request_send[2];
  MPI_Request request_recv[2];
  MPI_Status stat[2];
  //////////////FIRST LEFT AND RIGHT/////////////////
  if (!pc_->IsOut()) MPI_Irecv(Rbuf[0],bs*(Ny-2*bs)*(Nz-2*bs),MPI_DOUBLE,procRIGHT,1,MPI_COMM_WORLD,&request_recv[0]);
  if (!pc_->IsIn()) MPI_Irecv(Rbuf[1],bs*(Ny-2*bs)*(Nz-2*bs),MPI_DOUBLE,procLEFT,0,MPI_COMM_WORLD,&request_recv[1]);
  int c;
  if (!pc_->IsOut())
    {
      c=0;
      for (int i(Nx-2*bs);i<Nx-bs;i++)
	for (int j(bs);j<Ny-bs;j++)
	  for (int k(bs);k<Nz-bs;k++)
	    Sbuf[0][c++]=T(i,j,k);
      MPI_Isend(Sbuf[0],bs*(Ny-2*bs)*(Nz-2*bs),MPI_DOUBLE,procRIGHT,0,MPI_COMM_WORLD,&request_send[0]);
    }

  if (!pc_->IsIn())
    {
      c=0;
      for (int i(bs);i<2*bs;i++)
	for (int j(bs);j<Ny-bs;j++)
	  for (int k(bs);k<Nz-bs;k++)
	    Sbuf[1][c++]=T(i,j,k);
      MPI_Isend(Sbuf[1],bs*(Ny-2*bs)*(Nz-2*bs),MPI_DOUBLE,procLEFT,1,MPI_COMM_WORLD,&request_send[1]);
    }
  if ((!pc_->IsIn())&&(!pc_->IsOut()))  MPI_Waitall(2,&request_recv[0],&stat[0]);
  if ((pc_->IsIn())&&(!pc_->IsOut())) MPI_Waitall(1,&request_recv[0],&stat[0]);
  if ((!pc_->IsIn())&&(pc_->IsOut())) MPI_Waitall(1,&request_recv[1],&stat[1]);

  if (!pc_->IsOut())
    {
      c=0;
      for (int i(Nx-bs);i<Nx;i++)
	for (int j(bs);j<Ny-bs;j++)
	  for (int k(bs);k<Nz-bs;k++)
	    T(i,j,k)=Rbuf[0][c++];
    }
  if (!pc_->IsIn())
    {
      c=0;
      for (int i(0);i<bs;i++)
	for (int j(bs);j<Ny-bs;j++)
	  for (int k(bs);k<Nz-bs;k++)
	    T(i,j,k)= Rbuf[1][c++];
    }

  if ((!pc_->IsIn())&&(!pc_->IsOut()))  MPI_Waitall(2,&request_send[0],&stat[0]);
  if ((pc_->IsIn())&&(!pc_->IsOut())) MPI_Waitall(1,&request_send[0],&stat[0]);
  if ((!pc_->IsIn())&&(pc_->IsOut())) MPI_Waitall(1,&request_send[1],&stat[1]);
  /////////////////Second for TOP AND BOTTOM////////
  MPI_Irecv(Rbuf[0],(Nx)*bs*(Nz-2*bs),MPI_DOUBLE,procTOP,1,MPI_COMM_WORLD,&request_recv[0]);
  MPI_Irecv(Rbuf[1],(Nx)*bs*(Nz-2*bs),MPI_DOUBLE,procBOT,0,MPI_COMM_WORLD,&request_recv[1]);
  c=0;
  for (int i(0);i<Nx;i++)
    for (int j(Ny-2*bs);j<Ny-bs;j++)
      for (int k(bs);k<Nz-bs;k++)
	Sbuf[0][c++]=T(i,j,k);
  MPI_Isend(Sbuf[0],(Nx)*bs*(Nz-2*bs),MPI_DOUBLE,procTOP,0,MPI_COMM_WORLD,&request_send[0]);
  
  c=0;
  for (int i(0);i<Nx;i++)
    for (int j(bs);j<2*bs;j++)
      for (int k(bs);k<Nz-bs;k++)
	Sbuf[1][c++]=T(i,j,k);
  MPI_Isend(Sbuf[1],(Nx)*bs*(Nz-2*bs),MPI_DOUBLE,procBOT,1,MPI_COMM_WORLD,&request_send[1]);
  
  MPI_Waitall(2,&request_recv[0],&stat[0]);
  c=0;
  for (int i(0);i<Nx;i++)
    for (int j(Ny-bs);j<Ny;j++)
      for (int k(bs);k<Nz-bs;k++)
	T(i,j,k)= Rbuf[0][c++];

  c=0;
  for (int i(0);i<Nx;i++)
    for (int j(0);j<bs;j++)
      for (int k(bs);k<Nz-bs;k++)
	T(i,j,k)= Rbuf[1][c++];
  
  MPI_Waitall(2,&request_send[0],&stat[0]);
  
  
  ///////////////////////////THIRD FOR FRONT AND REAR/////////////////////
   MPI_Irecv(Rbuf[0],(Nx)*(Ny)*bs,MPI_DOUBLE,procFRONT,1,MPI_COMM_WORLD,&request_recv[0]);
   MPI_Irecv(Rbuf[1],(Nx)*(Ny)*bs,MPI_DOUBLE,procREAR,0,MPI_COMM_WORLD,&request_recv[1]);
   c=0;;
   for (int i(0);i<Nx;i++)
     for (int j(0);j<Ny;j++)
       for (int k(Nz-2*bs);k<Nz-bs;k++)
	 Sbuf[0][c++]=T(i,j,k);

   MPI_Isend(Sbuf[0],(Nx)*(Ny)*bs,MPI_DOUBLE,procFRONT,0,MPI_COMM_WORLD,&request_send[0]);
   
   c=0;
   for (int i(0);i<Nx;i++)
     for (int j(0);j<Ny;j++)
       for (int k(bs);k<2*bs;k++)
	 Sbuf[1][c++]=T(i,j,k);

   MPI_Isend(Sbuf[1],(Nx)*(Ny)*bs,MPI_DOUBLE,procREAR,1,MPI_COMM_WORLD,&request_send[1]);
   MPI_Waitall(2,&request_recv[0],&stat[0]);
   c=0;;
   for (int i(0);i<Nx;i++)
     for (int j(0);j<Ny;j++)
       for (int k(Nz-bs);k<Nz;k++)
	 T(i,j,k)=Rbuf[0][c++];

   c=0;
   for (int i(0);i<Nx;i++)
     for (int j(0);j<Ny;j++)
       for (int k(0);k<bs;k++)
	 T(i,j,k)= Rbuf[1][c++];
   MPI_Waitall(2,&request_send[0],&stat[0]);
}

void communicator::SEND_RECV_BLOCKING(tensor1 &T)
{
  SEND_RECV_BLOCKING(T.x);
  SEND_RECV_BLOCKING(T.y);
  SEND_RECV_BLOCKING(T.z);
}

void communicator::SEND_RECV_BLOCKING_INOUT(tensor1 &T)
{
  SEND_RECV_BLOCKING_INOUT(T.x);
  SEND_RECV_BLOCKING_INOUT(T.y);
  SEND_RECV_BLOCKING_INOUT(T.z);
}

void communicator::Sequential_Begin(int &x)
{
  if (pc_->RANK()>0) MPI_Recv(&x,1,MPI_INT,pc_->RANK()-1,pc_->RANK()-1,MPI_COMM_WORLD,&status);
}

void communicator::Sequential_End(int &x)
{
  if (pc_->RANK()<pc_->TOT()-1) MPI_Send(&x,1,MPI_INT,pc_->RANK()+1,pc_->RANK(),MPI_COMM_WORLD);
}

void communicator::SEND_RECV_BLOCKING_CUM(tensor0 &T)
{
  MPI_Request request_send[2];
  MPI_Request request_recv[2];
  MPI_Status stat[2];
  ///////////////////////////FIRST FOR FRONT AND REAR/////////////////////
   MPI_Irecv(Rbuf[0],(Nx)*(Ny)*bs,MPI_DOUBLE,procFRONT,1,MPI_COMM_WORLD,&request_recv[0]);
   MPI_Irecv(Rbuf[1],(Nx)*(Ny)*bs,MPI_DOUBLE,procREAR,0,MPI_COMM_WORLD,&request_recv[1]);
   int c=0;
   for (int i(0);i<Nx;i++)
     for (int j(0);j<Ny;j++)
       for (int k(Nz-bs);k<Nz;k++)
	 Sbuf[0][c++]=T(i,j,k);

   MPI_Isend(Sbuf[0],(Nx)*(Ny)*bs,MPI_DOUBLE,procFRONT,0,MPI_COMM_WORLD,&request_send[0]);

   c=0;
   for (int i(0);i<Nx;i++)
     for (int j(0);j<Ny;j++)
       for (int k(0);k<bs;k++)
	 Sbuf[1][c++]=T(i,j,k);

   MPI_Isend(Sbuf[1],(Nx)*(Ny)*bs,MPI_DOUBLE,procREAR,1,MPI_COMM_WORLD,&request_send[1]);
   MPI_Waitall(2,&request_recv[0],&stat[0]);
   c=0;
   for (int i(0);i<Nx;i++)
     for (int j(0);j<Ny;j++)
       for (int k(Nz-2*bs);k<Nz-bs;k++)
	 T(i,j,k)+=Rbuf[0][c++];
   c=0;
   for (int i(0);i<Nx;i++)
     for (int j(0);j<Ny;j++)
       for (int k(bs);k<2*bs;k++)
	 T(i,j,k)+=Rbuf[1][c++];
   MPI_Waitall(2,&request_send[0],&stat[0]);
  /////////////////Second for TOP AND BOTTOM////////
  MPI_Irecv(Rbuf[0],(Nx)*bs*(Nz-2*bs),MPI_DOUBLE,procTOP,1,MPI_COMM_WORLD,&request_recv[0]);
  MPI_Irecv(Rbuf[1],(Nx)*bs*(Nz-2*bs),MPI_DOUBLE,procBOT,0,MPI_COMM_WORLD,&request_recv[1]);
  c=0;
  for (int i(0);i<Nx;i++)
    for (int j(Ny-bs);j<Ny;j++)
      for (int k(bs);k<Nz-bs;k++)
	Sbuf[0][c++]=T(i,j,k);
  MPI_Isend(Sbuf[0],(Nx)*bs*(Nz-2*bs),MPI_DOUBLE,procTOP,0,MPI_COMM_WORLD,&request_send[0]);
  c=0;
  for (int i(0);i<Nx;i++)
    for (int j(0);j<bs;j++)
      for (int k(bs);k<Nz-bs;k++)
	Sbuf[1][c++]=T(i,j,k);
  MPI_Isend(Sbuf[1],(Nx)*bs*(Nz-2*bs),MPI_DOUBLE,procBOT,1,MPI_COMM_WORLD,&request_send[1]);
  MPI_Waitall(2,&request_recv[0],&stat[0]);
  c=0;
  for (int i(0);i<Nx;i++)
    for (int j(Ny-2*bs);j<Ny-bs;j++)
      for (int k(bs);k<Nz-bs;k++)
	T(i,j,k)+=Rbuf[0][c++];
  c=0;
  for (int i(0);i<Nx;i++)
    for (int j(bs);j<2*bs;j++)
      for (int k(bs);k<Nz-bs;k++)
	T(i,j,k)+=Rbuf[1][c++];
  MPI_Waitall(2,&request_send[0],&stat[0]);
  //////////////FIRST LEFT AND RIGHT/////////////////
  MPI_Irecv(Rbuf[0],bs*(Ny-2*bs)*(Nz-2*bs),MPI_DOUBLE,procRIGHT,1,MPI_COMM_WORLD,&request_recv[0]);
  MPI_Irecv(Rbuf[1],bs*(Ny-2*bs)*(Nz-2*bs),MPI_DOUBLE,procLEFT,0,MPI_COMM_WORLD,&request_recv[1]);
  c=0;
  for (int i(Nx-bs);i<Nx;i++)
    for (int j(bs);j<Ny-bs;j++)
      for (int k(bs);k<Nz-bs;k++)
	Sbuf[0][c++]=T(i,j,k);
  MPI_Isend(Sbuf[0],bs*(Ny-2*bs)*(Nz-2*bs),MPI_DOUBLE,procRIGHT,0,MPI_COMM_WORLD,&request_send[0]);
  c=0;
  for (int i(0);i<bs;i++)
    for (int j(bs);j<Ny-bs;j++)
      for (int k(bs);k<Nz-bs;k++)
	Sbuf[1][c++]=T(i,j,k);
  MPI_Isend(Sbuf[1],bs*(Ny-2*bs)*(Nz-2*bs),MPI_DOUBLE,procLEFT,1,MPI_COMM_WORLD,&request_send[1]);
  MPI_Waitall(2,&request_recv[0],&stat[0]);
  c=0;
  for (int i(Nx-2*bs);i<Nx-bs;i++)
    for (int j(bs);j<Ny-bs;j++)
      for (int k(bs);k<Nz-bs;k++)
	T(i,j,k)+=Rbuf[0][c++];
  c=0;
  for (int i(bs);i<2*bs;i++)
    for (int j(bs);j<Ny-bs;j++)
      for (int k(bs);k<Nz-bs;k++)
	T(i,j,k)+=Rbuf[1][c++];
  MPI_Waitall(2,&request_send[0],&stat[0]);
}

void communicator::SEND_RECV_BLOCKING_CUM(tensor1 &T)
{
  SEND_RECV_BLOCKING_CUM(T.x);
  SEND_RECV_BLOCKING_CUM(T.y);
  SEND_RECV_BLOCKING_CUM(T.z);
}

void communicator::SEND_RECV_BLOCKING_CUM_INOUT(tensor0 &T)
{
  MPI_Request request_send[2];
  MPI_Request request_recv[2];
  MPI_Status stat[2];
  ///////////////////////////FIRST FOR FRONT AND REAR/////////////////////
   MPI_Irecv(Rbuf[0],(Nx)*(Ny)*bs,MPI_DOUBLE,procFRONT,1,MPI_COMM_WORLD,&request_recv[0]);
   MPI_Irecv(Rbuf[1],(Nx)*(Ny)*bs,MPI_DOUBLE,procREAR,0,MPI_COMM_WORLD,&request_recv[1]);
   int c=0;
   for (int i(0);i<Nx;i++)
     for (int j(0);j<Ny;j++)
       for (int k(Nz-bs);k<Nz;k++)
	 Sbuf[0][c++]=T(i,j,k);

   MPI_Isend(Sbuf[0],(Nx)*(Ny)*bs,MPI_DOUBLE,procFRONT,0,MPI_COMM_WORLD,&request_send[0]);

   c=0;
   for (int i(0);i<Nx;i++)
     for (int j(0);j<Ny;j++)
       for (int k(0);k<bs;k++)
	 Sbuf[1][c++]=T(i,j,k);

   MPI_Isend(Sbuf[1],(Nx)*(Ny)*bs,MPI_DOUBLE,procREAR,1,MPI_COMM_WORLD,&request_send[1]);
   MPI_Waitall(2,&request_recv[0],&stat[0]);
   c=0;
   for (int i(0);i<Nx;i++)
     for (int j(0);j<Ny;j++)
       for (int k(Nz-2*bs);k<Nz-bs;k++)
	 T(i,j,k)+=Rbuf[0][c++];
   c=0;
   for (int i(0);i<Nx;i++)
     for (int j(0);j<Ny;j++)
       for (int k(bs);k<2*bs;k++)
	 T(i,j,k)+=Rbuf[1][c++];
   MPI_Waitall(2,&request_send[0],&stat[0]);
  /////////////////Second for TOP AND BOTTOM////////
  MPI_Irecv(Rbuf[0],(Nx)*bs*(Nz-2*bs),MPI_DOUBLE,procTOP,1,MPI_COMM_WORLD,&request_recv[0]);
  MPI_Irecv(Rbuf[1],(Nx)*bs*(Nz-2*bs),MPI_DOUBLE,procBOT,0,MPI_COMM_WORLD,&request_recv[1]);
  c=0;
  for (int i(0);i<Nx;i++)
    for (int j(Ny-bs);j<Ny;j++)
      for (int k(bs);k<Nz-bs;k++)
	Sbuf[0][c++]=T(i,j,k);
  MPI_Isend(Sbuf[0],(Nx)*bs*(Nz-2*bs),MPI_DOUBLE,procTOP,0,MPI_COMM_WORLD,&request_send[0]);
  c=0;
  for (int i(0);i<Nx;i++)
    for (int j(0);j<bs;j++)
      for (int k(bs);k<Nz-bs;k++)
	Sbuf[1][c++]=T(i,j,k);
  MPI_Isend(Sbuf[1],(Nx)*bs*(Nz-2*bs),MPI_DOUBLE,procBOT,1,MPI_COMM_WORLD,&request_send[1]);
  MPI_Waitall(2,&request_recv[0],&stat[0]);
  c=0;
  for (int i(0);i<Nx;i++)
    for (int j(Ny-2*bs);j<Ny-bs;j++)
      for (int k(bs);k<Nz-bs;k++)
	T(i,j,k)+=Rbuf[0][c++];
  c=0;
  for (int i(0);i<Nx;i++)
    for (int j(bs);j<2*bs;j++)
      for (int k(bs);k<Nz-bs;k++)
	T(i,j,k)+=Rbuf[1][c++];
  MPI_Waitall(2,&request_send[0],&stat[0]);
  //////////////FIRST LEFT AND RIGHT/////////////////
  if (!pc_->IsOut()) MPI_Irecv(Rbuf[0],bs*(Ny-2*bs)*(Nz-2*bs),MPI_DOUBLE,procRIGHT,1,MPI_COMM_WORLD,&request_recv[0]);
  if (!pc_->IsIn()) MPI_Irecv(Rbuf[1],bs*(Ny-2*bs)*(Nz-2*bs),MPI_DOUBLE,procLEFT,0,MPI_COMM_WORLD,&request_recv[1]);
  if (!pc_->IsOut())
    {
      c=0;
      for (int i(Nx-bs);i<Nx;i++)
	for (int j(bs);j<Ny-bs;j++)
	  for (int k(bs);k<Nz-bs;k++)
	    Sbuf[0][c++]=T(i,j,k);
      MPI_Isend(Sbuf[0],bs*(Ny-2*bs)*(Nz-2*bs),MPI_DOUBLE,procRIGHT,0,MPI_COMM_WORLD,&request_send[0]);
    }
  if (!pc_->IsIn())
    {
      c=0;
      for (int i(0);i<bs;i++)
	for (int j(bs);j<Ny-bs;j++)
	  for (int k(bs);k<Nz-bs;k++)
	    Sbuf[1][c++]=T(i,j,k);
      MPI_Isend(Sbuf[1],bs*(Ny-2*bs)*(Nz-2*bs),MPI_DOUBLE,procLEFT,1,MPI_COMM_WORLD,&request_send[1]);
    }
  
  if ((!pc_->IsIn())&&(!pc_->IsOut()))  MPI_Waitall(2,&request_recv[0],&stat[0]);
  if ((pc_->IsIn())&&(!pc_->IsOut())) MPI_Waitall(1,&request_recv[0],&stat[0]);
  if ((!pc_->IsIn())&&(pc_->IsOut())) MPI_Waitall(1,&request_recv[1],&stat[1]);
  if (!pc_->IsOut())
    {
      c=0;
      for (int i(Nx-2*bs);i<Nx-bs;i++)
	for (int j(bs);j<Ny-bs;j++)
	  for (int k(bs);k<Nz-bs;k++)
	    T(i,j,k)+=Rbuf[0][c++];
    }
  if (!pc_->IsIn())
    {
      c=0;
      for (int i(bs);i<2*bs;i++)
	for (int j(bs);j<Ny-bs;j++)
	  for (int k(bs);k<Nz-bs;k++)
	    T(i,j,k)+=Rbuf[1][c++];
    }
  if ((!pc_->IsIn())&&(!pc_->IsOut()))  MPI_Waitall(2,&request_send[0],&stat[0]);
  if ((pc_->IsIn())&&(!pc_->IsOut())) MPI_Waitall(1,&request_send[0],&stat[0]);
  if ((!pc_->IsIn())&&(pc_->IsOut())) MPI_Waitall(1,&request_send[1],&stat[1]);
}

void communicator::SEND_RECV_BLOCKING_CUM_INOUT(tensor1 &T)
{
  SEND_RECV_BLOCKING_CUM_INOUT(T.x);
  SEND_RECV_BLOCKING_CUM_INOUT(T.y);
  SEND_RECV_BLOCKING_CUM_INOUT(T.z);
}

void communicator::REDUCEX(tensor0 &T)
{
  for (int i(0);i<Nx;i++)
    {
      Sbuf[0][i]=0;
      for (int j(bs);j<Ny-bs;j++)
        for (int k(bs);k<Nz-bs;k++)
          Sbuf[0][i]+=T(i,j,k);
      Sbuf[0][i]/=(size_->Ny_tot()*size_->Nz_tot());
    }
  MPI_Allreduce(Sbuf[0],T.ptr_reduce_,Nx,MPI_DOUBLE,MPI_SUM,pc_->XMates());
}

