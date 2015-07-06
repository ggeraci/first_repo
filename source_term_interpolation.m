close all;
clear all;
clc;

u0=2.065; % m/s
T0=300;
rho_p=8900;
rho_f=1.2;
dp=1.02e-5;
mu=1.95e-5;
mp=rho_p*pi/6*dp^3;
Ap=pi/4*dp^2;
n0=3.297e10;
g=0;
Cvp=450;
Nu=2;K=0.0275;
h=Nu*K/dp;

Cv=723;
Cp=1012;
Pr=0.7;
eps=0.4;
I0=6.86e6;

Lx=0.16; Ly=0.04;Lz=0.04;
nx=512;ny=128;nz=128;
dx=Lx/nx; dy=Ly/ny;dz=Lz/nz;

x=dx/2:Lx/nx:Lx-dx/2;
y=dy/2:Ly/ny:Ly-dy/2;
z=dz/2:Lz/nz:Lz-dz/2;

path=['./data_IO/'];
Lx=0.16; Ly=0.04; Lz=0.04;
nf_f=120000;
Delta_t=0.00005;
nf=[2000:2000:nf_f];
t=nf*Delta_t;

for i=2:length(nf);
    i
    
    filename=[path 'Part_T_' num2str(nf(i)) '.bin'];
    fileid= fopen(filename);
    Tp=fread(fileid,'double');

    filename=[path 'Part_x_' num2str(nf(i)) '.bin'];
    fileid= fopen(filename);
    xp=fread(fileid,'double');
    
    filename=[path 'Part_y_' num2str(nf(i)) '.bin'];
    fileid= fopen(filename);
    yp=fread(fileid,'double');
    
    
    filename=[path 'Part_z_' num2str(nf(i)) '.bin'];
    fileid= fopen(filename);
    zp=fread(fileid,'double');
    
    P0=rho_f*287.058*T0;     
    filename=[path 'Rho_' num2str(nf(i)) '.bin'];
    fileid= fopen(filename);
    A=fread(fileid,'double');
    rho=reshape(A,nx,ny,nz);

    Tg=P0./rho/287.058;
    
    % building the new mesh with walls:
    x_n=[0 x Lx];
    y_n=[0 y Ly];
    z_n=[0 z Lz];
    
    T_nx0=(3*Tg(1,:,:)-Tg(2,:,:))/2;
    T_nx1=(3*Tg(end,:,:)-Tg(end-1,:,:))/2;
    T_n=cat(1,T_nx0,Tg,T_nx1);
    
    T_ny0=(3*T_n(:,1,:)-T_n(:,2,:))/2;
    T_ny1=(3*T_n(:,end,:)-T_n(:,end-1,:))/2;
    
    T_n=cat(2,T_ny0,T_n,T_ny1);
   
    T_nz0=(3*T_n(:,:,1)-T_n(:,:,2))/2;
    T_nz1=(3*T_n(:,:,end)-T_n(:,:,end-1))/2;
    
    T_n=cat(3,T_nz0,T_n,T_nz1);
   

   
    for j=1:length(xp)
        j; 
        ind=find(xp(j)<x_n);
        indx=ind(1);
        x_l=(xp(j)-x_n(indx-1))/(x_n(indx)-x_n(indx-1));
        ii=indx;
        
        ind=find(yp(j)<y_n);
        indy=ind(1);
        y_l=(yp(j)-y_n(indy-1))/(y_n(indy)-y_n(indy-1));
        jj=indy;
        
        ind=find(zp(j)<z_n);
        indz=ind(1);
        z_l=(zp(j)-z_n(indz-1))/(z_n(indz)-z_n(indz-1));
        kk=indz;
        
        Tg_int(j)=x_l*y_l*z_l*T_n(ii,jj,kk)+ ...
            (1-x_l)*y_l*z_l*T_n(ii-1,jj,kk)+x_l*(1-y_l)*z_l*T_n(ii,jj-1,kk)+x_l*y_l*(1-z_l)*T_n(ii,jj,kk-1)+...
            (1-x_l)*(1-y_l)*z_l*T_n(ii-1,jj-1,kk)+(1-x_l)*y_l*(1-z_l)*T_n(ii-1,jj,kk-1)+x_l*(y_l-1)*(z_l-1)*T_n(ii,jj-1,kk-1)+...
            (1-x_l)*(1-y_l)*(1-z_l)*T_n(ii-1,jj-1,kk-1);
    

    end
    
    sum_DT(i)=sum(Tg_int-Tp');
    clear Tg_int;
end

figure;
plot(sum_DT,'o');
save './sum_DT_st_1_new.mat' sum_DT t nf;