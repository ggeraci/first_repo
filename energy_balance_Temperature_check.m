close all;
clear all;

Lx=0.16; Ly=0.04;Lz=0.04;




%__________________________
%_________________________
% Sim1:
u0=0.445; % m/s
I0=1.498e6;
nf_f=62000
Delta_t=0.0004

S_1D=1.0796e+06*Lx*Ly*Lz;
H_flux_div_1D=1.0165e+06*Lx*Ly*Lz;
E_th_p_flux_grad_1D = 6.3689e+04*Lx*Ly*Lz;
%__________________________
%_________________________
% Sim3:
u0=1.1205; % m/s
I0=3.75e6;
nf_f=66000
Delta_t=0.0002

S_1D=2.7384e6*Lx*Ly*Lz;
H_flux_div_1D=2.5591e+06*Lx*Ly*Lz;
E_th_p_flux_grad_1D = 1.8040e+05*Lx*Ly*Lz;
%__________________________

%_________________________
% Sim5:
u0=2.065; % m/s
I0=6.68e6;
np0=3.297e10;
nf_f=70000;
Delta_t=0.00005;

S_1D=4.9936e+06*Lx*Ly*Lz;
H_flux_div_1D=4.6200e+06*Lx*Ly*Lz;
E_th_p_flux_grad_1D = 3.7566e+05*Lx*Ly*Lz;
T_p_1D=720.80;
T_g_1D=600.9;


T0=300;
rho_p=8900;
rho_f=1.2;
dp=1.02e-5;
mu=1.95e-5;
nu=mu/rho_f;
mp=rho_p*pi/6*dp^3;
Ap=pi/4*dp^2;
n0=3.297e10;
g=0;
Cvp=450;
Nu=2;
Cv=723;
Cp=1012;
Pr=0.7;
eps=0.4;
kappa=Pr*nu;
K=0.0275;
k=0.0275;
h=Nu*K/dp;
Lx=0.16; Ly=0.04;Lz=0.04;
nx=512;ny=128;nz=128;
dx=Lx/nx; dy=Ly/ny;dz=Lz/nz;

x=dx/2:Lx/nx:Lx-dx/2;
y=dy/2:Ly/ny:Ly-dy/2;
z=dz/2:Lz/nz:Lz-dz/2;

%__________________________________________________________________________

clc;
path=['./data_IO/'];
%path=['/home/mrahmani/Documents/test_cases_poly_paper/test/'];
%nf_f=32000;
nf=[2000:2000:nf_f];
t=nf*Delta_t;
tr=Lx/u0;

S2=load('sum_DT_st_1_new.mat');

source=-S2.sum_DT*pi*dp^2*h;

for i=1:length(nf);
i    
filename=[path 'numbers_' num2str(nf(i)) '.dat'];
n=load(filename);
Np(i)=n(4);
np(i)=Np(i)/(Lx*Ly*Lz);
S(i)=eps*Ap*I0*np(i)*Lx*Ly*Lz;
%_________________________________________
P0=rho_f*287.058*T0;     
filename=[path 'Rho_' num2str(nf(i)) '.bin'];
fileid= fopen(filename);
A=fread(fileid,'double');
rho=reshape(A,nx,ny,nz);

Tg=P0./rho/287.058;

%_________________________________________
filename=[path 'RU_' num2str(nf(i)) '.bin'];
fileid= fopen(filename);
A=fread(fileid,'double');
RU=reshape(A,nx,ny,nz,3);


u=squeeze(RU(:,:,:,1))./rho;
v=squeeze(RU(:,:,:,2))./rho;
w=squeeze(RU(:,:,:,3))./rho;

Delta_T_g_model(i)=source(i)/(rho_f*u0*Cp*Ly*Lz);
Delta_T_g_DNS(i)=mean(mean(squeeze(Tg(end,:,:))))-mean(mean(squeeze(Tg(1,:,:))));

%H_flux_div(i)=Cp*(sum(sum(Tg(end,:,:).*rho(end,:,:).*u(end,:,:)))-sum(sum(Tg(1,:,:).*rho(1,:,:).*u(1,:,:))))*dy*dz;
%_________________________________________
filename=[path 'Part_T_' num2str(nf(i)) '.bin'];
fileid= fopen(filename);
Tp=fread(fileid,'double');

filename=[path 'Part_x_' num2str(nf(i)) '.bin'];
fileid= fopen(filename);
px=fread(fileid,'double');

filename=[path 'Part_y_' num2str(nf(i)) '.bin'];
fileid= fopen(filename);
py=fread(fileid,'double');

filename=[path 'Part_z_' num2str(nf(i)) '.bin'];
fileid= fopen(filename);
pz=fread(fileid,'double');


ind1=find(px>(Lx-dx));
ind0=find(px<(dx));

filename=[path 'Part_u_' num2str(nf(i)) '.bin'];
fileid= fopen(filename);
up=fread(fileid,'double');

filename=[path 'Part_v_' num2str(nf(i)) '.bin'];
fileid= fopen(filename);
vp=fread(fileid,'double');

filename=[path 'Part_w_' num2str(nf(i)) '.bin'];
fileid= fopen(filename);
wp=fread(fileid,'double');


Delta_T_p_model(i)=(eps*Ap*I0*Np(i)-source(i))/(mp*Cvp)/(np0*u0);
Delta_T_p_DNS(i)=mean(Tp(ind1))-mean(Tp(ind0));


% ind1=find(py>(Ly-dy));
% ind0=find(py<(dy));
% E_th_p_flux_grad_y(i)=mp*Cvp*(sum(Tp(ind1).*vp(ind1))-sum(Tp(ind0).*vp(ind0)))/(dy);
% 
% ind1=find(pz>(Lz-dz));
% ind0=find(pz<(dz));
% E_th_p_flux_grad_z(i)=mp*Cvp*(sum(Tp(ind1).*wp(ind1))-sum(Tp(ind0).*wp(ind0)))/(dz);


% T_p_mean_u(i)=sum(Tp(ind1).*up(ind1))/sum(up(ind1));
% T_p_mean(i)=mean(Tp(ind1));
end

figure;
plot(t/tr,Delta_T_p_model,'b--','LineWidth',2.5); hold on;plot(t/tr,Delta_T_p_model,'bo','MarkerSize',6,'LineWidth',2.5); hold on; plot([0 100],[T_p_1D T_p_1D],'k','LineWidth',1.5); 
plot(t/tr,Delta_T_p_DNS,'r--','LineWidth',2.5); hold on;plot(t/tr,Delta_T_p_DNS,'ro','MarkerSize',6,'LineWidth',2.5);

set(gca,'FontSize',20); xlabel('t/t_{R}','FontSize',20);ylabel('T_{p} (K)','FontSize',20);
xlim([0 55]);

figure;
plot(t/tr,Delta_T_g_model,'b--','LineWidth',2.5); hold on;plot(t/tr,Delta_T_g_model,'bo','MarkerSize',6,'LineWidth',2.5); hold on; plot([0 100],[T_g_1D T_g_1D],'k','LineWidth',1.5); 
plot(t/tr,Delta_T_g_DNS,'r--','LineWidth',2.5); hold on;plot(t/tr,Delta_T_g_DNS,'ro','MarkerSize',6,'LineWidth',2.5);

set(gca,'FontSize',20); xlabel('t/t_{R}','FontSize',20);ylabel('T_{g} (K)','FontSize',20);
xlim([0 55]);
%print -depsc 'energy_flux_balance_st_1_new.eps';

save './Delta_T_comp.mat' Delta_T_p_DNS Delta_T_g_DNS Delta_T_p_model Delta_T_g_model;
