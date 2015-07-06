close all;
clear all;
clc;


Lx=0.16; Ly=0.04;Lz=0.04;
nx=512;ny=128;nz=128;
dx=Lx/nx; dy=Ly/ny;dz=Lz/nz;

x=dx/2:Lx/nx:Lx-dx/2;
y=dy/2:Ly/ny:Ly-dy/2;
z=dz/2:Lz/nz:Lz-dz/2;

%path=['../InOut_PSAAP_Nominal_Bidisperse_sim4/data_IO/'];
path=['./data_IO/'];
%path=['../InOut_PSAAP_Nominal_U0_2_St_0.4_sim3/data_IO/'];
%path=['../InOut_PSAAP_Nominal_Polydisperse_sim13/data_IO/'];

name_r=dir([path 'RU*']);
name_rho=dir([path 'Rho*']);

u_g_inlet=[];u_g_outlet=[];

for i=1:41%length(name_r)
    US(i)=str2num(name_r(i).name(4:end-4));
end

[N,M]=sort(US);
%name_T.name=name_T(ind).name


for j=11:41%length(name_T)
     j 

filename=[path name_r(M(j)).name]
fileid= fopen(filename);
A=fread(fileid,'double');
RU=reshape(A,nx,ny,nz,3);


filename=[path name_rho(M(j)).name]
fileid= fopen(filename);
A=fread(fileid,'double');
Rho=reshape(A,nx,ny,nz);

u=squeeze(RU(:,:,:,1))./Rho;

u_g(j,:)=mean(squeeze(mean(u,3)),2);


u_g_inlet=[u_g_inlet; reshape(squeeze(u(end,:,:)),ny*nz,1)];

u_g_outlet=[u_g_outlet; reshape(squeeze(u(end,:,:)),ny*nz,1)];

%figure
%plot(xc,temp,'LineWidth',2.5); hold on;

end   
figure;
plot(u_g(11:end,:)');

u_g_mean=mean(u_g(11:end,:),1);
hold on;plot(u_g_mean,'.');
save './gas_vel_mono_sim5.mat' u_g_mean u_g_inlet u_g_outlet;

