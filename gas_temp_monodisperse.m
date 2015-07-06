close all;
clear all;
clc;


Lx=0.16; Ly=0.04;Lz=0.04;
nx=512;ny=128;nz=128;
dx=Lx/nx; dy=Ly/ny;dz=Lz/nz;

x=dx/2:Lx/nx:Lx-dx/2;
y=dy/2:Ly/ny:Ly-dy/2;
z=dz/2:Lz/nz:Lz-dz/2;

path=['./data_IO/'];
%path=['../InOut_PSAAP_Nominal_U0_2_St_0.1_sim2/data_IO/'];
%path=['../InOut_PSAAP_Nominal_U0_2_St_0.4_sim3/data_IO/'];


name_r=dir([path 'Rho*']);

T_g_inlet=[];T_g_outlet=[];

for i=1:41%length(name_r)
    TS(i)=str2num(name_r(i).name(5:end-4));
end

[N,M]=sort(TS);
%name_T.name=name_T(ind).name


for j=11:41%length(name_T)
     j 
P0=1.2*287.058*300;     
filename=[path name_r(M(j)).name]
fileid= fopen(filename);
A=fread(fileid,'double');
rho=reshape(A,nx,ny,nz);

Tg=P0./rho/287.058;

T_g(j,:)=mean(squeeze(mean(Tg,3)),2);


T_g_inlet=[T_g_inlet; reshape(squeeze(Tg(1,:,:)),ny*nz,1)];

T_g_outlet=[T_g_outlet; reshape(squeeze(Tg(end,:,:)),ny*nz,1)];

%figure
%plot(xc,temp,'LineWidth',2.5); hold on;

end   
figure;

plot(T_g(11:end,:)');

T_g_mean=mean(T_g(11:end,:),1);
hold on;plot(T_g_mean,'.');

save './gas_temp_mono_sim5.mat' T_g_mean T_g_inlet T_g_outlet;

