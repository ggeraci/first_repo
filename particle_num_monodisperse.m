close all;
clear all;
clc;


Lx=0.16; Ly=0.04;Lz=0.04;
nx=511;ny=128;nz=128;
dx=Lx/nx; dy=Ly/ny;dz=Lz/nz;

x=dx/2:Lx/nx:Lx-dx/2;
y=dy/2:Ly/ny:Ly-dy/2;
z=dz/2:Lz/nz:Lz-dz/2;

path=['./data_IO/'];
%path=['../InOut_PSAAP_Nominal_U0_2_St_0.1_sim2/data_IO/'];
%path=['../InOut_PSAAP_Nominal_U0_2_St_0.4_sim3/data_IO/'];

name_u=dir([path 'Part_u*']);
name_x=dir([path 'Part_x*']);

u_p_inlet=[];u_p_outlet=[];

for i=1:41%length(name_u)
    US(i)=str2num(name_x(i).name(8:end-4));
end

[N,M]=sort(US);
%name_T.name=name_T(ind).name


for j=10:41%length(name_u)
     j 
filename=[path name_u(M(j)).name]
fileid= fopen(filename);
u=fread(fileid,'double');


filename=[path name_x(M(j)).name];
fileid= fopen(filename);
X=fread(fileid,'double');
[XS,ind]=sort(X);

US=u(ind);

clear X u;

for i=1:length(x)-1
    ind=find((XS>=x(i))&(XS<x(i+1)));
    N(i)=length(ind);
    xc(i)=(x(i)+x(i+1))/2;
    n(i)=N(i)/(dx*Ly*Lz);
end

n_p(j,:)=n;


end   
figure;
plot(n_p(16:end,:)');

n_p_mean=mean(n_p(16:end,:),1);
hold on;plot(n_p_mean,'.');
save './part_num_mono_sim5.mat' n_p_mean;

