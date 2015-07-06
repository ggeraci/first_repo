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


name_T=dir([path 'Part_T*']);
name_x=dir([path 'Part_x*']);

T_p_inlet=[];T_p_outlet=[];

for i=1:41%length(name_T)
    TS(i)=str2num(name_x(i).name(8:end-4));
end

[N,M]=sort(TS);
%name_T.name=name_T(ind).name


for j=16:41%length(name_T)
     j 
filename=[path name_T(M(j)).name]
fileid= fopen(filename);
T=fread(fileid,'double');


filename=[path name_x(M(j)).name];
fileid= fopen(filename);
X=fread(fileid,'double');
[XS,ind]=sort(X);

TS=T(ind);

clear X T;

for i=1:length(x)-1
    ind=find((XS>=x(i))&(XS<x(i+1)));
    temp(i)=mean(TS(ind));
    xc(i)=(x(i)+x(i+1))/2;
end

T_p(j,:)=temp;

ind=find((XS>=x(1))&(XS<x(2)));
T_p_inlet=[T_p_inlet; TS(ind)];

ind=find((XS>=x(end-1))&(XS<x(end)));
T_p_outlet=[T_p_outlet; TS(ind)];

%figure
%plot(xc,temp,'LineWidth',2.5); hold on;

end   
figure;
plot(T_p(16:end,:)');

T_p_mean=mean(T_p(16:end,:),1);
hold on;plot(T_p_mean,'.');
save './part_temp_mono_sim3.mat' T_p_mean T_p_inlet T_p_outlet;

