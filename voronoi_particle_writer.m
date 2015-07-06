close all;
clear all;

clc;

path=['/hpcc-thecus-mrahmani/test_cases_poly/Monodisperse_sim5/data_IO/'];
path=['./data_IO/'];
nf_f=120000;
nf=[80000:2000:nf_f];
Lx=0.16;

for i=1:length(nf);
i    
filename=[path 'Part_x_' num2str(nf(i)) '.bin'];
fileid= fopen(filename);
px=fread(fileid,'double');
Np=length(px);

ind=find(px>15/16*Lx);

filename=[path 'Part_y_' num2str(nf(i)) '.bin'];
fileid= fopen(filename);
py=fread(fileid,'double');

filename=[path 'Part_z_' num2str(nf(i)) '.bin'];
fileid= fopen(filename);
pz=fread(fileid,'double');

out_file=['Part_all_' num2str(nf(i)) '.txt'];
loc=[px(ind)'; py(ind)'; pz(ind)'];
%loc=[[1:Np]; px(1:end)'; py(1:end)'; pz(1:end)'];

fileID = fopen(out_file,'w');
fprintf(fileID,'%12.8f %12.8f %12.8f\n',loc);
fclose(fileID);
end
