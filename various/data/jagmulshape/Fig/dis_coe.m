clear; close all; clc;
addpath /home/zsh/documents/matlab/Toolbox/CurveLab-2.1.2/fdct_wrapping_matlab

filename='../synth12scale/blcrg1.dat';
fid=fopen(filename,'r');
blend=fread(fid,[470 100],'float');
fclose (fid);

filename='../synth12scale/crg1.dat';
fid=fopen(filename,'r');
data=fread(fid,[470 100],'float');
fclose (fid);


C1 = fdct_wrapping(blend,0,1,5,16);
C2 = fdct_wrapping(data,0,1,5,16);

for i=1:length(C1)
    k=1
    for j=1:length(C1{i})
        [m n] = size(C1{i}{j});
        for ii=1:m
            for jj=1:n;
                coe1(i,k)=real(C1{i}{j}(ii,jj));
                coe2(i,k)=real(C2{i}{j}(ii,jj));
                k=k+1;
            end
        end
    end
end

for i=1:5
temp1(i,:)= sort(coe1(i,:),2,'descend');
temp2(i,:)= sort(coe2(i,:),2,'descend');
end

figure;
k=3
plot (temp1(k,1:10000),'b.'); hold on
plot (temp2(k,1:10000),'r'); hold on



   