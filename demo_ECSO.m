%% demo_ECSO
%
% This demo illustrates the Entropy-Based Convex Set Optimization (ECSO) Endmember Extraction algorithm
% introduced in the paper 
%
% Shah, Dharambhai, Tanish Zaveri, Yogesh N. Trivedi, and Antonio Plaza. 
% "Entropy-based convex set optimization for spatial–spectral endmember extraction from hyperspectral images." 
% IECSO Journal of Selected Topics in Applied Earth Observations and Remote Sensing
% 13 (2020): 4200-4213.

clc;
close all;
clear all;

%% Image Read
s=load('Cuprite.mat');  % link for data source : https://rslab.ut.ac.ir/data
p=s.nRow;
q=s.nCol;
Bands=188;
Y=s.Y;
x=hyperConvert3d(Y,p,q,Bands);
for i=1:1:Bands
    R=x(:,:,i);
    minr=min(min(R));
    maxr=max(max(R));
    R1=(R-minr)/(maxr-minr);
    R1(R1==0)=0.0001;
    xp(:,:,i)=R1;
end

%% Virtual Dimension
VD=12;

%% ECSO algorithm
[endmemberindex] = ECSO(Y,VD);
endmemberindex_ECSO=change_index(endmemberindex,p,q);

%% VCA algorithm
[U_VCA,e_index,snrEstimate]=hyperVca(Y,VD);
endmemberindex_VCA=change_index(e_index,p,q);    

%% GT 
t1=load('groundTruth_Cuprite_nEnd12.mat');
gt=t1.M;
tit=t1.cood;
n1=gt(3:103,:);
n2=gt(114:147,:);
n3=gt(168:220,:);
gt=[n1;n2;n3];
[gt_m,gt_n]=size(gt);

%% Total Spectral Angle Mapper (TSAM) calculations
for i=1:gt_n
    for j=1:Bands
        extracted_ECSO(j,i)=xp(endmemberindex_ECSO(i,1),endmemberindex_ECSO(i,2),j);
        extracted_VCA(j,i)=xp(endmemberindex_VCA(i,1),endmemberindex_VCA(i,2),j);
    end
end

[ex_m,ex_n]=size(extracted_VCA);
store_ECSO=[0,0];
store_VCA=[0,0];
sam_VCA=0;
sam_ECSO=0;
sam_total_ECSO=0;
sam_total_VCA=0;

for i=1:gt_n
    for j=1:ex_n
            Mat_SAM_VCA(i,j)=real(acos(dot(gt(:,i),extracted_VCA(:,j))/(norm(gt(:,i)*norm(extracted_VCA(:,j))))));
            Mat_SAM_ECSO(i,j)=real(acos(dot(gt(:,i),extracted_ECSO(:,j))/(norm(gt(:,i)*norm(extracted_ECSO(:,j))))));
    end
end

for i=1:gt_n
    %VCA
    [max_value1,mrow]=min(Mat_SAM_VCA);
    [max_value,col_VCA]=min(max_value1);
    sam_total_VCA=sam_total_VCA+max_value;
    sam_VCA=[sam_VCA;max_value];
    row_VCA=mrow(col_VCA);
    s1=[row_VCA,col_VCA];
    store_VCA=[store_VCA;s1];
    save_VCA(row_VCA)=max_value;
    Mat_SAM_VCA(row_VCA,:)=[100*ones];
    Mat_SAM_VCA(:,col_VCA)=[100*ones];
    %ECSO
    [max_value1,mrow]=min(Mat_SAM_ECSO);
    [max_value,col_ECSO]=min(max_value1);
    sam_total_ECSO=sam_total_ECSO+max_value;
    sam_ECSO=[sam_ECSO;max_value];
    row_ECSO=mrow(col_ECSO);
    s1=[row_ECSO,col_ECSO];
    store_ECSO=[store_ECSO;s1];
    save_ECSO(row_ECSO)=max_value;
    Mat_SAM_ECSO(row_ECSO,:)=[100*ones];
    Mat_SAM_ECSO(:,col_ECSO)=[100*ones];
end

disp('TSAM of VCA');
disp(sum(save_VCA));
disp('TSAM of ECSO');
disp(sum(save_ECSO));