clc
clear
close all

nt=48;nr=36;
idx=1;
for Threshold=[0.4 0.3 0.2 0.1]
N=47;
n_e_val=[];
nr_after_pairing_all=[];
nr_after_pairing_all_RT=[];
paired_UE_list_all=[];
% H_SF=zeros(nr,nt,N);
for k=0:N
%% Reading binary files for recorded CSI from testbed
if(Threshold==0.4)
    file_name = strcat('channel_data/th_04_qpsk/ch',int2str(k),'.bkp');     %Locating the binary file recorded from testbed for CSI matrix (The file 'ch0.bkp' contains the 36×48 CSI matrix stored in array format, along with the paired UE list and the number of paired UEs.)
    r_ch_id=fopen(file_name,'r');                                           %Opening the binary file corresponding to above location
    r_ch_file=fread(r_ch_id,inf,'uint32');                                  %Reading the binary file of CSI matrix in unsigned int format
    fclose(r_ch_id);                                                        %Closing the binary file
    r_ch_bin_arr=de2bi(r_ch_file,32,'left-msb');                            %Converting the decimal datas to binary data
    r_ch_bin_big_endian_arr=little_to_big_endian_arr(r_ch_bin_arr);         %Performing little to big endian conversion
    r_ch_arr=bi_to_sp_float_arr(r_ch_bin_big_endian_arr);                   %Converting big endian binary to single precision float array, where the complex array is stored as interleaved real and imaginary samples.
    r_ch_c_arr=(r_ch_arr(1:2:2*nr*nt)+1i*r_ch_arr(2:2:2*nr*nt));            %Converting the real array to complex array
    H_ch=reshape(r_ch_c_arr,nt,nr).';                                       %Rearranging the elements of the complex array in matrix form to get the actual CSI matrix for kth subframe
    H=H_ch;                                                                 %CSI matrix of kth subframe
    nr_after_pairing_all_RT=[nr_after_pairing_all_RT,r_ch_arr(2*nr*nt+nr+1)]; % Number of paired UEs recorded from the testbed
    H(3:4,:)=zeros(2,nt);                                                   %UEs 3 and 4 did not receive any signal and are considered nulls. Since the recorded data contains the minimum data-format value at these positions, these values are set to zero in post-processing.
    H_SF(:,:,k+1)=H;                                                        %CSI matrix for all subframes (48 subframes recorded in each experiment)
elseif(Threshold==0.3)
    file_name = strcat('channel_data/th_03_qpsk/ch',int2str(k),'.bkp');
    r_ch_id=fopen(file_name,'r');
    r_ch_file=fread(r_ch_id,inf,'uint32');
    fclose(r_ch_id);
    r_ch_bin_arr=de2bi(r_ch_file,32,'left-msb');
    r_ch_bin_big_endian_arr=little_to_big_endian_arr(r_ch_bin_arr);
    r_ch_arr=bi_to_sp_float_arr(r_ch_bin_big_endian_arr);
    r_ch_c_arr=(r_ch_arr(1:2:2*nr*nt)+1i*r_ch_arr(2:2:2*nr*nt));
    H_ch=reshape(r_ch_c_arr,nt,nr).';
    H=H_ch;                         
    nr_after_pairing_all_RT=[nr_after_pairing_all_RT,r_ch_arr(2*nr*nt+nr+1)]; 
    H(27:28,:)=zeros(2,nt);                                                 %UEs 27 and 28 did not receive any signal, hence nulls.
    H_SF(:,:,k+1)=H;                
elseif(Threshold==0.2)
    file_name = strcat('channel_data/th_02_16qam/ch',int2str(k),'.bkp');
    r_ch_id=fopen(file_name,'r');
    r_ch_file=fread(r_ch_id,inf,'uint32');
    fclose(r_ch_id);
    r_ch_bin_arr=de2bi(r_ch_file,32,'left-msb');
    r_ch_bin_big_endian_arr=little_to_big_endian_arr(r_ch_bin_arr);
    r_ch_arr=bi_to_sp_float_arr(r_ch_bin_big_endian_arr);
    r_ch_c_arr=(r_ch_arr(1:2:2*nr*nt)+1i*r_ch_arr(2:2:2*nr*nt));
    H_ch=reshape(r_ch_c_arr,nt,nr).';
    H=H_ch;                         
    nr_after_pairing_all_RT=[nr_after_pairing_all_RT,r_ch_arr(2*nr*nt+nr+1)]; 
    H(31:36,:)=zeros(6,nt);                                                 %UEs 31 to 36 did not receive any signal, hence nulls.
    H_SF(:,:,k+1)=H;                
elseif(Threshold==0.1)
    file_name = strcat('channel_data/th_01_64qam/ch',int2str(k),'.bkp');
    r_ch_id=fopen(file_name,'r');
    r_ch_file=fread(r_ch_id,inf,'uint32');
    fclose(r_ch_id);
    r_ch_bin_arr=de2bi(r_ch_file,32,'left-msb');
    r_ch_bin_big_endian_arr=little_to_big_endian_arr(r_ch_bin_arr);
    r_ch_arr=bi_to_sp_float_arr(r_ch_bin_big_endian_arr);
    nr_after_pairing_all_RT=[nr_after_pairing_all_RT,r_ch_arr(2*nr*nt+nr+1)]; 
    
    file_name = strcat('channel_data/th_01_16qam/ch',int2str(k),'.bkp');
    r_ch_id=fopen(file_name,'r');
    r_ch_file=fread(r_ch_id,inf,'uint32');
    fclose(r_ch_id);
    r_ch_bin_arr=de2bi(r_ch_file,32,'left-msb');
    r_ch_bin_big_endian_arr=little_to_big_endian_arr(r_ch_bin_arr);
    r_ch_arr=bi_to_sp_float_arr(r_ch_bin_big_endian_arr);
    r_ch_c_arr=(r_ch_arr(1:2:2*nr*nt)+1i*r_ch_arr(2:2:2*nr*nt));
    H_ch=reshape(r_ch_c_arr,nt,nr).';
    H=H_ch;                         
    H(13:16,:)=zeros(4,nt);                                                 %UEs 13 to 16 did not receive any signal, hence nulls.
    H(31:36,:)=zeros(6,nt);                                                 %UEs 31 to 36 did not receive any signal, hence nulls.
    H_SF(:,:,k+1)=H;                
end
%% NULL removal
UE_list=zeros(1,nr);
count_row=1;
for row=1:nr
    nz=nnz(H(row,:));
    if(nz~=0)
        UE_list(count_row)=row;
        count_row=count_row+1;
    end
end
nr_after_NULL_check=count_row-1;                                            %Discarding the UEs who did not receive any signal
%% Row Nomalization
for r=1:nr_after_NULL_check
    H_norm(UE_list(r),:)=H(UE_list(r),:)/norm(H(UE_list(r),:));             %Normalizing the channel vector for each UE to remove path loss
end
%% Pairing algorithm
UE1=1;
paired_UE_list=[UE_list(UE1)];
for row=1:nr_after_NULL_check                                               %Finding the spatially separable UEs, whose channel vectors are nearly orthogonal
    if(row~=UE1)
        corr_count=0;
        for index=1:nnz(paired_UE_list)
            corr_metric=abs(H_norm(UE_list(row),:)*H_norm(paired_UE_list(index),:)'); %Correlation check between channel vector of UEs
            disp_corr_metric(row,index)=corr_metric;
            if(corr_metric<Threshold)
                corr_count=corr_count+1;
            end
        end
        if(corr_count==length(paired_UE_list))
            paired_UE_list=[paired_UE_list,UE_list(row)];                   %List of all spatially separable UEs 
        end
    end
end
nr_after_pairing=length(paired_UE_list);                                    %Number of paired UEs
nr_after_pairing_all=[nr_after_pairing_all,nr_after_pairing];               %Number of paired UEs for all subframes
paired_UE_list_all=[paired_UE_list_all,paired_UE_list];                     %Paired UE list for all subframes
%% channel matrix after pairing
H_paired=zeros(nr_after_pairing,nt);
for row=1:nr_after_pairing
    H_paired(row,:)=H_norm(paired_UE_list(row),:);                          %Constructing new CSI matrix for the list of paired UEs
end

W=H_paired'*(H_paired*H_paired');                                           %Generating precoder matrix
diagnoal_kth_channel=abs(diag(H_paired*W)).^2;
s_HW=size(H_paired*W);
n_diag=s_HW(1);

%% Eigen analysis
e_val_kth_channel=eig(H_paired*H_paired');                                  %Finding the eigen values of CSI matrix
p=size(H_paired*H_paired');
ein=p(1);
n_e_val_kth_channel=e_val_kth_channel/max(e_val_kth_channel);               %Computing normalized eigen values
n_e_val=[n_e_val;n_e_val_kth_channel];                                      %Normalized eigen values for 48 subframes

prod(sqrt(e_val_kth_channel));
end
n_e_val_Th{idx}=n_e_val;                                                    %Normalized eigen values for 4 different threshold values across all 48 subframes
nr_after_pairing_all_Th{idx}=nr_after_pairing_all_RT;                       %Number of paired UEs for 4 different threshold values across all 48 subframes
idx=idx+1;
end
% Downlink received constellations are captured for one particular subframe, the paired UEs in that subframe are listed below for different experiments
paired_UE_list_th_04_qpsk=[1,2,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,29,30,31,32,33,34,35,36];
paired_UE_list_th_03_qpsk=[1,2,3,4,5,8,9,10,11,12,15,16,18,19,22,23,25,29,31,32,33,34,35,36];
paired_UE_list_th_02_16qam=[1,3,6,7,12,13,14,15,18,20,23,26];
paired_UE_list_th_01_qpsk=[1,6,13,16,34];
paired_UE_list_th_01_16qam=[1,4,7,11,25];
paired_UE_list_th_01_64qam=[1,4,6,11,15];
%% Plotting Normalized Eigenvalues for different correlation thresholds
Normalized_Eigenvalue_CDF(n_e_val_Th);                                      %Plotting the CDF of normalized eigen values for different correlation thresholds
%% Plotting Sum spectral efficiency for different correlation thresholds
Sum_spectral_efficiency_CDF(nr_after_pairing_all_Th);                       %Plotting the CDF of sum spectral efficiency for different correlation thresholds
%% Downlink received constellation diagrams
REVM_th_04_qpsk=DL_received_constellation(paired_UE_list_th_04_qpsk,4,'qpsk'); %Downlink received constellation diagrams for all paired UEs decoded with QPSK for correlation threshold=0.4
disp('REVM for all UEs (Threshold = 0.4 and QPSK)');                           
disp(["UE", "REVM"]);
disp(REVM_th_04_qpsk);                                                         %Printing the REVM values for all paired UEs
REVM_th_03_qpsk=DL_received_constellation(paired_UE_list_th_03_qpsk,3,'qpsk'); %Downlink received constellation diagrams for all paired UEs decoded with QPSK for correlation threshold=0.3
disp('REVM for all UEs (Threshold = 0.3 and QPSK)');                           
disp(["UE", "REVM"]);
disp(REVM_th_03_qpsk);                                                         %Printing the REVM values for all paired UEs
REVM_th_02_16qam=DL_received_constellation(paired_UE_list_th_02_16qam,2,'16qam'); %Downlink received constellation diagrams for all paired UEs decoded with 16QAM for correlation threshold=0.2
disp('REVM for all UEs (Threshold = 0.2 and 16QAM)');
disp(["UE", "REVM"]);
disp(REVM_th_02_16qam);                                                           %Printing the REVM values for all paired UEs
REVM_th_01_qpsk=DL_received_constellation(paired_UE_list_th_01_qpsk,1,'qpsk'); %Downlink received constellation diagrams for all paired UEs decoded with QPSK for correlation threshold=0.1
disp('REVM for all UEs (Threshold = 0.1 and QPSK)');
disp(["UE", "REVM"]);
disp(REVM_th_01_qpsk);                                                         %Printing the REVM values for all paired UEs
REVM_th_01_16qam=DL_received_constellation(paired_UE_list_th_01_16qam,1,'16qam'); %Downlink received constellation diagrams for all paired UEs decoded with 16QAM for correlation threshold=0.1
disp('REVM for all UEs (Threshold = 0.1 and 16QAM)');
disp(["UE", "REVM"]);
disp(REVM_th_01_16qam);                                                           %Printing the REVM values for all paired UEs
REVM_th_01_64qam=DL_received_constellation(paired_UE_list_th_01_64qam,1,'64qam'); %Downlink received constellation diagrams for all paired UEs decoded with 64QAM for correlation threshold=0.1
disp('REVM for all UEs (Threshold = 0.1 and 64QAM)');
disp(["UE", "REVM"]);
disp(REVM_th_01_64qam);                                                           %Printing the REVM values for all paired UEs
%% Function definition for Plotting Normalized Eigenvalues for different correlation thresholds
function []=Normalized_Eigenvalue_CDF(n_e_val_Th) 

n_e_val_48x36_Th_0_4=n_e_val_Th{1};
n_e_val_48x36_Th_0_3=n_e_val_Th{2};
n_e_val_48x36_Th_0_2=n_e_val_Th{3};
n_e_val_48x36_Th_0_1=n_e_val_Th{4};

set(0,'defaultAxesFontSize',14)
set(0,'defaultAxesFontName','Times')

%New Expt
figure
h1=cdfplot(n_e_val_48x36_Th_0_4);
set(h1,'Marker','o','MarkerIndices',[1:200:3050]) %,'Color','r'
hold on
h2=cdfplot(n_e_val_48x36_Th_0_3);
set(h2,'Marker','s','MarkerIndices',[1:200:2200]) %,'Color','k'
h3=cdfplot(n_e_val_48x36_Th_0_2);
set(h3,'Color','k','Marker','d','MarkerIndices',[1:100:1000]) %,'Color','m'
% hold on
h4=cdfplot(n_e_val_48x36_Th_0_1);
set(h4,'Marker','^','MarkerIndices',[1:25:230]) %,'Color','g'

xlabel('Normalized eigenvalues');
ylabel('CDF');
title('');
legend("Th=0.4","Th=0.3","Th=0.2","Th=0.1",'Location','northwest');
end
%% Function definition for Plotting Sum spectral efficiency for different correlation thresholds
function []=Sum_spectral_efficiency_CDF(nr_after_pairing_all_Th)

% C_RT_48x36_Th_0_4=2*nr_after_pairing_all_Th{1};
C_RT_48x36_Th_0_3=2*nr_after_pairing_all_Th{2};
C_RT_48x36_Th_0_2=4*(nr_after_pairing_all_Th{3}-1);
C_RT_48x36_Th_0_1=6*nr_after_pairing_all_Th{4};
set(0,'defaultAxesFontSize',14)
set(0,'defaultAxesFontName','Times')

figure
% h1=cdfplot(C_RT_48x36_Th_0_4);
% set(h1,'Color','b','Marker','o','MarkerIndices',[1:2:10]) 
% hold on
h2=cdfplot(C_RT_48x36_Th_0_3);
set(h2,'Color','r','Marker','s','MarkerIndices',[1:1:10])
hold on
h3=cdfplot(C_RT_48x36_Th_0_2);
set(h3,'Color','k','Marker','d','MarkerIndices',[1:1:10])
h4=cdfplot(C_RT_48x36_Th_0_1);
set(h4,'Color','m','Marker','x','MarkerIndices',[1:1:10])

xlabel('Sum spectral efficiency (bits/sec/Hz)');
ylabel('CDF');
title('');
legend("Th=0.3","Th=0.2","Th=0.1",'Location','northwest');
axis([0,60,0,1])
end
%% Function definition for Downlink received constellation diagrams
function [REVM_UE]=DL_received_constellation(paired_UE_list,threshold,modulation)
% This script reads downlink constellation data from binary files.
% Prints REVM for each paired UE, and plots received constellations for QPSK, 16-QAM, and 64-QAM under different threshold settings recorded from testbed.
% Inputs: paired_UE_list- list of paired UEs recorded from a experiment
%         threshold- correlation threshold considered in the experiment multipled by 10
%         modulation- modulation scheme considered in the experiment
% Output: Plots constellation diagram for all paired UEs and returns their corresponding REVM values.

figure;
for n=1:length(paired_UE_list)
file_name = strcat('dl_received_constellation/th_0',int2str(threshold),'_',modulation,'/ue',int2str(paired_UE_list(n)),'.bkp');
ue_id=fopen(file_name,'r');
ue_file_arr=fread(ue_id,inf,'uint32');
fclose(ue_id);
ue_bin_arr=de2bi(ue_file_arr,32,'left-msb');
ue_bin_big_endian_arr=little_to_big_endian_arr(ue_bin_arr);
ue=bi_to_sp_float_arr(ue_bin_big_endian_arr);

if(threshold==3 && strcmp(modulation, 'qpsk'))
    REVM(n)=ue(385);  %03_qpsk
else
    REVM(n)=ue(433);    %01_qpsk, 01_16qam, 01_64qam, 02_16qam, 04_qpsk
end

if(strcmp(modulation, '64qam') || strcmp(modulation, '16qam'))
    sym=8;  %Total number of ofdm symbols +1
else
    sym=7;  %Total number of ofdm symbols +1
end

ue_1=(ue(1:(48*(sym-1))));
ue_cc=(ue_1(1:2:(48*(sym-1)))+1i*ue_1(2:2:(48*(sym-1))));
if(strcmp(modulation, '64qam') || strcmp(modulation, '16qam'))
    ue_c=ue_cc(1:end);
else
    ue_c=ue_cc(49:end);
end

if(threshold==1 && strcmp(modulation, '16qam'))
    ue_c=1.8*ue_cc(1:end);
    REVM(n)=REVM_computation(ue_c,paired_UE_list(n));
end

set(0,'defaultAxesFontSize',14)
set(0,'defaultAxesFontName','Times')
xstring='InPhase';
ystring='QPhase';

if(strcmp(modulation, 'qpsk'))
if(n<=12)
    if(threshold==3 || threshold==4)
        subplot(3,4,n);
    elseif(threshold==1)
        subplot(2,3,n);
    end
    scatter(real(ue_c)*sqrt(2),imag(ue_c)*sqrt(2));
    title(['UE',int2str(paired_UE_list(n))]);
    grid on;
    xlabel(xstring);
    ylabel(ystring);
    axis([-2 2 -2 2])
    if(threshold==4)
        sgtitle('QPSK (Th=0.4)');
    elseif(threshold==3)
        sgtitle('QPSK (Th=0.3)');
    elseif(threshold==1)
        sgtitle('QPSK (Th=0.1)');
    end
    box
elseif(n>12 && n<=24)
    if(n==13)
        figure;
    end
    subplot(3,4,n-12);
    scatter(real(ue_c)*sqrt(2),imag(ue_c)*sqrt(2));
    title(['UE',int2str(paired_UE_list(n))]);
    grid on;
    xlabel(xstring);
    ylabel(ystring);
    axis([-2 2 -2 2])
    if(threshold==4)
        sgtitle('QPSK (Th=0.4)');
    elseif(threshold==3)
        sgtitle('QPSK (Th=0.3)');
    elseif(threshold==1)
        sgtitle('QPSK (Th=0.1)');
    end
    box
elseif(n>24 && n<=36)
    if(n==25)
        figure;
    end
    subplot(3,4,n-24);
    scatter(real(ue_c)*sqrt(2),imag(ue_c)*sqrt(2));
    title(['UE',int2str(paired_UE_list(n))]);
    grid on;
    xlabel(xstring);
    ylabel(ystring);
    axis([-2 2 -2 2])
    if(threshold==4)
        sgtitle('QPSK (Th=0.4)');
    elseif(threshold==3)
        sgtitle('QPSK (Th=0.3)');
    elseif(threshold==1)
        sgtitle('QPSK (Th=0.1)');
    end
    box
end
elseif(strcmp(modulation, '16qam'))
    if(threshold==2)
        subplot(3,4,n);
    elseif(threshold==1)
        subplot(2,3,n);
    end
    scatter(real(ue_c)*sqrt(10),imag(ue_c)*sqrt(10));
    title(['UE',int2str(paired_UE_list(n))]);
    grid on;
    xlabel(xstring);
    ylabel(ystring);
    axis([-4 4 -4 4])
    if(threshold==2)
        sgtitle('16QAM (Th=0.2)');
    elseif(threshold==1)
        sgtitle('16QAM (Th=0.1)');
    end
    box
elseif(strcmp(modulation, '64qam'))
    subplot(2,3,n);
    scatter(real(ue_c)*sqrt(42),imag(ue_c)*sqrt(42));
    title(['UE',int2str(paired_UE_list(n))]);
    grid on;
    xlabel(xstring);
    ylabel(ystring);
    axis([-10 10 -10 10])
    sgtitle('64QAM (Th=0.1)');
    box
end
end
REVM_UE=[paired_UE_list.' REVM.'];
end
%% REVM Computation
function [EVM_ue]=REVM_computation(ue_c,UE_num)
sym_num=9;      %number of ofdm symbols for data
sym_num_1=7;    %number of ofdm symbols for which EVM is calculated
%Actual transmitted data is denoted as ap1,ap2,..,ap6, which will be used for REVM calculation
%48x36 9 ofdm 16QAM data symbol
ap1=[0.3162,0.9486,0.3162,0.9486,-0.3162,0.9486,0.9486,-0.3162,0.9486,-0.9486,0.3162,-0.9486,-0.3162,0.9486,0.9486,0.3162,-0.9486,-0.9486,-0.9486,0.3162,-0.9486,-0.3162,0.9486,-0.9486,-0.3162,0.3162,0.3162,0.9486,0.3162,0.3162,0.3162,-0.3162,0.9486,0.9486,0.3162,-0.9486,-0.3162,-0.9486,0.9486,0.9486,0.9486,-0.3162,-0.3162,0.3162,-0.3162,0.9486,0.9486,-0.9486,-0.3162,0.9486,-0.9486,0.9486,-0.9486,-0.9486,0.9486,0.3162,-0.9486,0.3162,0.3162,-0.9486,-0.3162,0.9486,0.3162,0.3162,0.3162,-0.3162,0.3162,-0.9486,0.3162,0.9486,-0.9486,0.9486,-0.9486,-0.3162,0.3162,0.9486,0.3162,0.9486,0.9486,0.9486,0.3162,0.9486,0.3162,-0.9486,-0.9486,0.3162,-0.3162,-0.3162,-0.3162,-0.3162,-0.9486,0.9486,0.9486,-0.9486,0.9486,0.9486,-0.3162,-0.9486,-0.9486,-0.3162,-0.9486,0.9486,-0.9486,-0.9486,0.3162,0.3162,0.3162,0.9486,-0.3162,0.9486,0.9486,0.9486,0.9486,-0.3162,0.3162,-0.3162,0.9486,-0.3162,0.3162,-0.9486,0.3162,0.3162,0.9486,0.3162,-0.9486,-0.9486,0.9486,0.9486,-0.3162,0.9486,0.9486,-0.9486,-0.3162,0.3162,0.9486,-0.3162,0.3162,-0.3162,-0.9486,0.3162,0.9486,-0.9486,-0.3162,0.9486,0.9486,-0.3162,-0.9486,0.9486,0.9486,-0.3162,-0.9486,-0.9486,0.9486,0.3162,-0.9486,-0.3162,0.3162,-0.9486,0.3162,-0.9486,-0.3162,-0.3162,-0.3162,-0.9486,-0.9486,0.3162,0.9486,-0.3162,-0.3162,0.3162,0.3162,0.3162,0.9486,0.3162,-0.3162,0.3162,0.3162,-0.9486,-0.9486,-0.9486,-0.9486,-0.3162,0.9486,-0.3162,-0.3162,0.3162,0.3162,0.3162,-0.9486,-0.3162,0.9486,-0.9486,-0.9486,-0.9486,-0.9486,-0.3162,0.3162,-0.3162,0.3162,-0.9486,-0.3162,0.3162,-0.3162,-0.3162,-0.9486,-0.9486,0.9486,-0.3162,-0.3162,-0.3162,0.3162,-0.3162,-0.3162,0.9486,-0.9486,0.3162,0.3162,-0.3162,-0.9486,-0.3162,0.9486,0.3162,0.3162,-0.9486,-0.9486,-0.9486,-0.3162,0.3162,-0.3162,-0.3162,-0.3162,0.3162,0.3162,0.3162,-0.3162,0.9486,-0.3162,0.3162,0.3162,-0.9486,-0.9486,-0.9486,-0.3162,-0.3162,0.3162,-0.9486,-0.9486,-0.9486,0.9486,0.3162,0.9486,0.9486,-0.9486,0.9486,0.3162,-0.9486,-0.3162,0.9486,0.3162,-0.3162,-0.9486,-0.9486,0.9486,-0.9486,-0.9486,0.9486,-0.3162,0.9486,0.9486,0.9486,0.3162,0.9486,-0.9486,-0.3162,0.3162,0.9486,-0.9486,-0.9486,-0.9486,0.3162,0.9486,-0.3162,0.9486,-0.3162,0.3162,0.3162,0.3162,-0.3162,-0.3162,-0.3162,-0.9486,0.9486,-0.3162,0.9486,-0.9486,0.3162,-0.9486,-0.3162,-0.9486,0.3162,0.3162,-0.9486,-0.3162,-0.3162,0.9486,0.3162,-0.3162,-0.9486,-0.9486,-0.3162,0.3162,-0.3162,-0.9486,-0.9486,0.9486,-0.9486,-0.9486,0.9486,-0.9486,0.3162,-0.3162,0.9486,0.9486,-0.3162,0.9486,0.3162,-0.9486,-0.9486,0.9486,0.9486,0.9486,0.3162,-0.3162,0.9486,-0.3162,-0.9486,0.9486,-0.9486,0.9486,-0.9486,0.9486,0.3162,0.9486,-0.9486,-0.3162,0.9486,0.3162,-0.9486,0.9486,0.9486,0.9486,0.9486,0.9486,0.3162,0.3162,0.9486,0.3162,-0.3162,-0.3162,0.9486,0.9486,0.9486,-0.3162,0.3162,0.9486,-0.3162,-0.9486,0.3162,-0.9486,-0.3162,0.3162,0.9486,0.9486,-0.9486,-0.3162,0.9486,-0.9486,0.3162,-0.9486,0.3162,0.9486,0.3162,0.9486,-0.3162,0.9486,-0.3162,-0.3162,0.3162,0.3162,0.3162,0.9486,-0.3162,0.3162,-0.3162,-0.9486,0.3162,0.9486,-0.3162,-0.3162,-0.9486,-0.3162,-0.3162,0.9486,0.3162,0.9486,0.3162,-0.3162,0.3162,-0.3162,0.9486,-0.9486,0.9486,0.9486,0.9486,0.9486,0.9486,-0.9486,0.9486,-0.3162,0.3162,-0.3162,-0.9486,-0.3162,-0.3162,-0.3162,-0.3162,-0.9486,0.9486,0.3162,0.9486,-0.3162,-0.9486];
ap1_c=ap1(1:2:end)+1i*ap1(2:2:end);
ap1_rearranged = reshape(ap1_c,sym_num,24);
ap1_array = reshape(ap1_rearranged.',1,sym_num*24).';
ap1_Tx = ap1_array(1:sym_num_1*24);

%48x36 9 ofdm 16QAM data symbol
ap2=[0.3162,0.3162,-0.9486,0.3162,0.3162,-0.9486,0.3162,-0.9486,0.3162,-0.3162,0.3162,-0.3162,-0.3162,-0.9486,-0.9486,0.3162,0.9486,-0.3162,-0.3162,-0.3162,0.9486,0.3162,-0.9486,-0.3162,-0.3162,-0.3162,0.9486,-0.3162,0.9486,-0.3162,0.9486,-0.9486,-0.9486,0.9486,0.9486,-0.3162,-0.3162,-0.9486,-0.9486,-0.3162,-0.9486,0.9486,0.9486,0.3162,0.3162,-0.3162,0.9486,-0.9486,0.3162,-0.9486,-0.9486,0.9486,-0.3162,0.9486,0.9486,0.3162,-0.9486,-0.3162,-0.9486,0.9486,0.9486,-0.3162,0.3162,-0.3162,0.9486,-0.3162,-0.3162,-0.3162,0.3162,-0.9486,0.9486,0.9486,0.3162,-0.9486,-0.9486,-0.9486,0.9486,-0.9486,0.3162,0.3162,0.9486,0.9486,0.9486,-0.9486,-0.3162,-0.9486,0.9486,0.3162,0.9486,0.9486,-0.9486,-0.9486,0.9486,-0.3162,0.3162,-0.9486,0.3162,-0.9486,-0.3162,0.9486,0.9486,0.3162,-0.3162,-0.9486,-0.9486,-0.9486,-0.3162,0.3162,0.9486,0.9486,-0.3162,-0.3162,0.3162,-0.9486,-0.3162,-0.9486,0.9486,0.9486,0.3162,-0.3162,0.3162,0.9486,0.3162,-0.3162,-0.3162,0.9486,-0.9486,-0.3162,-0.9486,0.3162,0.9486,-0.3162,0.3162,-0.3162,-0.9486,-0.9486,-0.9486,-0.9486,-0.9486,0.3162,0.9486,0.3162,0.9486,-0.3162,0.3162,0.9486,-0.9486,-0.3162,-0.9486,0.9486,0.3162,0.9486,0.3162,0.3162,0.3162,-0.9486,0.3162,-0.3162,-0.9486,-0.3162,-0.3162,0.3162,0.3162,0.3162,0.9486,0.9486,0.9486,-0.9486,0.9486,0.9486,0.9486,0.9486,-0.9486,0.9486,0.3162,0.9486,0.9486,0.9486,0.9486,0.3162,-0.9486,0.9486,-0.9486,0.9486,-0.9486,0.9486,0.3162,-0.9486,0.3162,0.9486,0.9486,0.9486,0.9486,-0.3162,-0.9486,0.3162,0.3162,-0.9486,-0.3162,0.9486,0.3162,-0.3162,0.3162,-0.3162,-0.9486,-0.3162,0.3162,-0.3162,-0.3162,-0.9486,0.3162,-0.3162,-0.3162,0.9486,0.9486,0.3162,0.3162,0.3162,0.9486,-0.3162,0.9486,-0.9486,-0.9486,0.3162,0.9486,0.3162,-0.3162,-0.9486,-0.3162,-0.3162,-0.3162,-0.9486,0.3162,-0.9486,-0.3162,0.3162,-0.9486,-0.9486,0.9486,-0.3162,-0.3162,0.3162,-0.9486,-0.3162,0.3162,0.9486,0.9486,0.3162,0.9486,-0.9486,-0.9486,-0.3162,0.9486,-0.9486,-0.3162,0.3162,0.3162,-0.9486,-0.3162,-0.9486,0.9486,0.9486,0.3162,-0.9486,0.3162,0.9486,-0.9486,-0.9486,0.3162,-0.9486,0.9486,-0.3162,-0.3162,0.9486,-0.3162,-0.3162,0.9486,0.9486,0.9486,-0.9486,-0.3162,-0.3162,-0.9486,0.9486,0.3162,-0.3162,-0.3162,0.3162,0.9486,-0.3162,-0.9486,-0.9486,0.9486,-0.9486,-0.9486,-0.3162,0.3162,0.3162,-0.9486,-0.9486,0.3162,0.3162,0.3162,0.3162,0.9486,0.3162,-0.3162,0.9486,-0.9486,-0.9486,-0.3162,-0.3162,-0.9486,0.3162,0.9486,0.9486,-0.9486,-0.9486,-0.3162,0.3162,-0.3162,-0.9486,0.9486,0.9486,-0.9486,-0.9486,0.9486,0.3162,0.9486,0.9486,-0.9486,-0.9486,0.3162,-0.3162,0.9486,-0.9486,0.9486,-0.3162,0.9486,-0.3162,0.9486,0.3162,-0.3162,-0.9486,-0.9486,-0.3162,0.3162,0.3162,0.3162,-0.9486,-0.3162,0.3162,-0.9486,0.3162,0.9486,-0.9486,0.9486,0.9486,0.3162,0.9486,0.3162,0.3162,0.3162,0.3162,-0.9486,0.3162,0.9486,-0.3162,0.9486,0.9486,-0.9486,0.9486,0.3162,-0.9486,-0.3162,-0.9486,0.9486,-0.9486,0.3162,0.9486,-0.9486,0.3162,0.9486,0.9486,-0.3162,0.3162,0.3162,-0.3162,-0.9486,0.9486,-0.9486,0.9486,-0.9486,0.3162,-0.3162,-0.9486,0.3162,-0.9486,-0.3162,-0.3162,-0.3162,-0.3162,-0.3162,-0.9486,-0.3162,-0.3162,0.3162,-0.9486,0.9486,0.9486,-0.3162,0.3162,-0.9486,-0.9486,-0.3162,-0.3162,-0.3162,-0.3162,0.3162,0.9486,0.9486,0.3162,0.3162,0.3162,-0.9486,0.9486,-0.9486,-0.9486,0.3162,-0.3162,0.3162,0.9486];
ap2_c=ap2(1:2:end)+1i*ap2(2:2:end);
ap2_rearranged = reshape(ap2_c,sym_num,24);
ap2_array = reshape(ap2_rearranged.',1,sym_num*24).';
ap2_Tx = ap2_array(1:sym_num_1*24);

%48x36 9 ofdm 16QAM data symbol
ap3=[0.3162,0.9486,0.9486,-0.9486,0.9486,-0.3162,0.3162,-0.9486,0.9486,0.3162,0.3162,-0.9486,-0.3162,0.9486,-0.3162,-0.3162,0.3162,0.9486,-0.3162,-0.3162,-0.9486,0.3162,-0.3162,-0.3162,-0.3162,-0.9486,0.9486,0.3162,0.9486,-0.9486,-0.9486,0.9486,0.9486,0.3162,0.3162,0.3162,-0.9486,-0.3162,-0.9486,-0.3162,0.9486,0.3162,-0.3162,-0.9486,0.3162,-0.9486,-0.9486,0.9486,-0.3162,-0.3162,0.9486,-0.3162,-0.3162,-0.3162,0.3162,0.3162,0.9486,0.9486,-0.9486,0.9486,-0.9486,-0.3162,-0.9486,-0.3162,-0.3162,-0.9486,-0.3162,0.3162,-0.3162,-0.9486,-0.9486,0.9486,-0.3162,0.9486,-0.3162,0.3162,0.9486,0.3162,-0.9486,0.9486,0.3162,-0.9486,-0.9486,-0.3162,-0.9486,0.9486,-0.9486,-0.3162,-0.9486,0.9486,-0.9486,0.9486,0.3162,0.9486,-0.3162,-0.9486,0.9486,0.9486,0.3162,0.3162,0.3162,0.3162,-0.3162,-0.9486,-0.9486,0.9486,0.3162,-0.3162,-0.9486,0.9486,-0.9486,-0.3162,-0.9486,-0.9486,0.3162,-0.9486,-0.3162,0.3162,-0.9486,-0.3162,0.9486,-0.3162,-0.3162,-0.9486,0.3162,0.3162,-0.3162,-0.9486,-0.3162,0.3162,-0.3162,-0.9486,0.3162,-0.9486,-0.3162,-0.3162,0.9486,-0.9486,0.9486,-0.3162,0.3162,-0.9486,0.3162,-0.3162,0.3162,-0.3162,0.3162,-0.9486,-0.9486,0.3162,0.9486,-0.9486,-0.3162,-0.3162,-0.9486,0.9486,-0.3162,-0.9486,0.3162,-0.3162,0.9486,0.9486,0.9486,-0.3162,-0.3162,0.3162,-0.3162,-0.9486,-0.9486,0.3162,0.3162,0.3162,0.9486,-0.3162,-0.9486,0.9486,0.9486,-0.9486,-0.3162,0.3162,0.3162,-0.3162,-0.3162,-0.3162,0.9486,0.3162,-0.9486,-0.3162,0.3162,0.3162,0.9486,-0.9486,-0.3162,0.3162,-0.3162,-0.9486,-0.9486,-0.3162,-0.9486,0.9486,0.3162,0.9486,-0.9486,0.9486,-0.9486,-0.3162,-0.9486,-0.9486,0.3162,0.3162,-0.3162,-0.9486,0.9486,-0.9486,-0.3162,0.3162,0.3162,0.9486,0.9486,-0.9486,0.3162,-0.3162,-0.3162,-0.9486,-0.3162,0.3162,0.9486,-0.9486,-0.9486,-0.3162,0.9486,0.9486,-0.9486,-0.9486,-0.9486,0.9486,0.9486,-0.9486,-0.9486,-0.9486,-0.3162,-0.9486,0.3162,-0.9486,0.9486,0.9486,-0.9486,-0.9486,-0.3162,-0.3162,-0.9486,0.9486,0.9486,0.9486,-0.3162,-0.3162,0.3162,0.9486,-0.3162,-0.9486,0.9486,-0.9486,0.9486,0.3162,-0.9486,-0.3162,0.9486,-0.9486,-0.9486,0.9486,0.9486,-0.9486,0.9486,-0.3162,-0.9486,-0.3162,0.3162,-0.3162,0.3162,-0.9486,-0.3162,-0.3162,-0.9486,-0.9486,-0.3162,0.3162,0.3162,-0.3162,0.3162,-0.3162,0.3162,0.9486,-0.3162,0.3162,-0.9486,-0.3162,0.9486,0.3162,-0.9486,0.3162,0.9486,0.9486,0.9486,-0.3162,-0.3162,0.9486,-0.9486,-0.3162,-0.3162,-0.9486,-0.3162,0.9486,0.3162,0.9486,0.9486,-0.9486,0.9486,0.9486,-0.9486,-0.9486,-0.9486,-0.9486,-0.9486,-0.3162,-0.3162,0.9486,0.3162,0.3162,-0.3162,0.3162,-0.3162,0.9486,0.9486,0.9486,-0.9486,-0.9486,0.3162,-0.9486,0.9486,0.3162,-0.3162,0.9486,-0.3162,-0.9486,0.3162,0.9486,-0.9486,-0.3162,0.9486,0.3162,0.3162,0.3162,0.9486,-0.9486,0.3162,-0.9486,0.3162,-0.9486,0.9486,-0.3162,0.3162,0.3162,0.9486,0.3162,0.3162,0.3162,0.3162,-0.9486,0.9486,-0.9486,-0.9486,0.3162,-0.3162,-0.3162,0.3162,-0.3162,0.9486,0.9486,0.3162,0.3162,0.3162,0.9486,-0.3162,0.9486,-0.3162,-0.3162,-0.9486,0.9486,-0.9486,0.3162,-0.3162,-0.3162,0.9486,-0.3162,0.9486,0.3162,0.3162,0.3162,-0.3162,0.3162,0.3162,0.3162,-0.3162,-0.9486,0.3162,0.9486,-0.3162,0.3162,0.9486,0.9486,-0.9486,-0.9486,-0.3162,-0.9486,0.3162,0.3162,0.3162,0.9486,-0.9486,0.3162,0.9486,-0.9486,-0.3162,0.3162,0.3162,-0.9486,0.9486,0.9486,-0.3162,-0.3162,-0.9486,0.9486];
ap3_c=ap3(1:2:end)+1i*ap3(2:2:end);
ap3_rearranged = reshape(ap3_c,sym_num,24);
ap3_array = reshape(ap3_rearranged.',1,sym_num*24).';
ap3_Tx = ap3_array(1:sym_num_1*24);

%48x36 9 ofdm 16QAM data symbol
ap4=[0.3162,0.9486,0.3162,0.9486,-0.3162,0.9486,0.9486,-0.3162,0.9486,-0.9486,0.3162,-0.9486,-0.3162,0.9486,0.9486,0.3162,-0.9486,-0.9486,-0.9486,0.3162,-0.9486,-0.3162,0.9486,-0.9486,-0.3162,0.3162,0.3162,0.9486,0.3162,0.3162,0.3162,-0.3162,0.9486,0.9486,0.3162,-0.9486,-0.3162,-0.9486,0.9486,0.9486,0.9486,-0.3162,-0.3162,0.3162,-0.3162,0.9486,0.9486,-0.9486,-0.3162,0.9486,-0.9486,0.9486,-0.9486,-0.9486,0.9486,0.3162,-0.9486,0.3162,0.3162,-0.9486,-0.3162,0.9486,0.3162,0.3162,0.3162,-0.3162,0.3162,-0.9486,0.3162,0.9486,-0.9486,0.9486,-0.9486,-0.3162,0.3162,0.9486,0.3162,0.9486,0.9486,0.9486,0.3162,0.9486,0.3162,-0.9486,-0.9486,0.3162,-0.3162,-0.3162,-0.3162,-0.3162,-0.9486,0.9486,0.9486,-0.9486,0.9486,0.9486,-0.3162,-0.9486,-0.9486,-0.3162,-0.9486,0.9486,-0.9486,-0.9486,0.3162,0.3162,0.3162,0.9486,-0.3162,0.9486,0.9486,0.9486,0.9486,-0.3162,0.3162,-0.3162,0.9486,-0.3162,0.3162,-0.9486,0.3162,0.3162,0.9486,0.3162,-0.9486,-0.9486,0.9486,0.9486,-0.3162,0.9486,0.9486,-0.9486,-0.3162,0.3162,0.9486,-0.3162,0.3162,-0.3162,-0.9486,0.3162,0.9486,-0.9486,-0.3162,0.9486,0.9486,-0.3162,-0.9486,0.9486,0.9486,-0.3162,-0.9486,-0.9486,0.9486,0.3162,-0.9486,-0.3162,0.3162,-0.9486,0.3162,-0.9486,-0.3162,-0.3162,-0.3162,-0.9486,-0.9486,0.3162,0.9486,-0.3162,-0.3162,0.3162,0.3162,0.3162,0.9486,0.3162,-0.3162,0.3162,0.3162,-0.9486,-0.9486,-0.9486,-0.9486,-0.3162,0.9486,-0.3162,-0.3162,0.3162,0.3162,0.3162,-0.9486,-0.3162,0.9486,-0.9486,-0.9486,-0.9486,-0.9486,-0.3162,0.3162,-0.3162,0.3162,-0.9486,-0.3162,0.3162,-0.3162,-0.3162,-0.9486,-0.9486,0.9486,-0.3162,-0.3162,-0.3162,0.3162,-0.3162,-0.3162,0.9486,-0.9486,0.3162,0.3162,-0.3162,-0.9486,-0.3162,0.9486,0.3162,0.3162,-0.9486,-0.9486,-0.9486,-0.3162,0.3162,-0.3162,-0.3162,-0.3162,0.3162,0.3162,0.3162,-0.3162,0.9486,-0.3162,0.3162,0.3162,-0.9486,-0.9486,-0.9486,-0.3162,-0.3162,0.3162,-0.9486,-0.9486,-0.9486,0.9486,0.3162,0.9486,0.9486,-0.9486,0.9486,0.3162,-0.9486,-0.3162,0.9486,0.3162,-0.3162,-0.9486,-0.9486,0.9486,-0.9486,-0.9486,0.9486,-0.3162,0.9486,0.9486,0.9486,0.3162,0.9486,-0.9486,-0.3162,0.3162,0.9486,-0.9486,-0.9486,-0.9486,0.3162,0.9486,-0.3162,0.9486,-0.3162,0.3162,0.3162,0.3162,-0.3162,-0.3162,-0.3162,-0.9486,0.9486,-0.3162,0.9486,-0.9486,0.3162,-0.9486,-0.3162,-0.9486,0.3162,0.3162,-0.9486,-0.3162,-0.3162,0.9486,0.3162,-0.3162,-0.9486,-0.9486,-0.3162,0.3162,-0.3162,-0.9486,-0.9486,0.9486,-0.9486,-0.9486,0.9486,-0.9486,0.3162,-0.3162,0.9486,0.9486,-0.3162,0.9486,0.3162,-0.9486,-0.9486,0.9486,0.9486,0.9486,0.3162,-0.3162,0.9486,-0.3162,-0.9486,0.9486,-0.9486,0.9486,-0.9486,0.9486,0.3162,0.9486,-0.9486,-0.3162,0.9486,0.3162,-0.9486,0.9486,0.9486,0.9486,0.9486,0.9486,0.3162,0.3162,0.9486,0.3162,-0.3162,-0.3162,0.9486,0.9486,0.9486,-0.3162,0.3162,0.9486,-0.3162,-0.9486,0.3162,-0.9486,-0.3162,0.3162,0.9486,0.9486,-0.9486,-0.3162,0.9486,-0.9486,0.3162,-0.9486,0.3162,0.9486,0.3162,0.9486,-0.3162,0.9486,-0.3162,-0.3162,0.3162,0.3162,0.3162,0.9486,-0.3162,0.3162,-0.3162,-0.9486,0.3162,0.9486,-0.3162,-0.3162,-0.9486,-0.3162,-0.3162,0.9486,0.3162,0.9486,0.3162,-0.3162,0.3162,-0.3162,0.9486,-0.9486,0.9486,0.9486,0.9486,0.9486,0.9486,-0.9486,0.9486,-0.3162,0.3162,-0.3162,-0.9486,-0.3162,-0.3162,-0.3162,-0.3162,-0.9486,0.9486,0.3162,0.9486,-0.3162,-0.9486];
ap4_c=ap4(1:2:end)+1i*ap4(2:2:end);
ap4_rearranged = reshape(ap4_c,sym_num,24);
ap4_array = reshape(ap4_rearranged.',1,sym_num*24).';
ap4_Tx = ap4_array(1:sym_num_1*24);

%48x36 9 ofdm 16QAM data symbol
ap5=[0.3162,0.3162,-0.9486,0.3162,0.3162,-0.9486,0.3162,-0.9486,0.3162,-0.3162,0.3162,-0.3162,-0.3162,-0.9486,-0.9486,0.3162,0.9486,-0.3162,-0.3162,-0.3162,0.9486,0.3162,-0.9486,-0.3162,-0.3162,-0.3162,0.9486,-0.3162,0.9486,-0.3162,0.9486,-0.9486,-0.9486,0.9486,0.9486,-0.3162,-0.3162,-0.9486,-0.9486,-0.3162,-0.9486,0.9486,0.9486,0.3162,0.3162,-0.3162,0.9486,-0.9486,0.3162,-0.9486,-0.9486,0.9486,-0.3162,0.9486,0.9486,0.3162,-0.9486,-0.3162,-0.9486,0.9486,0.9486,-0.3162,0.3162,-0.3162,0.9486,-0.3162,-0.3162,-0.3162,0.3162,-0.9486,0.9486,0.9486,0.3162,-0.9486,-0.9486,-0.9486,0.9486,-0.9486,0.3162,0.3162,0.9486,0.9486,0.9486,-0.9486,-0.3162,-0.9486,0.9486,0.3162,0.9486,0.9486,-0.9486,-0.9486,0.9486,-0.3162,0.3162,-0.9486,0.3162,-0.9486,-0.3162,0.9486,0.9486,0.3162,-0.3162,-0.9486,-0.9486,-0.9486,-0.3162,0.3162,0.9486,0.9486,-0.3162,-0.3162,0.3162,-0.9486,-0.3162,-0.9486,0.9486,0.9486,0.3162,-0.3162,0.3162,0.9486,0.3162,-0.3162,-0.3162,0.9486,-0.9486,-0.3162,-0.9486,0.3162,0.9486,-0.3162,0.3162,-0.3162,-0.9486,-0.9486,-0.9486,-0.9486,-0.9486,0.3162,0.9486,0.3162,0.9486,-0.3162,0.3162,0.9486,-0.9486,-0.3162,-0.9486,0.9486,0.3162,0.9486,0.3162,0.3162,0.3162,-0.9486,0.3162,-0.3162,-0.9486,-0.3162,-0.3162,0.3162,0.3162,0.3162,0.9486,0.9486,0.9486,-0.9486,0.9486,0.9486,0.9486,0.9486,-0.9486,0.9486,0.3162,0.9486,0.9486,0.9486,0.9486,0.3162,-0.9486,0.9486,-0.9486,0.9486,-0.9486,0.9486,0.3162,-0.9486,0.3162,0.9486,0.9486,0.9486,0.9486,-0.3162,-0.9486,0.3162,0.3162,-0.9486,-0.3162,0.9486,0.3162,-0.3162,0.3162,-0.3162,-0.9486,-0.3162,0.3162,-0.3162,-0.3162,-0.9486,0.3162,-0.3162,-0.3162,0.9486,0.9486,0.3162,0.3162,0.3162,0.9486,-0.3162,0.9486,-0.9486,-0.9486,0.3162,0.9486,0.3162,-0.3162,-0.9486,-0.3162,-0.3162,-0.3162,-0.9486,0.3162,-0.9486,-0.3162,0.3162,-0.9486,-0.9486,0.9486,-0.3162,-0.3162,0.3162,-0.9486,-0.3162,0.3162,0.9486,0.9486,0.3162,0.9486,-0.9486,-0.9486,-0.3162,0.9486,-0.9486,-0.3162,0.3162,0.3162,-0.9486,-0.3162,-0.9486,0.9486,0.9486,0.3162,-0.9486,0.3162,0.9486,-0.9486,-0.9486,0.3162,-0.9486,0.9486,-0.3162,-0.3162,0.9486,-0.3162,-0.3162,0.9486,0.9486,0.9486,-0.9486,-0.3162,-0.3162,-0.9486,0.9486,0.3162,-0.3162,-0.3162,0.3162,0.9486,-0.3162,-0.9486,-0.9486,0.9486,-0.9486,-0.9486,-0.3162,0.3162,0.3162,-0.9486,-0.9486,0.3162,0.3162,0.3162,0.3162,0.9486,0.3162,-0.3162,0.9486,-0.9486,-0.9486,-0.3162,-0.3162,-0.9486,0.3162,0.9486,0.9486,-0.9486,-0.9486,-0.3162,0.3162,-0.3162,-0.9486,0.9486,0.9486,-0.9486,-0.9486,0.9486,0.3162,0.9486,0.9486,-0.9486,-0.9486,0.3162,-0.3162,0.9486,-0.9486,0.9486,-0.3162,0.9486,-0.3162,0.9486,0.3162,-0.3162,-0.9486,-0.9486,-0.3162,0.3162,0.3162,0.3162,-0.9486,-0.3162,0.3162,-0.9486,0.3162,0.9486,-0.9486,0.9486,0.9486,0.3162,0.9486,0.3162,0.3162,0.3162,0.3162,-0.9486,0.3162,0.9486,-0.3162,0.9486,0.9486,-0.9486,0.9486,0.3162,-0.9486,-0.3162,-0.9486,0.9486,-0.9486,0.3162,0.9486,-0.9486,0.3162,0.9486,0.9486,-0.3162,0.3162,0.3162,-0.3162,-0.9486,0.9486,-0.9486,0.9486,-0.9486,0.3162,-0.3162,-0.9486,0.3162,-0.9486,-0.3162,-0.3162,-0.3162,-0.3162,-0.3162,-0.9486,-0.3162,-0.3162,0.3162,-0.9486,0.9486,0.9486,-0.3162,0.3162,-0.9486,-0.9486,-0.3162,-0.3162,-0.3162,-0.3162,0.3162,0.9486,0.9486,0.3162,0.3162,0.3162,-0.9486,0.9486,-0.9486,-0.9486,0.3162,-0.3162,0.3162,0.9486];
ap5_c=ap5(1:2:end)+1i*ap5(2:2:end);
ap5_rearranged = reshape(ap5_c,sym_num,24);
ap5_array = reshape(ap5_rearranged.',1,sym_num*24).';
ap5_Tx = ap5_array(1:sym_num_1*24);

%48x36 9 ofdm 16QAM data symbol
ap6=[0.3162,0.9486,0.9486,-0.9486,0.9486,-0.3162,0.3162,-0.9486,0.9486,0.3162,0.3162,-0.9486,-0.3162,0.9486,-0.3162,-0.3162,0.3162,0.9486,-0.3162,-0.3162,-0.9486,0.3162,-0.3162,-0.3162,-0.3162,-0.9486,0.9486,0.3162,0.9486,-0.9486,-0.9486,0.9486,0.9486,0.3162,0.3162,0.3162,-0.9486,-0.3162,-0.9486,-0.3162,0.9486,0.3162,-0.3162,-0.9486,0.3162,-0.9486,-0.9486,0.9486,-0.3162,-0.3162,0.9486,-0.3162,-0.3162,-0.3162,0.3162,0.3162,0.9486,0.9486,-0.9486,0.9486,-0.9486,-0.3162,-0.9486,-0.3162,-0.3162,-0.9486,-0.3162,0.3162,-0.3162,-0.9486,-0.9486,0.9486,-0.3162,0.9486,-0.3162,0.3162,0.9486,0.3162,-0.9486,0.9486,0.3162,-0.9486,-0.9486,-0.3162,-0.9486,0.9486,-0.9486,-0.3162,-0.9486,0.9486,-0.9486,0.9486,0.3162,0.9486,-0.3162,-0.9486,0.9486,0.9486,0.3162,0.3162,0.3162,0.3162,-0.3162,-0.9486,-0.9486,0.9486,0.3162,-0.3162,-0.9486,0.9486,-0.9486,-0.3162,-0.9486,-0.9486,0.3162,-0.9486,-0.3162,0.3162,-0.9486,-0.3162,0.9486,-0.3162,-0.3162,-0.9486,0.3162,0.3162,-0.3162,-0.9486,-0.3162,0.3162,-0.3162,-0.9486,0.3162,-0.9486,-0.3162,-0.3162,0.9486,-0.9486,0.9486,-0.3162,0.3162,-0.9486,0.3162,-0.3162,0.3162,-0.3162,0.3162,-0.9486,-0.9486,0.3162,0.9486,-0.9486,-0.3162,-0.3162,-0.9486,0.9486,-0.3162,-0.9486,0.3162,-0.3162,0.9486,0.9486,0.9486,-0.3162,-0.3162,0.3162,-0.3162,-0.9486,-0.9486,0.3162,0.3162,0.3162,0.9486,-0.3162,-0.9486,0.9486,0.9486,-0.9486,-0.3162,0.3162,0.3162,-0.3162,-0.3162,-0.3162,0.9486,0.3162,-0.9486,-0.3162,0.3162,0.3162,0.9486,-0.9486,-0.3162,0.3162,-0.3162,-0.9486,-0.9486,-0.3162,-0.9486,0.9486,0.3162,0.9486,-0.9486,0.9486,-0.9486,-0.3162,-0.9486,-0.9486,0.3162,0.3162,-0.3162,-0.9486,0.9486,-0.9486,-0.3162,0.3162,0.3162,0.9486,0.9486,-0.9486,0.3162,-0.3162,-0.3162,-0.9486,-0.3162,0.3162,0.9486,-0.9486,-0.9486,-0.3162,0.9486,0.9486,-0.9486,-0.9486,-0.9486,0.9486,0.9486,-0.9486,-0.9486,-0.9486,-0.3162,-0.9486,0.3162,-0.9486,0.9486,0.9486,-0.9486,-0.9486,-0.3162,-0.3162,-0.9486,0.9486,0.9486,0.9486,-0.3162,-0.3162,0.3162,0.9486,-0.3162,-0.9486,0.9486,-0.9486,0.9486,0.3162,-0.9486,-0.3162,0.9486,-0.9486,-0.9486,0.9486,0.9486,-0.9486,0.9486,-0.3162,-0.9486,-0.3162,0.3162,-0.3162,0.3162,-0.9486,-0.3162,-0.3162,-0.9486,-0.9486,-0.3162,0.3162,0.3162,-0.3162,0.3162,-0.3162,0.3162,0.9486,-0.3162,0.3162,-0.9486,-0.3162,0.9486,0.3162,-0.9486,0.3162,0.9486,0.9486,0.9486,-0.3162,-0.3162,0.9486,-0.9486,-0.3162,-0.3162,-0.9486,-0.3162,0.9486,0.3162,0.9486,0.9486,-0.9486,0.9486,0.9486,-0.9486,-0.9486,-0.9486,-0.9486,-0.9486,-0.3162,-0.3162,0.9486,0.3162,0.3162,-0.3162,0.3162,-0.3162,0.9486,0.9486,0.9486,-0.9486,-0.9486,0.3162,-0.9486,0.9486,0.3162,-0.3162,0.9486,-0.3162,-0.9486,0.3162,0.9486,-0.9486,-0.3162,0.9486,0.3162,0.3162,0.3162,0.9486,-0.9486,0.3162,-0.9486,0.3162,-0.9486,0.9486,-0.3162,0.3162,0.3162,0.9486,0.3162,0.3162,0.3162,0.3162,-0.9486,0.9486,-0.9486,-0.9486,0.3162,-0.3162,-0.3162,0.3162,-0.3162,0.9486,0.9486,0.3162,0.3162,0.3162,0.9486,-0.3162,0.9486,-0.3162,-0.3162,-0.9486,0.9486,-0.9486,0.3162,-0.3162,-0.3162,0.9486,-0.3162,0.9486,0.3162,0.3162,0.3162,-0.3162,0.3162,0.3162,0.3162,-0.3162,-0.9486,0.3162,0.9486,-0.3162,0.3162,0.9486,0.9486,-0.9486,-0.9486,-0.3162,-0.9486,0.3162,0.3162,0.3162,0.9486,-0.9486,0.3162,0.9486,-0.9486,-0.3162,0.3162,0.3162,-0.9486,0.9486,0.9486,-0.3162,-0.3162,-0.9486,0.9486];
ap6_c=ap6(1:2:end)+1i*ap6(2:2:end);
ap6_rearranged = reshape(ap6_c,sym_num,24);
ap6_array = reshape(ap6_rearranged.',1,sym_num*24).';
ap6_Tx = ap6_array(1:sym_num_1*24);

ap_Tx(:,1)=ap1_Tx;
ap_Tx(:,2)=ap2_Tx;
ap_Tx(:,3)=ap3_Tx;
ap_Tx(:,4)=ap4_Tx;
ap_Tx(:,5)=ap5_Tx;
ap_Tx(:,6)=ap6_Tx;

if(UE_num>=1 && UE_num<=6)
    ap_index=0;
elseif(UE_num>=7 && UE_num<=12)
    ap_index=1;
elseif(UE_num>=13 && UE_num<=18)
    ap_index=2;
elseif(UE_num>=19 && UE_num<=24)
    ap_index=3;
elseif(UE_num>=25 && UE_num<=30)
    ap_index=4;
end
EVM_ue = sqrt(sum(abs(ue_c-ap_Tx(:,UE_num-6*ap_index)).^2)/length(ap1_Tx))*100;
end
%% Little to big endian conversion for array of 32 bit binary datas
function [bi_big_endian_mat]=little_to_big_endian_arr(bi_little_endian_mat)
%For 32bit format
bi_big_endian_mat(:,1:8)=bi_little_endian_mat(:,end-7:end);
bi_big_endian_mat(:,9:16)=bi_little_endian_mat(:,end-15:end-8);
bi_big_endian_mat(:,17:24)=bi_little_endian_mat(:,end-23:end-16);
bi_big_endian_mat(:,25:32)=bi_little_endian_mat(:,end-31:end-24);
end
%% binary matrix to single presion (32 bit) floating point array conversion as per IEEE 754
function [data_out_arr]=bi_to_sp_float_arr(bin_big_endian_mat)
%input : bin_big_endian_mat = binary matrix in big endian mode, 
%        row=index of 32 bit data , column = index of binary of each data
%output : data_out_arr = single presion (32 bit) floating point data array
S=bin_big_endian_mat(:,1);
E=bi2de(bin_big_endian_mat(:,2:9),'left-msb');
F=1*ones(size(bin_big_endian_mat,1),1);
for k=1:23
    F=F+bin_big_endian_mat(:,k+9)*(2^-k);
end
data_out_arr=((-1).^S).*F.*2.^(E-127);
end