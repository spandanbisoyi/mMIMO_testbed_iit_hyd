clc
clear
close all

nt=48;nr=36;
idx=1;
for Threshold=0.3
N=23;
n_e_val=[];
nr_after_pairing_all=[];
nr_after_pairing_all_RT=[];
paired_UE_list_all=[];
% H_SF=zeros(nr,nt,N);
for k=0:N
%% Reading binary files for recorded CSI from testbed
    file_name = strcat('channel_data/ch',int2str(k),'.bkp');                %Locating the binary file recorded from testbed for CSI matrix (The file 'ch0.bkp' contains the 36×48 CSI matrix stored in array format, along with the paired UE list and the number of paired UEs.)
    r_ch_id=fopen(file_name,'r');                                           %Opening the binary file corresponding to above location
    r_ch_file=fread(r_ch_id,inf,'uint32');                                  %Reading the binary file of CSI matrix in unsigned int format
    fclose(r_ch_id);                                                        %Closing the binary file
    r_ch_bin_arr=de2bi(r_ch_file,32,'left-msb');                            %Converting the decimal datas to binary data
    r_ch_bin_big_endian_arr=little_to_big_endian_arr(r_ch_bin_arr);         %Performing little to big endian conversion
    r_ch_arr=bi_to_sp_float_arr(r_ch_bin_big_endian_arr);                   %Converting big endian binary to single precision float array, where the complex array is stored as interleaved real and imaginary samples.
    r_ch_c_arr=(r_ch_arr(1:2:2*nr*nt)+1i*r_ch_arr(2:2:2*nr*nt));            %Converting the real array to complex array
    H_ch=reshape(r_ch_c_arr,nt,nr).';                                       %Rearranging the elements of the complex array in matrix form to get the actual CSI matrix for kth subframe
    H=H_ch;                                                                 %CSI matrix of kth subframe
    H(7:8,:)=zeros(2,nt);                                                   %UEs 7 and 8 did not receive any signal and are considered nulls. Since the recorded data contains the minimum data-format value at these positions, these values are set to zero in post-processing.
    H(11:12,:)=zeros(2,nt);                                                 %UEs 11 and 12 did not receive any signal
    H(17:18,:)=zeros(2,nt);                                                 %UEs 17 to 18 did not receive any signal
    H(25:36,:)=zeros(12,nt);
    H_SF(:,:,k+1)=H;                                                        %CSI matrix for all subframes (48 subframes recorded in each experiment)
    nr_after_pairing_all_RT=[nr_after_pairing_all_RT,r_ch_arr(4*nr*nt+nr+1)]; %Number of paired UEs recorded from the testbed for threshold=0.3
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
%% Plotting Normalized Eigenvalues for different correlation thresholds
Normalized_Eigenvalue_CDF(n_e_val_Th);                                      %Plotting the CDF of normalized eigen values for different correlation thresholds
%% Plotting Sum spectral efficiency for different correlation thresholds
Sum_spectral_efficiency_CDF(nr_after_pairing_all_Th);                       %Plotting the CDF of sum spectral efficiency for different correlation thresholds
%% Downlink received constellation diagrams
REVM_th_03_qpsk_indoor=DL_received_constellation(3,'qpsk',23);              %Downlink received constellation diagrams for all paired UEs decoded with QPSK for correlation threshold=0.3
disp('REVM for all UEs (Threshold = 0.3, QPSK, User pairing)');                           
disp(["UE", "REVM"]);
disp(REVM_th_03_qpsk_indoor);                                               %Printing the REVM values for all paired UEs
%% Function definition for Plotting Normalized Eigenvalues for different correlation thresholds
function []=Normalized_Eigenvalue_CDF(n_e_val_Th) 

n_e_val_48x24_Th_0_3=n_e_val_Th{1};

set(0,'defaultAxesFontSize',14)
set(0,'defaultAxesFontName','Times')

figure
h1=cdfplot(n_e_val_48x24_Th_0_3);
set(h1,'Color','b','Marker','s','MarkerIndices',[1:80:500])

xlabel('Normalized eigenvalues');
ylabel('CDF');
title('');
legend("Indoor: Th=0.3",'Location','northwest');
% legend("Th=0.3","Th=0.2","Th=0.1",'Location','northwest');
end
%% Function definition for Plotting Sum spectral efficiency for different correlation thresholds
function []=Sum_spectral_efficiency_CDF(nr_after_pairing_all_Th)

C_RT_48x36_Th_0_3=2*nr_after_pairing_all_Th{1};

set(0,'defaultAxesFontSize',14)
set(0,'defaultAxesFontName','Times')

figure
h1=cdfplot(C_RT_48x36_Th_0_3);
set(h1,'Color','b','Marker','s','MarkerIndices',[1:1:10])

xlabel('Sum spectral efficiency (bits/sec/Hz)');
ylabel('CDF');
title('');
legend("Indoor: Th=0.3",'Location','northwest');
% legend("Th=0.3","Th=0.2","Th=0.1",'Location','northwest');
axis([0,30,0,1])
end
%% Function definition for Downlink received constellation diagrams
function [REVM_UE]=DL_received_constellation(threshold,modulation,subframe_num)
% This script reads downlink constellation data from binary files.
% Prints REVM for each paired UE, and plots received constellations for QPSK, 16-QAM, and 64-QAM under different threshold settings recorded from testbed.
% Inputs: threshold- correlation threshold considered in the experiment multipled by 10
%         modulation- modulation scheme considered in the experiment
%         subframe_num- Enter a subframe number between 1 and 36 to plot the corresponding data.
% Output: Plots constellation diagram for all paired UEs and returns their
%         corresponding REVM values.

file_name = strcat('dl_received_constellation/dat_SF',int2str(subframe_num),'.bkp');
dat_SF_id=fopen(file_name,'r');
dat_SF_file_arr=fread(dat_SF_id,inf,'uint32');
fclose(dat_SF_id);
dat_SF_bin_arr=de2bi(dat_SF_file_arr,32,'left-msb');
dat_SF_bin_big_endian_arr=little_to_big_endian_arr(dat_SF_bin_arr);
dat_SF=bi_to_sp_float_arr(dat_SF_bin_big_endian_arr);
%% Paired UE list
PUEL=[dat_SF(6951);dat_SF(6915:6915+dat_SF(6951)-1)];
if(subframe_num==23)
    paired_UE_list=[PUEL(2:6);PUEL(8:end)];
else
    paired_UE_list=PUEL(2:end);
end
np=length(paired_UE_list);
%% Modulation Index
for PUEL_ind=2:dat_SF(6951)+1
    MI(1,1)=dat_SF(6989)-1;
    MI(PUEL_ind,1)=dat_SF(6952+PUEL(PUEL_ind));
end
%% Receiver Error Vector Magnitude
for PUEL_ind=2:dat_SF(6951)+1
    EVM(1,1)=dat_SF(6952);
    EVM(PUEL_ind,1)=dat_SF(6990+PUEL(PUEL_ind));
end
% F_data=[PUEL,MI,EVM]
%% DL Received Constellation
sym=4;  %Total number of ofdm symbols +1
k1=7027;
k2=191;

figure;
for ue_index=1:np
ue=dat_SF(k1+k2*(paired_UE_list(ue_index)-1)+(paired_UE_list(ue_index)-1):k1+k2*(paired_UE_list(ue_index))+(paired_UE_list(ue_index)-1));
ue_1=(ue(1:(48*(sym-1))));
ue_cc=(ue_1(1:2:(48*(sym-1)))+1i*ue_1(2:2:(48*(sym-1))));
ue_c=ue_cc(1:end);
REVM(ue_index,1)=REVM_computation(ue_c,modulation);

set(0,'defaultAxesFontSize',14)
set(0,'defaultAxesFontName','Times')
xstring='InPhase';
ystring='QPhase';

if(strcmp(modulation, 'qpsk'))
    subplot(3,5,ue_index);
    scatter(real(ue_c)*sqrt(2),imag(ue_c)*sqrt(2));
    title(['UE',int2str(paired_UE_list(ue_index))]);
    grid on;
    xlabel(xstring);
    ylabel(ystring);
    axis([-2 2 -2 2])
    sgtitle('QPSK (Th=0.3) Indoor');
    box

elseif(strcmp(modulation, '16qam'))
    subplot(2,3,ue_index);
    scatter(real(ue_c)*sqrt(10),imag(ue_c)*sqrt(10));
    title(['UE',int2str(paired_UE_list(ue_index))]);
    grid on;
    xlabel(xstring);
    ylabel(ystring);
    axis([-4 4 -4 4])
    sgtitle('16QAM (Th=0.2) User pairing and Fair scheduling');
    box
elseif(strcmp(modulation, '64qam'))
    subplot(2,3,ue_index);
    scatter(real(ue_c)*sqrt(42),imag(ue_c)*sqrt(42));
    title(['UE',int2str(paired_UE_list(ue_index))]);
    grid on;
    xlabel(xstring);
    ylabel(ystring);
    axis([-10 10 -10 10])
    sgtitle('64QAM (Th=0.1) User pairing and Fair scheduling');
    box
end
end
REVM_UE=[paired_UE_list,REVM];
end
%% REVM Computation
function [EVM_ue]=REVM_computation(ue_c,modulation)
sym_num=6;      %number of ofdm symbols for data
sym_num_1=3;    %number of ofdm symbols for which EVM is calculated

%Actual QPSK transmitted data is denoted as ap_qpsk_9_sym, which will be used for REVM calculation
%48x36 9 ofdm data symbol
ap_qpsk_9_sym=[0.7071,0.7071,0.7071,-0.7071,-0.7071,0.7071,0.7071,0.7071,-0.7071,-0.7071,-0.7071,-0.7071,0.7071,-0.7071,-0.7071,-0.7071,-0.7071,-0.7071,-0.7071,0.7071,0.7071,-0.7071,0.7071,-0.7071,0.7071,-0.7071,-0.7071,0.7071,0.7071,-0.7071,-0.7071,-0.7071,-0.7071,-0.7071,0.7071,0.7071,-0.7071,-0.7071,0.7071,-0.7071,-0.7071,0.7071,0.7071,0.7071,0.7071,0.7071,0.7071,-0.7071,0.7071,-0.7071,0.7071,0.7071,-0.7071,0.7071,-0.7071,0.7071,0.7071,-0.7071,0.7071,-0.7071,0.7071,-0.7071,0.7071,-0.7071,-0.7071,-0.7071,-0.7071,0.7071,-0.7071,-0.7071,-0.7071,0.7071,-0.7071,-0.7071,0.7071,0.7071,0.7071,0.7071,0.7071,-0.7071,0.7071,-0.7071,0.7071,-0.7071,-0.7071,-0.7071,-0.7071,0.7071,-0.7071,0.7071,-0.7071,-0.7071,0.7071,-0.7071,0.7071,-0.7071,-0.7071,-0.7071,-0.7071,0.7071,-0.7071,-0.7071,0.7071,0.7071,-0.7071,0.7071,0.7071,0.7071,0.7071,-0.7071,-0.7071,0.7071,-0.7071,0.7071,-0.7071,-0.7071,-0.7071,0.7071,0.7071,-0.7071,-0.7071,0.7071,-0.7071,-0.7071,-0.7071,0.7071,0.7071,0.7071,-0.7071,-0.7071,-0.7071,0.7071,0.7071,-0.7071,-0.7071,0.7071,0.7071,0.7071,-0.7071,-0.7071,0.7071,-0.7071,-0.7071,-0.7071,0.7071,-0.7071,0.7071,-0.7071,-0.7071,-0.7071,0.7071,-0.7071,-0.7071,-0.7071,0.7071,-0.7071,0.7071,0.7071,-0.7071,0.7071,-0.7071,0.7071,0.7071,0.7071,-0.7071,0.7071,0.7071,0.7071,-0.7071,0.7071,0.7071,-0.7071,0.7071,-0.7071,0.7071,0.7071,0.7071,0.7071,0.7071,-0.7071,-0.7071,-0.7071,0.7071,-0.7071,0.7071,-0.7071,0.7071,0.7071,0.7071,0.7071,-0.7071,0.7071,-0.7071,-0.7071,0.7071,0.7071,-0.7071,-0.7071,-0.7071,0.7071,0.7071,0.7071,-0.7071,-0.7071,0.7071,-0.7071,-0.7071,-0.7071,-0.7071,-0.7071,0.7071,0.7071,-0.7071,-0.7071,-0.7071,0.7071,-0.7071,-0.7071,-0.7071,0.7071,-0.7071,0.7071,-0.7071,0.7071,0.7071,0.7071,-0.7071,-0.7071,0.7071,0.7071,-0.7071,0.7071,0.7071,0.7071,0.7071,-0.7071,0.7071,-0.7071,0.7071,0.7071,-0.7071,0.7071,0.7071,0.7071,-0.7071,-0.7071,0.7071,0.7071,-0.7071,-0.7071,0.7071,0.7071,0.7071,0.7071,0.7071,-0.7071,0.7071,-0.7071,-0.7071,-0.7071,-0.7071,-0.7071,-0.7071,-0.7071,0.7071,0.7071,0.7071,0.7071,0.7071,0.7071,-0.7071,0.7071,0.7071,-0.7071,-0.7071,-0.7071,-0.7071,-0.7071,-0.7071,0.7071,-0.7071,0.7071,-0.7071,-0.7071,0.7071,-0.7071,-0.7071,0.7071,0.7071,0.7071,-0.7071,-0.7071,-0.7071,0.7071,-0.7071,0.7071,0.7071,-0.7071,0.7071,0.7071,-0.7071,-0.7071,0.7071,0.7071,0.7071,0.7071,0.7071,-0.7071,0.7071,0.7071,-0.7071,0.7071,-0.7071,0.7071,0.7071,-0.7071,0.7071,-0.7071,-0.7071,0.7071,-0.7071,0.7071,0.7071,-0.7071,-0.7071,-0.7071,-0.7071,0.7071,-0.7071,-0.7071,0.7071,0.7071,-0.7071,0.7071,0.7071,0.7071,0.7071,0.7071,0.7071,-0.7071,-0.7071,-0.7071,0.7071,-0.7071,0.7071,-0.7071,0.7071,0.7071,-0.7071,-0.7071,0.7071,0.7071,0.7071,0.7071,0.7071,0.7071,-0.7071,0.7071,0.7071,0.7071,-0.7071,0.7071,0.7071,-0.7071,-0.7071,0.7071,0.7071,0.7071,0.7071,0.7071,0.7071,-0.7071,-0.7071,-0.7071,-0.7071,-0.7071,0.7071,0.7071,0.7071,-0.7071,0.7071,0.7071,0.7071,-0.7071,-0.7071,0.7071,-0.7071,0.7071,0.7071,-0.7071,-0.7071,0.7071,0.7071,0.7071,0.7071,0.7071,0.7071,-0.7071,-0.7071,-0.7071,-0.7071,0.7071,0.7071,0.7071,0.7071,0.7071,-0.7071,-0.7071,-0.7071,-0.7071,0.7071,0.7071,0.7071,0.7071,0.7071,-0.7071,0.7071,0.7071,-0.7071,-0.7071,-0.7071,-0.7071,-0.7071,-0.7071,0.7071,0.7071,0.7071,-0.7071,-0.7071,-0.7071,0.7071,-0.7071];
ap_qpsk=ap_qpsk_9_sym(1:24*sym_num*2);
ap_qpsk_c=ap_qpsk(1:2:end)+1i*ap_qpsk(2:2:end);
ap_qpsk_rearranged = reshape(ap_qpsk_c,sym_num,24);
ap_qpsk_array = reshape(ap_qpsk_rearranged.',1,sym_num*24).';
ap_qpsk_Tx = ap_qpsk_array(1:sym_num_1*24);

%Actual 16QAM transmitted data is denoted as ap_16qam_9_sym, which will be used for REVM calculation
%48x36 9 ofdm data symbol
ap_16qam_9_sym=[0.3162,0.9486,0.3162,-0.9486,0.3162,0.3162,0.3162,-0.3162,0.9486,0.9486,0.9486,-0.9486,0.9486,0.3162,0.9486,-0.3162,-0.9486,0.9486,-0.9486,-0.9486,-0.9486,0.3162,-0.9486,-0.3162,-0.3162,0.9486,-0.3162,-0.9486,-0.3162,0.3162,-0.3162,-0.3162,0.3162,0.9486,0.3162,-0.9486,0.3162,0.3162,0.3162,-0.3162,0.9486,0.9486,0.9486,-0.9486,0.9486,0.3162,0.9486,-0.3162,-0.9486,0.9486,-0.9486,-0.9486,-0.9486,0.3162,-0.9486,-0.3162,-0.3162,0.9486,-0.3162,-0.9486,-0.3162,0.3162,-0.3162,-0.3162,0.3162,0.9486,0.3162,-0.9486,0.3162,0.3162,0.3162,-0.3162,0.9486,0.9486,0.9486,-0.9486,0.9486,0.3162,0.9486,-0.3162,-0.9486,0.9486,-0.9486,-0.9486,-0.9486,0.3162,-0.9486,-0.3162,-0.3162,0.9486,-0.3162,-0.9486,-0.3162,0.3162,-0.3162,-0.3162,0.3162,0.9486,0.3162,-0.9486,0.3162,0.3162,0.3162,-0.3162,0.9486,0.9486,0.9486,-0.9486,0.9486,0.3162,0.9486,-0.3162,-0.9486,0.9486,-0.9486,-0.9486,-0.9486,0.3162,-0.9486,-0.3162,-0.3162,0.9486,-0.3162,-0.9486,-0.3162,0.3162,-0.3162,-0.3162,0.3162,0.9486,0.3162,-0.9486,0.3162,0.3162,0.3162,-0.3162,0.9486,0.9486,0.9486,-0.9486,0.9486,0.3162,0.9486,-0.3162,-0.9486,0.9486,-0.9486,-0.9486,-0.9486,0.3162,-0.9486,-0.3162,-0.3162,0.9486,-0.3162,-0.9486,-0.3162,0.3162,-0.3162,-0.3162,0.3162,0.9486,0.3162,-0.9486,0.3162,0.3162,0.3162,-0.3162,0.9486,0.9486,0.9486,-0.9486,0.9486,0.3162,0.9486,-0.3162,-0.9486,0.9486,-0.9486,-0.9486,-0.9486,0.3162,-0.9486,-0.3162,-0.3162,0.9486,-0.3162,-0.9486,-0.3162,0.3162,-0.3162,-0.3162,0.3162,0.9486,0.3162,-0.9486,0.3162,0.3162,0.3162,-0.3162,0.9486,0.9486,0.9486,-0.9486,0.9486,0.3162,0.9486,-0.3162,-0.9486,0.9486,-0.9486,-0.9486,-0.9486,0.3162,-0.9486,-0.3162,-0.3162,0.9486,-0.3162,-0.9486,-0.3162,0.3162,-0.3162,-0.3162,0.3162,0.9486,0.3162,-0.9486,0.3162,0.3162,0.3162,-0.3162,0.9486,0.9486,0.9486,-0.9486,0.9486,0.3162,0.9486,-0.3162,-0.9486,0.9486,-0.9486,-0.9486,-0.9486,0.3162,-0.9486,-0.3162,-0.3162,0.9486,-0.3162,-0.9486,-0.3162,0.3162,-0.3162,-0.3162,0.3162,0.9486,0.3162,-0.9486,0.3162,0.3162,0.3162,-0.3162,0.9486,0.9486,0.9486,-0.9486,0.9486,0.3162,0.9486,-0.3162,-0.9486,0.9486,-0.9486,-0.9486,-0.9486,0.3162,-0.9486,-0.3162,-0.3162,0.9486,-0.3162,-0.9486,-0.3162,0.3162,-0.3162,-0.3162,0.3162,0.9486,0.3162,-0.9486,0.3162,0.3162,0.3162,-0.3162,0.9486,0.9486,0.9486,-0.9486,0.9486,0.3162,0.9486,-0.3162,-0.9486,0.9486,-0.9486,-0.9486,-0.9486,0.3162,-0.9486,-0.3162,-0.3162,0.9486,-0.3162,-0.9486,-0.3162,0.3162,-0.3162,-0.3162,0.3162,0.9486,0.3162,-0.9486,0.3162,0.3162,0.3162,-0.3162,0.9486,0.9486,0.9486,-0.9486,0.9486,0.3162,0.9486,-0.3162,-0.9486,0.9486,-0.9486,-0.9486,-0.9486,0.3162,-0.9486,-0.3162,-0.3162,0.9486,-0.3162,-0.9486,-0.3162,0.3162,-0.3162,-0.3162,0.3162,0.9486,0.3162,-0.9486,0.3162,0.3162,0.3162,-0.3162,0.9486,0.9486,0.9486,-0.9486,0.9486,0.3162,0.9486,-0.3162,-0.9486,0.9486,-0.9486,-0.9486,-0.9486,0.3162,-0.9486,-0.3162,-0.3162,0.9486,-0.3162,-0.9486,-0.3162,0.3162,-0.3162,-0.3162,0.3162,0.9486,0.3162,-0.9486,0.3162,0.3162,0.3162,-0.3162,0.9486,0.9486,0.9486,-0.9486,0.9486,0.3162,0.9486,-0.3162,-0.9486,0.9486,-0.9486,-0.9486,-0.9486,0.3162,-0.9486,-0.3162,-0.3162,0.9486,-0.3162,-0.9486,-0.3162,0.3162,-0.3162,-0.3162,0.3162,0.9486,0.3162,-0.9486,0.3162,0.3162,0.3162,-0.3162,0.9486,0.9486,0.9486,-0.9486,0.9486,0.3162,0.9486,-0.3162];
ap_16qam=ap_16qam_9_sym(1:24*sym_num*2);
ap_16qam_c=ap_16qam(1:2:end)+1i*ap_16qam(2:2:end);
ap_16qam_rearranged = reshape(ap_16qam_c,sym_num,24);
ap_16qam_array = reshape(ap_16qam_rearranged.',1,sym_num*24).';
ap_16qam_Tx = ap_16qam_array(1:sym_num_1*24);

%Actual 64QAM transmitted data is denoted as ap_64qam_9_sym, which will be used for REVM calculation
%48x36 9 ofdm data symbol
ap_64qam_9_sym=[-1.0801,1.0801,-1.0801,0.7715,-1.0801,0.1543,-1.0801,0.4629,-1.0801,-1.0801,-1.0801,-0.7715,-1.0801,-0.1543,-1.0801,-0.4629,-0.7715,1.0801,-0.7715,0.7715,-0.7715,0.1543,-0.7715,0.4629,-0.7715,-1.0801,-0.7715,-0.7715,-0.7715,-0.1543,-0.7715,-0.4629,-0.1543,1.0801,-0.1543,0.7715,-0.1543,0.1543,-0.1543,0.4629,-0.1543,-1.0801,-0.1543,-0.7715,-0.1543,-0.1543,-0.1543,-0.4629,-0.4629,1.0801,-0.4629,0.7715,-0.4629,0.1543,-0.4629,0.4629,-0.4629,-1.0801,-0.4629,-0.7715,-0.4629,-0.1543,-0.4629,-0.4629,1.0801,1.0801,1.0801,0.7715,1.0801,0.1543,1.0801,0.4629,1.0801,-1.0801,1.0801,-0.7715,1.0801,-0.1543,1.0801,-0.4629,0.7715,1.0801,0.7715,0.7715,0.7715,0.1543,0.7715,0.4629,0.7715,-1.0801,0.7715,-0.7715,0.7715,-0.1543,0.7715,-0.4629,0.1543,1.0801,0.1543,0.7715,0.1543,0.1543,0.1543,0.4629,0.1543,-1.0801,0.1543,-0.7715,0.1543,-0.1543,0.1543,-0.4629,0.4629,1.0801,0.4629,0.7715,0.4629,0.1543,0.4629,0.4629,0.4629,-1.0801,0.4629,-0.7715,0.4629,-0.1543,0.4629,-0.4629,-1.0801,1.0801,-1.0801,0.7715,-1.0801,0.1543,-1.0801,0.4629,-1.0801,-1.0801,-1.0801,-0.7715,-1.0801,-0.1543,-1.0801,-0.4629,-0.7715,1.0801,-0.7715,0.7715,-0.7715,0.1543,-0.7715,0.4629,-0.7715,-1.0801,-0.7715,-0.7715,-0.7715,-0.1543,-0.7715,-0.4629,-0.1543,1.0801,-0.1543,0.7715,-0.1543,0.1543,-0.1543,0.4629,-0.1543,-1.0801,-0.1543,-0.7715,-0.1543,-0.1543,-0.1543,-0.4629,-0.4629,1.0801,-0.4629,0.7715,-0.4629,0.1543,-0.4629,0.4629,-0.4629,-1.0801,-0.4629,-0.7715,-0.4629,-0.1543,-0.4629,-0.4629,1.0801,1.0801,1.0801,0.7715,1.0801,0.1543,1.0801,0.4629,1.0801,-1.0801,1.0801,-0.7715,1.0801,-0.1543,1.0801,-0.4629,0.7715,1.0801,0.7715,0.7715,0.7715,0.1543,0.7715,0.4629,0.7715,-1.0801,0.7715,-0.7715,0.7715,-0.1543,0.7715,-0.4629,0.1543,1.0801,0.1543,0.7715,0.1543,0.1543,0.1543,0.4629,0.1543,-1.0801,0.1543,-0.7715,0.1543,-0.1543,0.1543,-0.4629,0.4629,1.0801,0.4629,0.7715,0.4629,0.1543,0.4629,0.4629,0.4629,-1.0801,0.4629,-0.7715,0.4629,-0.1543,0.4629,-0.4629,-1.0801,1.0801,-1.0801,0.7715,-1.0801,0.1543,-1.0801,0.4629,-1.0801,-1.0801,-1.0801,-0.7715,-1.0801,-0.1543,-1.0801,-0.4629,-0.7715,1.0801,-0.7715,0.7715,-0.7715,0.1543,-0.7715,0.4629,-0.7715,-1.0801,-0.7715,-0.7715,-0.7715,-0.1543,-0.7715,-0.4629,-0.1543,1.0801,-0.1543,0.7715,-0.1543,0.1543,-0.1543,0.4629,-0.1543,-1.0801,-0.1543,-0.7715,-0.1543,-0.1543,-0.1543,-0.4629,-0.4629,1.0801,-0.4629,0.7715,-0.4629,0.1543,-0.4629,0.4629,-0.4629,-1.0801,-0.4629,-0.7715,-0.4629,-0.1543,-0.4629,-0.4629,1.0801,1.0801,1.0801,0.7715,1.0801,0.1543,1.0801,0.4629,1.0801,-1.0801,1.0801,-0.7715,1.0801,-0.1543,1.0801,-0.4629,0.7715,1.0801,0.7715,0.7715,0.7715,0.1543,0.7715,0.4629,0.7715,-1.0801,0.7715,-0.7715,0.7715,-0.1543,0.7715,-0.4629,0.1543,1.0801,0.1543,0.7715,0.1543,0.1543,0.1543,0.4629,0.1543,-1.0801,0.1543,-0.7715,0.1543,-0.1543,0.1543,-0.4629,0.4629,1.0801,0.4629,0.7715,0.4629,0.1543,0.4629,0.4629,0.4629,-1.0801,0.4629,-0.7715,0.4629,-0.1543,0.4629,-0.4629,-1.0801,1.0801,-1.0801,0.7715,-1.0801,0.1543,-1.0801,0.4629,-1.0801,-1.0801,-1.0801,-0.7715,-1.0801,-0.1543,-1.0801,-0.4629,-0.7715,1.0801,-0.7715,0.7715,-0.7715,0.1543,-0.7715,0.4629,-0.7715,-1.0801,-0.7715,-0.7715,-0.7715,-0.1543,-0.7715,-0.4629,-0.1543,1.0801,-0.1543,0.7715,-0.1543,0.1543,-0.1543,0.4629,-0.1543,-1.0801,-0.1543,-0.7715,-0.1543,-0.1543,-0.1543,-0.4629];
ap_64qam=ap_64qam_9_sym(1:24*sym_num*2);
ap_64qam_c=ap_64qam(1:2:end)+1i*ap_64qam(2:2:end);
ap_64qam_rearranged = reshape(ap_64qam_c,sym_num,24);
ap_64qam_array = reshape(ap_64qam_rearranged.',1,sym_num*24).';
ap_64qam_Tx = ap_64qam_array(1:sym_num_1*24);

if(strcmp(modulation, 'qpsk'))
    EVM_ue = sqrt(sum(abs(ue_c-ap_qpsk_Tx).^2)/length(ap_qpsk_Tx))*100;
elseif(strcmp(modulation, '16qam'))
    EVM_ue = sqrt(sum(abs(ue_c-ap_16qam_Tx).^2)/length(ap_16qam_Tx))*100;
elseif(strcmp(modulation, '64qam'))
    EVM_ue = sqrt(sum(abs(ue_c-ap_64qam_Tx).^2)/length(ap_64qam_Tx))*100;
end
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