clc
clear
close all

nt=48;nr=36;
Threshold=0.3;
idx=1;
for expt=[1 2 3]
N=47;
n_e_val=[];
nr_after_pairing_all=[];
nr_after_pairing_all_RT=[];
paired_UE_list_all=[];
% H_SF=zeros(nr,nt,N);
for k=0:N
%% Reading binary files for recorded CSI from testbed
if(expt==1)
    file_name = strcat('channel_data/Expt_',int2str(expt),'/ch',int2str(k),'.bkp');     %Locating the binary file recorded from testbed for CSI matrix (The file 'ch0.bkp' contains the 36×48 CSI matrix stored in array format, along with the paired UE list and the number of paired UEs.)
    r_ch_id=fopen(file_name,'r');                                           %Opening the binary file corresponding to above location
    r_ch_file=fread(r_ch_id,inf,'uint32');                                  %Reading the binary file of CSI matrix in unsigned int format
    fclose(r_ch_id);                                                        %Closing the binary file
    r_ch_bin_arr=de2bi(r_ch_file,32,'left-msb');                            %Converting the decimal datas to binary data
    r_ch_bin_big_endian_arr=little_to_big_endian_arr(r_ch_bin_arr);         %Performing little to big endian conversion
    r_ch_arr=bi_to_sp_float_arr(r_ch_bin_big_endian_arr);                   %Converting big endian binary to single precision float array, where the complex array is stored as interleaved real and imaginary samples.
    r_ch_c_arr=(r_ch_arr(1:2:2*nr*nt)+1i*r_ch_arr(2:2:2*nr*nt));            %Converting the real array to complex array
    H_ch=reshape(r_ch_c_arr,nt,nr).';                                       %Rearranging the elements of the complex array in matrix form to get the actual CSI matrix for kth subframe
    H=H_ch;                                                                 %CSI matrix of kth subframe
    H(5:6,:)=zeros(2,nt);                                                   %UEs 5 and 6 did not receive any signal and are considered nulls. Since the recorded data contains the minimum data-format value at these positions, these values are set to zero in post-processing.
    H(23:24,:)=zeros(2,nt);                                                 %UEs 23 and 24 did not receive any signal
    H(25:36,:)=zeros(12,nt);                                                %UEs 25 to 36 did not receive any signal
    nr_after_pairing_all_RT=[nr_after_pairing_all_RT,r_ch_arr(4*nr*nt+nr+1)]; %Number of paired UEs recorded from the testbed for threshold=0.3
elseif(expt==2)
    file_name = strcat('channel_data/Expt_',int2str(expt),'/ch',int2str(k),'.bkp');     
    r_ch_id=fopen(file_name,'r');                                           
    r_ch_file=fread(r_ch_id,inf,'uint32');                                  
    fclose(r_ch_id);                                                        
    r_ch_bin_arr=de2bi(r_ch_file,32,'left-msb');                            
    r_ch_bin_big_endian_arr=little_to_big_endian_arr(r_ch_bin_arr);         
    r_ch_arr=bi_to_sp_float_arr(r_ch_bin_big_endian_arr);                   
    r_ch_c_arr=(r_ch_arr(1:2:2*nr*nt)+1i*r_ch_arr(2:2:2*nr*nt));            
    H_ch=reshape(r_ch_c_arr,nt,nr).';                                       %Rearranging the elements of the complex array in matrix form to get the actual CSI matrix for kth subframe
    H=H_ch;
    H(7:8,:)=zeros(2,nt);                                                   %UEs 7 and 8 did not receive any signal and are considered nulls. Since the recorded data contains the minimum data-format value at these positions, these values are set to zero in post-processing.
    H(23:24,:)=zeros(2,nt);                                                 %UEs 23 and 24 did not receive any signal
    H(25:36,:)=zeros(12,nt);                                                %UEs 25 to 36 did not receive any signal
    nr_after_pairing_all_RT=[nr_after_pairing_all_RT,r_ch_arr(4*nr*nt+nr+1)]; %Number of paired UEs recorded from the testbed for threshold=0.2
elseif(expt==3)
    file_name = strcat('channel_data/Expt_',int2str(expt),'/ch',int2str(k),'.bkp');     
    r_ch_id=fopen(file_name,'r');                                           
    r_ch_file=fread(r_ch_id,inf,'uint32');                                  
    fclose(r_ch_id);                                                        
    r_ch_bin_arr=de2bi(r_ch_file,32,'left-msb');                            
    r_ch_bin_big_endian_arr=little_to_big_endian_arr(r_ch_bin_arr);         
    r_ch_arr=bi_to_sp_float_arr(r_ch_bin_big_endian_arr);                   
    r_ch_c_arr=(r_ch_arr(1:2:2*nr*nt)+1i*r_ch_arr(2:2:2*nr*nt));            
    H_ch=reshape(r_ch_c_arr,nt,nr).';                                       %Rearranging the elements of the complex array in matrix form to get the actual CSI matrix for kth subframe
    H=H_ch;
    H(11:12,:)=zeros(2,nt);                                                 %UEs 11 and 12 did not receive any signal and are considered nulls. Since the recorded data contains the minimum data-format value at these positions, these values are set to zero in post-processing.
    H(23:24,:)=zeros(2,nt);                                                 %UEs 23 and 24 did not receive any signal
    H(25:36,:)=zeros(12,nt);                                                %UEs 25 to 36 did not receive any signal
    nr_after_pairing_all_RT=[nr_after_pairing_all_RT,r_ch_arr(4*nr*nt+nr+1)]; %Number of paired UEs recorded from the testbed for threshold=0.1
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
%% Plotting Normalized Eigenvalues for different correlation thresholds
% Normalized_Eigenvalue_CDF(n_e_val_Th);                                    %Plotting the CDF of normalized eigen values for different correlation thresholds
%% Plotting Sum spectral efficiency for different correlation thresholds
Sum_spectral_efficiency_CDF(nr_after_pairing_all_Th);                       %Plotting the CDF of sum spectral efficiency for different correlation thresholds
%% Downlink received constellation diagrams
REVM_th_03_qpsk_E1=DL_received_constellation(3,'qpsk',21,1);                %Downlink received constellation diagrams for all paired UEs decoded with QPSK for correlation threshold=0.3
disp('REVM for all UEs in Expt 1');                           
disp(["UE", "REVM"]);
disp(REVM_th_03_qpsk_E1);                                                   %Printing the REVM values for all paired UEs
REVM_th_03_qpsk_E2=DL_received_constellation(3,'qpsk',24,2);                %Downlink received constellation diagrams for all paired UEs decoded with QPSK for correlation threshold=0.3
disp('REVM for all UEs in Expt 2');                           
disp(["UE", "REVM"]);
disp(REVM_th_03_qpsk_E2);                                                   %Printing the REVM values for all paired UEs
REVM_th_02_16qam_E3=DL_received_constellation(3,'qpsk',24,3);               %Downlink received constellation diagrams for all paired UEs decoded with 16QAM for correlation threshold=0.2
disp('REVM for all UEs in Expt 3');
disp(["UE", "REVM"]);
disp(REVM_th_02_16qam_E3);                                                  %Printing the REVM values for all paired UEs
%% Number of paired UEs for different experiments
Np_bar();
%% Plotting REVM for different experiments
REVM_CDF();
%% Normalized pairing opportunity for all UEs in different experiments
Pairing_opportunity();
%% Function definition for Plotting Normalized Eigenvalues for different correlation thresholds
function []=Normalized_Eigenvalue_CDF(n_e_val_Th) 

n_e_val_48x24_Th_E1=n_e_val_Th{1};
n_e_val_48x24_Th_E2=n_e_val_Th{2};
n_e_val_48x24_Th_E3=n_e_val_Th{3};

set(0,'defaultAxesFontSize',14)
set(0,'defaultAxesFontName','Times')

figure
h1=cdfplot(n_e_val_48x24_Th_E1);
set(h1,'Color','b','Marker','s','MarkerIndices',[1:80:800])
hold on
h2=cdfplot(n_e_val_48x24_Th_E2);
set(h2,'Color','k','Marker','d','MarkerIndices',[1:50:500])
h3=cdfplot(n_e_val_48x24_Th_E3);
set(h3,'Color','r','Marker','x','MarkerIndices',[1:20:200])

xlabel('Normalized eigenvalues');
ylabel('CDF');
title('');
legend("Th=0.3","Th=0.2","Th=0.1",'Location','northwest');
end
%% Function definition for Plotting Sum spectral efficiency for different correlation thresholds
function []=Sum_spectral_efficiency_CDF(nr_after_pairing_all_Th)

C_RT_48x24_Th_E1=2*nr_after_pairing_all_Th{1};
C_RT_48x24_Th_E2=2*(nr_after_pairing_all_Th{2}-2);
C_RT_48x24_Th_E3=2*nr_after_pairing_all_Th{3};

set(0,'defaultAxesFontSize',14)
set(0,'defaultAxesFontName','Times')

figure
h1=cdfplot(C_RT_48x24_Th_E1);
set(h1,'Color','b','Marker','s','MarkerIndices',[1:1:10])
hold on
h2=cdfplot(C_RT_48x24_Th_E2);
set(h2,'Color','k','Marker','d','MarkerIndices',[1:1:10])
h3=cdfplot(C_RT_48x24_Th_E3);
set(h3,'Color','r','Marker','x','MarkerIndices',[1:1:10])

xlabel('Sum spectral efficiency (bits/sec/Hz)');
ylabel('CDF');
title('');
legend("Expt 1","Expt 2","Expt 3",'Location','northwest');
axis([0,30,0,1])
end
%% Function definition for Downlink received constellation diagrams
function [REVM_UE]=DL_received_constellation(threshold,modulation,subframe_num,expt_num)
% This script reads downlink constellation data from binary files.
% Prints REVM for each paired UE, and plots received constellations for QPSK, 16-QAM, and 64-QAM under different threshold settings recorded from testbed.
% Inputs: threshold- correlation threshold considered in the experiment multipled by 10
%         modulation- modulation scheme considered in the experiment
%         subframe_num- Enter a subframe number between 1 and 36 to plot the corresponding data.
%         expt_num- Enter the experiment number between 1 and 4
% Output: Plots constellation diagram for all paired UEs and returns their
%         corresponding REVM values.

file_name = strcat('dl_received_constellation/Expt_',int2str(expt_num),'/dat_SF',int2str(subframe_num),'.bkp');
dat_SF_id=fopen(file_name,'r');
dat_SF_file_arr=fread(dat_SF_id,inf,'uint32');
fclose(dat_SF_id);
dat_SF_bin_arr=de2bi(dat_SF_file_arr,32,'left-msb');
dat_SF_bin_big_endian_arr=little_to_big_endian_arr(dat_SF_bin_arr);
dat_SF=bi_to_sp_float_arr(dat_SF_bin_big_endian_arr);
%% Paired UE list
PUEL=[dat_SF(6951);dat_SF(6915:6915+dat_SF(6951)-1)];
if(expt_num==2 && subframe_num==24)
    paired_UE_list=PUEL(2:end-2);
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
k1=7027;    %Starting index of dl received data
k2=191;     %length of data for a UE

figure;
for ue_index=1:np
ue=dat_SF(k1+k2*(paired_UE_list(ue_index)-1)+(paired_UE_list(ue_index)-1):k1+k2*(paired_UE_list(ue_index))+(paired_UE_list(ue_index)-1));
ue_1=(ue(1:(48*(sym-1))));
ue_cc=(ue_1(1:2:(48*(sym-1)))+1i*ue_1(2:2:(48*(sym-1))));
ue_c=ue_cc(1:end);
REVM(ue_index,1)=REVM_computation(ue_c,modulation);                         %Computing REVM for a UE

set(0,'defaultAxesFontSize',14)
set(0,'defaultAxesFontName','Times')
xstring='InPhase';
ystring='QPhase';

if(strcmp(modulation, 'qpsk'))
    if(expt_num==1)
        subplot(3,4,ue_index);
    elseif(expt_num==2)
        subplot(3,4,ue_index);
    elseif(expt_num==3)
        subplot(2,3,ue_index);
    end
    scatter(real(ue_c)*sqrt(2),imag(ue_c)*sqrt(2));
    title(['UE',int2str(paired_UE_list(ue_index))]);
    grid on;
    xlabel(xstring);
    ylabel(ystring);
    axis([-2 2 -2 2])
    if(expt_num==1)
        sgtitle('Expt 1');
    elseif(expt_num==2)
        sgtitle('Expt 2');
    elseif(expt_num==3)
        sgtitle('Expt 3');
    end
    box

elseif(strcmp(modulation, '16qam'))
    if(expt_num==1)
        subplot(3,4,ue_index);
    elseif(expt_num==2)
        subplot(2,4,ue_index);
    elseif(expt_num==3)
        subplot(2,3,ue_index);
    end
    scatter(real(ue_c)*sqrt(10),imag(ue_c)*sqrt(10));
    title(['UE',int2str(paired_UE_list(ue_index))]);
    grid on;
    xlabel(xstring);
    ylabel(ystring);
    axis([-4 4 -4 4])
    if(expt_num==1)
        sgtitle('Expt 1');
    elseif(expt_num==2)
        sgtitle('Expt 2');
    elseif(expt_num==3)
        sgtitle('Expt 3');
    end
    box
elseif(strcmp(modulation, '64qam'))
    if(expt_num==1)
        subplot(3,4,ue_index);
    elseif(expt_num==2)
        subplot(2,4,ue_index);
    elseif(expt_num==3)
        subplot(2,3,ue_index);
    end
    scatter(real(ue_c)*sqrt(42),imag(ue_c)*sqrt(42));
    title(['UE',int2str(paired_UE_list(ue_index))]);
    grid on;
    xlabel(xstring);
    ylabel(ystring);
    axis([-10 10 -10 10])
    if(expt_num==1)
        sgtitle('Expt 1');
    elseif(expt_num==2)
        sgtitle('Expt 2');
    elseif(expt_num==3)
        sgtitle('Expt 3');
    end
    box
end
end
REVM_UE=[paired_UE_list(1:np),REVM];
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
%% REVM CDF
function []=REVM_CDF()
load("dl_received_constellation/Expt_1/REVM_E1.mat")
load("dl_received_constellation/Expt_2/REVM_E2.mat")
load("dl_received_constellation/Expt_3/REVM_E3.mat")

set(0,'defaultAxesFontSize',14)
set(0,'defaultAxesFontName','Times')

figure
R_E1=cdfplot(REVM_E1);
set(R_E1,'Color','b','Marker','s','MarkerIndices',[1:20:360])
hold on
R_E2=cdfplot(REVM_E2);
set(R_E2,'Color','k','Marker','d','MarkerIndices',[1:20:280])
R_E3=cdfplot(REVM_E3);
set(R_E3,'Color','r','Marker','x','MarkerIndices',[1:10:120])

xlabel('REVM for all paired UEs');
ylabel('CDF');
title('');
legend("Expt-1","Expt-2","Expt-3",'Location','southeast');
end
%% Bar Diagram for number of paired UEs
function []=Np_bar()
figure
set(0,'defaultAxesFontSize',16)
set(0,'defaultAxesFontName','Times')

%The number of paired UEs observed across all subframes varies between 9 to 11 for experiment 1 denoted as N_p_E1 
%and pairing_freq_E1 denotes their occurance frequency respectively
N_p_E1=[12,11,10,9];
pairing_freq_E1=[1,5,10,4]/20;
hold on

N_p_E2=[8,7,6,5];
pairing_freq_E2=[5,7,7,1]/20;

N_p_E3=[6,5,4,0];
pairing_freq_E3=[3,7,3,0]/13;

% N_p_E4=[3,2,1];
% pairing_freq_E4=[7,6,1]/14;

N_p=[N_p_E1; N_p_E2; N_p_E3];
pairing_freq=[pairing_freq_E1; pairing_freq_E2; pairing_freq_E3];
bar(N_p.', pairing_freq.')
box
axis([3,12,0,1])
grid on
xlabel("Number of users paired in each time instant");
ylabel("Probability");
legend("Expt-1","Expt-2","Expt-3")
end
%% Bar Diagram for normalized pairing opportunity for all UEs in different experiments
function []=Pairing_opportunity()
%The list of all 24 UEs and their respective pairing opportunity recorded during the experiment are captured below for all four experiments
UE_list_E1=[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22];
Pairing_opportunity_E1=[1,18,17,3,0,0,15,1,20,19,2,17,19,19,5,5,11,13,17,3,1,5]/22;
UE_list_E2=[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22];
Pairing_opportunity_E2=[9,8,10,10,10,9,0,0,9,5,10,10,10,10,10,11,10,6,8,4,10,10]/22;
UE_list_E3=[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22];
Pairing_opportunity_E3=[6,7,0,0,7,5,4,6,2,2,0,0,0,0,0,3,0,0,6,7,8,1]/13;

figure;
set(0,'defaultAxesFontSize',14)
set(0,'defaultAxesFontName','Times')

subplot(2,2,1)
bar(UE_list_E1, Pairing_opportunity_E1)
axis([0,25,0,1])
xlabel("UE numbers","FontSize",20);
ylabel({"Normalized"; "Pairing opportunity"},"FontSize",20);
grid on
title("Expt-1")

subplot(2,2,2)
bar(UE_list_E2, Pairing_opportunity_E2)
axis([0,25,0,1])
xlabel("UE numbers","FontSize",20);
ylabel({"Normalized"; "Pairing opportunity"},"FontSize",20);
grid on
title("Expt-2")

subplot(2,2,3)
bar(UE_list_E3, Pairing_opportunity_E3)
axis([0,25,0,1])
xlabel("UE numbers","FontSize",20);
ylabel({"Normalized"; "Pairing opportunity"},"FontSize",20);
grid on
title("Expt-3")

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