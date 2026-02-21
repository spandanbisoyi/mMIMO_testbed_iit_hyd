clc
clear
close all

nt=24;nr=24;
idx=1;
Threshold=0.3;
for scenario=1:2
N=47;
n_e_val=[];
nr_after_pairing_all=[];
nr_after_pairing_all_RT=[];
paired_UE_list_all=[];
% H_SF=zeros(nr,nt,N);
for k=0:N
%% Reading binary files for recorded CSI from testbed
if(scenario==1)                                                             %scenario 1 refers to the indoor scenario
    file_name = strcat('channel_data/indoor/ch',int2str(k),'.bkp');         %Locating the binary file recorded from testbed for CSI matrix (The file 'ch0.bkp' contains the 36×48 CSI matrix stored in array format, along with the paired UE list and the number of paired UEs.)
    r_ch_id=fopen(file_name,'r');                                           %Opening the binary file corresponding to above location
    r_ch_file=fread(r_ch_id,inf,'uint32');                                  %Reading the binary file of CSI matrix in unsigned int format
    fclose(r_ch_id);                                                        %Closing the binary file
    r_ch_bin_arr=de2bi(r_ch_file,32,'left-msb');                            %Converting the decimal datas to binary data
    r_ch_bin_big_endian_arr=little_to_big_endian_arr(r_ch_bin_arr);         %Performing little to big endian conversion
    r_ch_arr=bi_to_sp_float_arr(r_ch_bin_big_endian_arr);                   %Converting big endian binary to single precision float array, where the complex array is stored as interleaved real and imaginary samples.
    r_ch_c_arr=(r_ch_arr(1:2:2*nr*nt)+1i*r_ch_arr(2:2:2*nr*nt));            %Converting the real array to complex array
    H_ch=reshape(r_ch_c_arr,nt,nr).';                                       %Rearranging the elements of the complex array in matrix form to get the actual CSI matrix for kth subframe
    H=H_ch;                                                                 %CSI matrix of kth subframe
elseif(scenario==2)                                                         %scenario 2 refers to the outdoor scenario
    file_name = strcat('channel_data/outdoor/ch',int2str(k),'.bkp');
    r_ch_id=fopen(file_name,'r');
    r_ch_file=fread(r_ch_id,inf,'uint32');
    fclose(r_ch_id);
    r_ch_bin_arr=de2bi(r_ch_file,32,'left-msb');
    r_ch_bin_big_endian_arr=little_to_big_endian_arr(r_ch_bin_arr);
    r_ch_arr=bi_to_sp_float_arr(r_ch_bin_big_endian_arr);
    r_ch_c_arr=(r_ch_arr(1:2:2*nr*nt)+1i*r_ch_arr(2:2:2*nr*nt));
    H_ch=reshape(r_ch_c_arr,nt,nr).';
    H=H_ch;                         
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
paired_UE_list_indoor=[1,2,5,6,9,10,11,13,15,17,18,20,22,23,24];
paired_UE_list_outdoor=[1,2,3,4,9,11,12,20,22];
%% Plotting Normalized Eigenvalues for different scenarios with correlation threshold 0.3
Normalized_Eigenvalue_CDF(n_e_val_Th);                                      %Plotting the CDF of normalized eigen values for different correlation thresholds
%% Sum spectral efficiency for different scenarios with correlation threshold 0.3
cidx=1;
for scenario=1:2
N_p_all=[];
for k=0:23
    if(scenario==1)
        file_name = strcat('frame_data/indoor/frame_data',int2str(k),'.bkp');
    elseif(scenario==2)
        file_name = strcat('frame_data/outdoor/frame_data',int2str(k),'.bkp');
    end
    frame_data_id=fopen(file_name,'r');
    frame_data_file=fread(frame_data_id,inf,'uint32');
    fclose(frame_data_id);
    frame_data_bin_arr=de2bi(frame_data_file,32,'left-msb');
    frame_data_bin_big_endian_arr=little_to_big_endian_arr(frame_data_bin_arr);
    frame_data_arr=bi_to_sp_float_arr(frame_data_bin_big_endian_arr);
    for FNum=0
        %Reading the paired UE list recorded in frame_data_arr for multiple subframes. FNum is the subframe index. frame_data_arr(37+nr+1+110*FNum) is number of paired UEs.
        PUEL=[frame_data_arr(37+nr+1+110*FNum);frame_data_arr(38+110*FNum:38+110*FNum+frame_data_arr(37+nr+1+110*FNum)-1)];
        %Reading the modulation index for all UEs across multiple subframes. frame_data_arr(1+110*FNum) is the subframe number
        for PUEL_ind=2:frame_data_arr(37+nr+1+110*FNum)+1
            MI(1,1)=frame_data_arr(1+110*FNum);
            MI(PUEL_ind,1)=frame_data_arr(37*2+110*FNum+PUEL(PUEL_ind));
        end
        for PUEL_ind=2:frame_data_arr(37+nr+1+110*FNum)+1
            EVM(1,1)=0;
            EVM(PUEL_ind,1)=frame_data_arr(111+110*FNum+PUEL(PUEL_ind));
        end
        % F_data=[PUEL,MI,EVM];
    end
    N_p_ind=1;
    for FNum=0:9
        N_p(N_p_ind)=frame_data_arr(37+nr+1+110*FNum);                      %Reading the number of paired UEs for all subframes.
        frame_data_arr(37+nr+1+110*FNum);
        frame_data_arr(1+110*FNum);                                         %Subframe number
        N_p_ind=N_p_ind+1;
    end
    N_p_all=[N_p_all,N_p];
end
C_RT_cMIMO{cidx}=2*N_p_all;                                                 %Computing sum spectral efficiency
cidx=cidx+1;
end
Sum_spectral_efficiency_CDF(C_RT_cMIMO);                                    %Plotting the CDF of sum spectral efficiency for different scenarios with correlation threshold 0.3
%% Downlink received constellation diagrams
REVM_indoor=DL_received_constellation(paired_UE_list_indoor,1,'qpsk');      %Downlink received constellation diagrams for all paired UEs decoded with QPSK in indoor scenario
disp('REVM for all paired UEs (Indoor)');                           
disp(["UE", "REVM"]);
disp(REVM_indoor);                                                          %Printing the REVM values for all paired UEs
REVM_outdoor=DL_received_constellation(paired_UE_list_outdoor,2,'qpsk');    %Downlink received constellation diagrams for all paired UEs decoded with QPSK in outdoor scenario
disp('REVM for all paired UEs (Outdoor)');                           
disp(["UE", "REVM"]);
disp(REVM_outdoor);                                                         %Printing the REVM values for all paired UEs
%% Function definition for Plotting Normalized Eigenvalues for different correlation thresholds
function []=Normalized_Eigenvalue_CDF(n_e_val_Th) 

n_e_val_48x36_Th_0_3_indoor=n_e_val_Th{1};
n_e_val_48x36_Th_0_3_outdoor=n_e_val_Th{2};

set(0,'defaultAxesFontSize',14)
set(0,'defaultAxesFontName','Times')

figure
h1=cdfplot(n_e_val_48x36_Th_0_3_indoor);
set(h1,'Marker','o','MarkerIndices',[1:100:1000]) %,'Color','r'
hold on
h2=cdfplot(n_e_val_48x36_Th_0_3_outdoor);
set(h2,'Marker','s','MarkerIndices',[1:100:800]) %,'Color','k'

xlabel('Normalized eigenvalues');
ylabel('CDF');
title('');
legend("Indoor: Th=0.3","Outdoor: Th=0.3",'Location','northwest');
end
%% Function definition for Plotting Sum spectral efficiency for different correlation thresholds
function []=Sum_spectral_efficiency_CDF(C_RT_cMIMO)

C_RT_cMIMO_I=C_RT_cMIMO{1};
C_RT_cMIMO_O=C_RT_cMIMO{2};

set(0,'defaultAxesFontSize',14)
set(0,'defaultAxesFontName','Times')

figure
h1=cdfplot(C_RT_cMIMO_I);
set(h1,'Color','b','Marker','o','MarkerIndices',[1:2:18]) %,'Color','red'
hold on
h2=cdfplot(C_RT_cMIMO_O);
set(h2,'Color','r','Marker','s','MarkerIndices',[1:1:10]) %,'Color','red'

xlabel('Sum spectral efficiency (bits/sec/Hz)');
ylabel('CDF');
title('');
legend("Indoor","Outdoor",'Location','southeast');
axis([0,30,0,1])
end
%% Function definition for Downlink received constellation diagrams
function [REVM_UE]=DL_received_constellation(paired_UE_list,scenario,modulation)
% This script reads downlink constellation data from binary files.
% Prints REVM for each paired UE, and plots received constellations for QPSK, 16-QAM, and 64-QAM under different threshold settings recorded from testbed.
% Inputs: paired_UE_list- list of paired UEs recorded from a experiment
%         scenario- 1 for indoor and 2 for outdoor
%         modulation- modulation scheme considered in the experiment
% Output: Plots constellation diagram for all paired UEs and returns their
%         corresponding REVM values.
figure;
for n=1:length(paired_UE_list)
    if(scenario==1)
        file_name = strcat('dl_received_constellation/indoor/ue',int2str(paired_UE_list(n)),'.bkp');
    elseif(scenario==2)
        file_name = strcat('dl_received_constellation/outdoor/ue',int2str(paired_UE_list(n)),'.bkp');
    end
ue_id=fopen(file_name,'r');
ue_file_arr=fread(ue_id,inf,'uint32');
fclose(ue_id);
ue_bin_arr=de2bi(ue_file_arr,32,'left-msb');
ue_bin_big_endian_arr=little_to_big_endian_arr(ue_bin_arr);
ue=bi_to_sp_float_arr(ue_bin_big_endian_arr);

sym=10;  %Total number of ofdm symbols +1

ue_1=(ue(1:(48*(sym-1))));
ue_cc=(ue_1(1:2:(48*(sym-1)))+1i*ue_1(2:2:(48*(sym-1))));
ue_c=ue_cc(49:end);

REVM(n)=REVM_computation(ue_c);

set(0,'defaultAxesFontSize',14)
set(0,'defaultAxesFontName','Times')
xstring='InPhase';
ystring='QPhase';

if(strcmp(modulation, 'qpsk'))
    if(scenario==1)
        subplot(3,5,n);
    elseif(scenario==2)
        subplot(3,4,n);
    end
    scatter(real(ue_c)*sqrt(2),imag(ue_c)*sqrt(2));
    title(['UE',int2str(paired_UE_list(n))]);
    grid on;
    xlabel(xstring);
    ylabel(ystring);
    axis([-2 2 -2 2])
    if(scenario==1)
        sgtitle('Indoor (Th=0.3)');
    elseif(scenario==2)
        sgtitle('Outdoor (Th=0.3)');
    end
    box
elseif(strcmp(modulation, '16qam'))
    if(scenario==1)
        subplot(3,5,n);
    elseif(scenario==2)
        subplot(3,4,n);
    end
    scatter(real(ue_c)*sqrt(10),imag(ue_c)*sqrt(10));
    title(['UE',int2str(paired_UE_list(n))]);
    grid on;
    xlabel(xstring);
    ylabel(ystring);
    axis([-4 4 -4 4])
    if(scenario==1)
        sgtitle('Indoor (Th=0.3)');
    elseif(scenario==2)
        sgtitle('Outdoor (Th=0.3)');
    end
    box
elseif(strcmp(modulation, '64qam'))
    if(scenario==1)
        subplot(3,5,n);
    elseif(scenario==2)
        subplot(3,4,n);
    end
    scatter(real(ue_c)*sqrt(42),imag(ue_c)*sqrt(42));
    title(['UE',int2str(paired_UE_list(n))]);
    grid on;
    xlabel(xstring);
    ylabel(ystring);
    axis([-10 10 -10 10])
    if(scenario==1)
        sgtitle('Indoor (Th=0.3)');
    elseif(scenario==2)
        sgtitle('Outdoor (Th=0.3)');
    end
    box
end
end
REVM_UE=[paired_UE_list.' REVM.'];
end
%% REVM Computation
function [EVM_ue]=REVM_computation(ue_c)
sym_num=9;      %number of ofdm symbols for data
sym_num_1=7;    %number of ofdm symbols for which EVM is calculated

%Actual QPSK transmitted data is denoted as ap1_9_sym, which will be used for REVM calculation
ap1_9_sym=[0.7071,0.7071,0.7071,-0.7071,-0.7071,0.7071,0.7071,0.7071,-0.7071,-0.7071,-0.7071,-0.7071,0.7071,-0.7071,-0.7071,-0.7071,-0.7071,-0.7071,-0.7071,0.7071,0.7071,-0.7071,0.7071,-0.7071,0.7071,-0.7071,-0.7071,0.7071,0.7071,-0.7071,-0.7071,-0.7071,-0.7071,-0.7071,0.7071,0.7071,-0.7071,-0.7071,0.7071,-0.7071,-0.7071,0.7071,0.7071,0.7071,0.7071,0.7071,0.7071,-0.7071,0.7071,-0.7071,0.7071,0.7071,-0.7071,0.7071,-0.7071,0.7071,0.7071,-0.7071,0.7071,-0.7071,0.7071,-0.7071,0.7071,-0.7071,-0.7071,-0.7071,-0.7071,0.7071,-0.7071,-0.7071,-0.7071,0.7071,-0.7071,-0.7071,0.7071,0.7071,0.7071,0.7071,0.7071,-0.7071,0.7071,-0.7071,0.7071,-0.7071,-0.7071,-0.7071,-0.7071,0.7071,-0.7071,0.7071,-0.7071,-0.7071,0.7071,-0.7071,0.7071,-0.7071,-0.7071,-0.7071,-0.7071,0.7071,-0.7071,-0.7071,0.7071,0.7071,-0.7071,0.7071,0.7071,0.7071,0.7071,-0.7071,-0.7071,0.7071,-0.7071,0.7071,-0.7071,-0.7071,-0.7071,0.7071,0.7071,-0.7071,-0.7071,0.7071,-0.7071,-0.7071,-0.7071,0.7071,0.7071,0.7071,-0.7071,-0.7071,-0.7071,0.7071,0.7071,-0.7071,-0.7071,0.7071,0.7071,0.7071,-0.7071,-0.7071,0.7071,-0.7071,-0.7071,-0.7071,0.7071,-0.7071,0.7071,-0.7071,-0.7071,-0.7071,0.7071,-0.7071,-0.7071,-0.7071,0.7071,-0.7071,0.7071,0.7071,-0.7071,0.7071,-0.7071,0.7071,0.7071,0.7071,-0.7071,0.7071,0.7071,0.7071,-0.7071,0.7071,0.7071,-0.7071,0.7071,-0.7071,0.7071,0.7071,0.7071,0.7071,0.7071,-0.7071,-0.7071,-0.7071,0.7071,-0.7071,0.7071,-0.7071,0.7071,0.7071,0.7071,0.7071,-0.7071,0.7071,-0.7071,-0.7071,0.7071,0.7071,-0.7071,-0.7071,-0.7071,0.7071,0.7071,0.7071,-0.7071,-0.7071,0.7071,-0.7071,-0.7071,-0.7071,-0.7071,-0.7071,0.7071,0.7071,-0.7071,-0.7071,-0.7071,0.7071,-0.7071,-0.7071,-0.7071,0.7071,-0.7071,0.7071,-0.7071,0.7071,0.7071,0.7071,-0.7071,-0.7071,0.7071,0.7071,-0.7071,0.7071,0.7071,0.7071,0.7071,-0.7071,0.7071,-0.7071,0.7071,0.7071,-0.7071,0.7071,0.7071,0.7071,-0.7071,-0.7071,0.7071,0.7071,-0.7071,-0.7071,0.7071,0.7071,0.7071,0.7071,0.7071,-0.7071,0.7071,-0.7071,-0.7071,-0.7071,-0.7071,-0.7071,-0.7071,-0.7071,0.7071,0.7071,0.7071,0.7071,0.7071,0.7071,-0.7071,0.7071,0.7071,-0.7071,-0.7071,-0.7071,-0.7071,-0.7071,-0.7071,0.7071,-0.7071,0.7071,-0.7071,-0.7071,0.7071,-0.7071,-0.7071,0.7071,0.7071,0.7071,-0.7071,-0.7071,-0.7071,0.7071,-0.7071,0.7071,0.7071,-0.7071,0.7071,0.7071,-0.7071,-0.7071,0.7071,0.7071,0.7071,0.7071,0.7071,-0.7071,0.7071,0.7071,-0.7071,0.7071,-0.7071,0.7071,0.7071,-0.7071,0.7071,-0.7071,-0.7071,0.7071,-0.7071,0.7071,0.7071,-0.7071,-0.7071,-0.7071,-0.7071,0.7071,-0.7071,-0.7071,0.7071,0.7071,-0.7071,0.7071,0.7071,0.7071,0.7071,0.7071,0.7071,-0.7071,-0.7071,-0.7071,0.7071,-0.7071,0.7071,-0.7071,0.7071,0.7071,-0.7071,-0.7071,0.7071,0.7071,0.7071,0.7071,0.7071,0.7071,-0.7071,0.7071,0.7071,0.7071,-0.7071,0.7071,0.7071,-0.7071,-0.7071,0.7071,0.7071,0.7071,0.7071,0.7071,0.7071,-0.7071,-0.7071,-0.7071,-0.7071,-0.7071,0.7071,0.7071,0.7071,-0.7071,0.7071,0.7071,0.7071,-0.7071,-0.7071,0.7071,-0.7071,0.7071,0.7071,-0.7071,-0.7071,0.7071,0.7071,0.7071,0.7071,0.7071,0.7071,-0.7071,-0.7071,-0.7071,-0.7071,0.7071,0.7071,0.7071,0.7071,0.7071,-0.7071,-0.7071,-0.7071,-0.7071,0.7071,0.7071,0.7071,0.7071,0.7071,-0.7071,0.7071,0.7071,-0.7071,-0.7071,-0.7071,-0.7071,-0.7071,-0.7071,0.7071,0.7071,0.7071,-0.7071,-0.7071,-0.7071,0.7071,-0.7071];
ap1=ap1_9_sym(1:24*sym_num*2);
ap1_c=ap1(1:2:end)+1i*ap1(2:2:end);
ap1_rearranged = reshape(ap1_c,sym_num,24);
ap1_array = reshape(ap1_rearranged.',1,sym_num*24).';
ap1_Tx = ap1_array(1:sym_num_1*24);

ap_Tx(:,1)=ap1_Tx;

EVM_ue = sqrt(sum(abs(ue_c-ap_Tx(:,1)).^2)/length(ap1_Tx))*100;
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