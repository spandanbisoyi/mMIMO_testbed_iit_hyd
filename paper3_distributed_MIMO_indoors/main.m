clc
clear
close all

nt=48;nr=17;
N=54;
n_e_val=[];
for k=1:N
%% Reading binary files for recorded CSI from testbed
file_name = strcat('channel_data/ch',int2str(k),'.bkp');                    %Locating the binary file recorded from testbed for CSI matrix (The file 'ch0.bkp' contains the 36×48 CSI matrix stored in array format, along with the paired UE list and the number of paired UEs.)
r_ch_id=fopen(file_name,'r');                                               %Opening the binary file corresponding to above location
r_ch_file=fread(r_ch_id,inf,'uint32');                                      %Reading the binary file of CSI matrix in unsigned int format
fclose(r_ch_id);                                                            %Closing the binary file
r_ch_bin_arr=de2bi(r_ch_file,32,'left-msb');                                %Converting the decimal datas to binary data
r_ch_bin_big_endian_arr=little_to_big_endian_arr(r_ch_bin_arr);             %Performing little to big endian conversion
r_ch_arr=bi_to_sp_float_arr(r_ch_bin_big_endian_arr);                       %Converting big endian binary to single precision float array, where the complex array is stored as interleaved real and imaginary samples.
r_ch_c_arr=(r_ch_arr(1:2:2*nr*nt)+1i*r_ch_arr(2:2:2*nr*nt));                %Converting the real array to complex array
H_ch=reshape(r_ch_c_arr,nt,nr).';                                           %Rearranging the elements of the complex array in matrix form to get the actual CSI matrix for kth subframe
H=H_ch;                                                                     %CSI matrix of kth subframe
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
%% Eigen analysis
e_val_kth_channel=eig(H_norm*H_norm');                                      %Finding the eigen values of CSI matrix
p=size(H_norm*H_norm');
ein=p(1);
n_e_val_kth_channel=e_val_kth_channel/max(e_val_kth_channel);               %Computing normalized eigen values
n_e_val=[n_e_val;n_e_val_kth_channel];                                      %Normalized eigen values for 48 subframes

prod(sqrt(e_val_kth_channel));
end
UE_list=[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17];                        %UE_list denotes the list of all UEs in the experiment
%% Plotting Normalized Eigenvalues for different correlation thresholds
Normalized_Eigenvalue_CDF(n_e_val);                                         %Plotting the CDF of normalized eigen values
%% Plotting Sum spectral efficiency for different correlation thresholds
Sum_spectral_efficiency_CDF();                                              %Plotting the CDF of sum spectral efficiency
%% Downlink received constellation diagrams
REVM=DL_received_constellation(UE_list);                                    %Downlink received constellation diagrams for all paired UEs decoded with QPSK for correlation threshold=0.4
disp('REVM and modulation scheme for all UEs');                           
disp(["UE", "modulation", "REVM"]);
disp(REVM);                                                                 %Printing the REVM values for all paired UEs
%% Function definition for Plotting Normalized Eigenvalues
function []=Normalized_Eigenvalue_CDF(n_e_val) 

set(0,'defaultAxesFontSize',14)
set(0,'defaultAxesFontName','Times')

figure
h=cdfplot(n_e_val);
set(h,'Marker','o','MarkerIndices',[1:200:1800]) 
xlabel('Normalized eigenvalues');
ylabel('CDF');
title('');
legend("dMIMO",'Location','southeast');
end
%% Function definition for Plotting Sum spectral efficiency
function []=Sum_spectral_efficiency_CDF()
N=47;
C_RT_dMIMO=[];
for k=0:N
    file_name = strcat('frame_data/frame_data',int2str(k),'.bkp');
    frame_data_id=fopen(file_name,'r');
    frame_data_file=fread(frame_data_id,inf,'uint32');
    fclose(frame_data_id);
    frame_data_bin_arr=de2bi(frame_data_file,32,'left-msb');
    frame_data_bin_big_endian_arr=little_to_big_endian_arr(frame_data_bin_arr);
    frame_data_arr=bi_to_sp_float_arr(frame_data_bin_big_endian_arr);
    %Reading datas for all subframes 
    for FNum=0:8
        %Reading the paired UE list recorded in frame_data_arr for multiple subframes. FNum is the subframe index. frame_data_arr(56+110*FNum) is number of paired UEs.
        PUEL=[frame_data_arr(56+110*FNum);frame_data_arr(38+110*FNum:38+110*FNum+frame_data_arr(56+110*FNum)-1)]; 
        %Reading the modulation index for all UEs across multiple subframes. frame_data_arr(1+110*FNum) is the subframe number
        for PUEL_ind=2:frame_data_arr(56+110*FNum)+1
            MI(1,1)=frame_data_arr(1+110*FNum);
            MI(PUEL_ind,1)=frame_data_arr(37*2+110*FNum+PUEL(PUEL_ind));
        end
        for PUEL_ind=2:frame_data_arr(56+110*FNum)+1
            EVM(1,1)=0;
            EVM(PUEL_ind,1)=frame_data_arr(111+110*FNum+2*PUEL(PUEL_ind)-1);
            % 2*PUEL(PUEL_ind)-1
        end
        % disp(strcat('SF',int2str(FNum)))
        F_data=[PUEL,MI,EVM];
        count_QPSK=0;
        count_16QAM=0;
        count_64QAM=0;
        for MI_ind=2:frame_data_arr(56+110*FNum)
            if(MI(MI_ind)==1)
                count_QPSK=count_QPSK+1;
            elseif(MI(MI_ind)==2)
                count_16QAM=count_16QAM+1;
            elseif(MI(MI_ind)==3)
                count_64QAM=count_64QAM+1;
            end
        end
        C_RT_kth_SF=2*count_QPSK+4*count_16QAM+6*count_64QAM;               %Computing the sum spectral efficiency for a subframe
        C_RT_dMIMO=[C_RT_dMIMO;C_RT_kth_SF];                                %Sum spectral efficiency for all subframes
    end
end

set(0,'defaultAxesFontSize',14)
set(0,'defaultAxesFontName','Times')

figure
h=cdfplot(C_RT_dMIMO);
set(h,'Color','b','Marker','o','MarkerIndices',[1:3:48]) %,'Color','red'
xlabel('Sum Spectral Efficiency (bits/sec/Hz)');
ylabel('CDF');
title('');
legend("dMIMO",'Location','southeast');
end
%% Function definition for Downlink received constellation diagrams
function [REVM_UE]=DL_received_constellation(UE_list)
% This script reads downlink constellation data from binary files.
% Prints REVM for each paired UE, and plots received constellations for QPSK, 16-QAM, and 64-QAM under different threshold settings recorded from testbed.
% Inputs: UE_list- list of all UEs
% Output: Plots constellation diagram for all UEs and returns their
%         corresponding REVM values and modulation schemes.
np=length(UE_list);
MI=[1,2,2,2,3,2,2,2,2,2,2,2,1,1,1,1,2];
figure;
for UE_index=1:np
file_name = strcat('dl_received_constellation/ue',int2str(UE_list(UE_index)),'.bkp');
ue_id=fopen(file_name,'r');
ue_file_arr=fread(ue_id,inf,'uint32');
fclose(ue_id);
ue_bin_arr=de2bi(ue_file_arr,32,'left-msb');
ue_bin_big_endian_arr=little_to_big_endian_arr(ue_bin_arr);
ue=bi_to_sp_float_arr(ue_bin_big_endian_arr);

% REVM(UE_index)=ue(385);
sym=7;  %Total number of ofdm symbols +1

ue_1=(ue(1:(48*(sym-1))));
ue_cc=(ue_1(1:2:(48*(sym-1)))+1i*ue_1(2:2:(48*(sym-1))));
ue_c=ue_cc(49:end);

set(0,'defaultAxesFontSize',14)
set(0,'defaultAxesFontName','Times')
xstring='InPhase';
ystring='QPhase';

if(UE_index<=9)
subplot(3,3,UE_index);
if(MI(UE_index)==1)
    modulation(UE_index)="qpsk";
    REVM(UE_index)=REVM_computation(ue_c,modulation(UE_index));
    scatter(real(ue_c)*sqrt(2),imag(ue_c)*sqrt(2));
    title(['UE',int2str(UE_list(UE_index))]);
    grid on;
    xlabel(xstring);
    ylabel(ystring);
    axis([-2 2 -2 2])
    box
elseif(MI(UE_index)==2)
    modulation(UE_index)="16qam";
    REVM(UE_index)=REVM_computation(ue_c,modulation(UE_index));
    scatter(real(ue_c)*sqrt(10),imag(ue_c)*sqrt(10));
    title(['UE',int2str(UE_list(UE_index))]);
    grid on;
    xlabel(xstring);
    ylabel(ystring);
    axis([-4 4 -4 4])
elseif(MI(UE_index)==3)
    modulation(UE_index)="64qam";
    REVM(UE_index)=REVM_computation(ue_c,modulation(UE_index));
    scatter(real(ue_c)*sqrt(42),imag(ue_c)*sqrt(42));
    title(['UE',int2str(UE_list(UE_index))]);
    grid on;
    xlabel(xstring);
    ylabel(ystring);
    axis([-10 10 -10 10])
end
elseif(UE_index>9 && UE_index<=18)
    if(UE_index==10)
        figure;
    end
    subplot(3,3,UE_index-9);
    if(MI(UE_index)==1)
        modulation(UE_index)="qpsk";
        REVM(UE_index)=REVM_computation(ue_c,modulation(UE_index));
        scatter(real(ue_c)*sqrt(2),imag(ue_c)*sqrt(2));
        title(['UE',int2str(UE_list(UE_index))]);
        grid on;
        xlabel(xstring);
        ylabel(ystring);
        axis([-2 2 -2 2])
        box
    elseif(MI(UE_index)==2)
        modulation(UE_index)="16qam";
        REVM(UE_index)=REVM_computation(ue_c,modulation(UE_index));
        scatter(real(ue_c)*sqrt(10),imag(ue_c)*sqrt(10));
        title(['UE',int2str(UE_list(UE_index))]);
        grid on;
        xlabel(xstring);
        ylabel(ystring);
        axis([-4 4 -4 4])
    elseif(MI(UE_index)==3)
        modulation(UE_index)="64qam";
        REVM(UE_index)=REVM_computation(ue_c,modulation(UE_index));
        scatter(real(ue_c)*sqrt(42),imag(ue_c)*sqrt(42));
        title(['UE',int2str(UE_list(UE_index))]);
        grid on;
        xlabel(xstring);
        ylabel(ystring);
        axis([-10 10 -10 10])
    end
end
end
REVM_UE=[UE_list.' modulation.' REVM.'];
end
%% REVM Computation
function [EVM_ue]=REVM_computation(ue_c,modulation)
sym_num=6;      %number of ofdm symbols for data
sym_num_1=4;    %number of ofdm symbols for which EVM is calculated

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