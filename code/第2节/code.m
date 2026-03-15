function [MIMO_Channel,ArrayResponse_TX,ArrayResponse_RX,Alpha] = ChannelGenereationMIMO(Nt,Nr,NumCluster,NumRay,AS)
%%
% Input : 
%         Nt            : Number of transmit antennas 
%         Nr            : Number of receive antennas
%         NumCluster    : Number of clusters
%         NumRay        : Number of rays per cluster
%         AS            : Fixed angular spread at both the transmitter and receiver
% Output : 
%         MIMO_Chan         : Nr x Nt, MIMO channel matrix
%         ArrayResponse_TX  : Nt x (NumCluster x NumRay), transmit steering vectors
%         ArrayResponse_RX  : Nr x (NumCluster x NumRay), receive steering vectors
%         Alpha             : 1 x (NumCluster x NumRay), complex path gain
%
%% Generate random mean AOD, EOD, AOA, EOA of each cluster

% Azimuth   : [-180,180] degree
% Elevation : [0,180] degree

% Transmitter(sectorized)
% Azimuth   : 60 degree wide w.r.t 0 degree
% Elevation : 20 degree wide w.r.t 90 degree
minAOD = -30;
maxAOD = 30; 
minEOD = 80;
maxEOD = 100;
Cluster_AOD = rand(NumCluster,1)*(maxAOD-minAOD)+minAOD; % [-30,30] degree
Cluster_EOD = rand(NumCluster,1)*(maxEOD-minEOD)+minEOD; % [80,100] degree

% Receiver(omni-directional)
Cluster_AOA = (rand(NumCluster,1)-0.5)*360; % [-180,180] degree
Cluster_EOA = rand(NumCluster,1)*180; % [0,180] degree


%% Generate random AOD, EOD, AOA, EOA of rays per cluster (Laplacian distribution)

b = AS/sqrt(2); % Scaling parameter, degree

Randomness = rand(NumRay*NumCluster,1)-0.5;

% Dimension of AOD, EOD, AOA, EOA : (NumRay x NumCluster) x 1 
Ray_AOD = repelem(Cluster_AOD,NumRay,1)-b*sign(Randomness).*log(1-2.*abs(Randomness));
Ray_EOD = repelem(Cluster_EOD,NumRay,1)-b*sign(Randomness).*log(1-2.*abs(Randomness));
Ray_AOA = repelem(Cluster_AOA,NumRay,1)-b*sign(Randomness).*log(1-2.*abs(Randomness));
Ray_EOA = repelem(Cluster_EOA,NumRay,1)-b*sign(Randomness).*log(1-2.*abs(Randomness));

%% Obtain antenna element position vectors (normalized by half of the wavelength)

% Transmitter
Nt_H = sqrt(Nt);
Nt_V = sqrt(Nt);

X_Tx = zeros(1,Nt);
[Y_Tx,Z_Tx] = meshgrid(0:Nt_H-1,0:Nt_V-1);
TxPos = [X_Tx;Y_Tx(:).';Z_Tx(:).']; % 3 x Nt

% Receiver
Nr_H = sqrt(Nr);
Nr_V = sqrt(Nr);

X_Rx = zeros(1,Nr);
[Y_Rx,Z_Rx] = meshgrid(0:Nr_H-1,0:Nr_V-1);
RxPos = [X_Rx;Y_Rx(:).';Z_Rx(:).']; % 3 x Nr

%% Obtain array response vectors at the transmitter and receiver

SphericalUnitVecTx = getSphericalUnitVector(Ray_EOD,Ray_AOD); % 3 x NumRay*NumCluster
SphericalUnitVecRx = getSphericalUnitVector(Ray_EOA,Ray_AOA); % 3 x NumRay*NumCluster

ArrayResponse_TX = (1/sqrt(Nt))*exp(1i*pi*TxPos.'*SphericalUnitVecTx); % Nt x NumRay*NumCluster
ArrayResponse_RX = (1/sqrt(Nr))*exp(1i*pi*RxPos.'*SphericalUnitVecRx); % Nr x NumRay*NumCluster

%% Generate complex path gain

Alpha = sqrt(1/2)*(randn(1,NumRay*NumCluster)+1i*randn(1,NumRay*NumCluster));

%% Generate MIMO channel matrix 

MIMO_Channel = sqrt((Nt*Nr)/(NumRay*NumCluster))*ArrayResponse_RX*diag(Alpha)*ArrayResponse_TX';

end

function SphericalUnitVector = getSphericalUnitVector(theta,phi)
%% 
% Input:
%       theta, phi : M x 1
% Output:
%       SphericalUnitVector: 3 x M
% 
SphericalUnitVector = [(sind(theta).*cosd(phi)).';...
                       (sind(theta).*sind(phi)).';...
                       (cosd(theta)).']; 
end

function [IndexAt,IndexAr] = findSteeringVector(H,At,Ar,Ns)
%%% This is the function to obtain the beam steering vectors at the transmitter and receiver 
%%% by finding the array response vectors corresponding to the largest effective channel gain (an exhaustive search).
%%% Only suitable for single data stream case.
% 
% Input
%       H     : channel matrix, Nr x Nt
%       At    : the collection of transmit steering vectors, Nt x NumPath
%       Ar    : the collection of receive steering vectors, Nr x NumPath
% Output
%      IndexAt : index of the path selected at the transmitter
%      IndexAr : index of the path selected at the receiver
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Nt = size(At,1);
Nr = size(Ar,1);

NumPath = size(At,2);

EffGain = zeros(NumPath);

if Ns == 1    
    
    for t = 1:NumPath
        
        for r = 1:NumPath
        
            EffGain(t,r) = abs(Ar(:,r)'*H*At(:,t));
        
        end
    
    end
else
    
    error('Not single data stream.');
    
end


[row,col] = find(EffGain == max(max(EffGain)));

IndexAt = row;
IndexAr = col;
      
end

function [SpectralEfficiency_Optimal,SpectralEfficiency_Hybrid,SpectralEfficiency_BeamSteering] = SpatiallySparsePrecoding(Nt,Nr,Ns,NumRF,NumCluster,NumRay,AS,SNR_dB,ITER)
%% 
% Input:
%       Nt          :   Number of transmit antennas
%       Nr          :   Number of receive antennas
%       Ns          :   Number of transmitted data stream
%       NumRF       :   Number of RF chains at both the transmitter and receiver
%       NumCluster  :   Number of scattering clusters
%       NumRay      :   Number of rays per cluster
%       AS          :   Angular spread in the cluster
%       SNR_dB      :   SNR values in decibel
%       ITER        :   Number of random channel generations
% Output:
%       SpectralEfficiency_Optimal      :    Spectral efficiency in bits/s/Hz, obtained by optimal unconstrained precoding 
%       SpectralEfficiency_Hybrid       :   Spectral efficiency in bits/s/Hz, obtained by hybrid precoding 
%       SpectralEfficiency_BeamSteering :   Spectral efficiency in bits/s/Hz, obtained by beam steering 
%
%% Initialize variables

% Initialize spectral efficiency (bits/s/Hz) obtained by 
% A. optimal unconstrained precoding
% B. hybrid precoding and combining
% C. beam steering

SpectralEfficiency_Optimal = zeros(ITER,length(SNR_dB)); 
SpectralEfficiency_Hybrid = zeros(ITER,length(SNR_dB)); 
SpectralEfficiency_BeamSteering = zeros(ITER,length(SNR_dB)); 

%% Obtain simulation parameters

% SNR in linear scale
SNR = 10.^(SNR_dB./10);

%% Precoding and combining algorithms

for s = 1:length(SNR) % Loop over SNR
           
    for i = 1:ITER
                
        % Generate MIMO channel matrix, array response vectors at the transmitter and receiver and complex path gain
        % H     :   Nr x Nt
        % At    :   Nt x (NumRay x NumCluster)
        % Ar    :   Nr x (NumRay x NumCluster)
        % Alpha :   1 x (NumRay x NumCluster)
        
        [H,At,Ar,Alpha] = ChannelGenereationMIMO(Nt,Nr,NumCluster,NumRay,AS);
               
        %%%% A. Optimal unconstrained precoding and combining
        
        [U,S,V] = svd(H);
        F_Opt = V(:,1:Ns); % Nt x Ns, optimal precoder
        W_Opt = ((1/sqrt(SNR(s)))*(F_Opt'*H'*H*F_Opt+Ns/SNR(s)*eye(Ns))\(F_Opt'*H'))'; % Nr x Ns, optimal combiner

        %%%% B. Hybrid sparse precoding and combining via orthogonal matching pursuit (OMP)
        
        % Algorithm 1 - Hybrid precoding        
        F_RF = [];
        F_Res = F_Opt;
        
        for r = 1:NumRF             
            Psi = At'*F_Res; % Step 4, the projection of all paths on the optimal precoder
            [~,k] = max(diag(Psi*Psi')); % Step 5, select the path that has the largest projection
            F_RF = [F_RF At(:,k)]; % Step 6, append the selected vector to the RF precoder
            F_BB = (F_RF'*F_RF)\(F_RF'*F_Opt); % Step7, baseband precoder calculated by least squares solution
            F_Res = (F_Opt-F_RF*F_BB)/norm(F_Opt-F_RF*F_BB,'fro'); % Step 8, remove the contribution of the selected vector from F_res        
        end       
        F_BB = sqrt(Ns)*(F_BB/norm(F_RF*F_BB,'fro')); % Step 10, normalize F_BB to meet total power constraint
        
        % Algorithm 2 - Hybrid combining
        CovRx = (SNR(s)/Ns)*H*F_RF*F_BB*F_BB'*F_RF'*H'+eye(Nr); % Nr x Nr, covariance matrix of received signals at receive antennas
        W_MMSE = ((1/sqrt(SNR(s)))*(F_BB'*F_RF'*H'*H*F_RF*F_BB+(Ns/SNR(s))*eye(Ns))\(F_BB'*F_RF'*H'))';
                
        W_RF = [];
        W_Res = W_MMSE;
        
        for r = 1:NumRF
            Psi = Ar'*CovRx*W_Res; 
            [~,k] = max(diag(Psi*Psi'));
            W_RF = [W_RF Ar(:,k)];
            W_BB = (W_RF'*CovRx*W_RF)\(W_RF'*CovRx*W_MMSE);
            W_Res = (W_MMSE-W_RF*W_BB)/norm(W_MMSE-W_RF*W_BB,'fro');
        end
        
        %%%% C. Beam steering based on the highest effective channel gain
        
        [IndexAt,IndexAr] = findSteeringVector(H,At,Ar,Ns);
        
        F_BS = [];
        W_BS = [];
        
        for n = 1:Ns 
            F_BS = [F_BS At(:,IndexAt(n))];
            W_BS = [W_BS Ar(:,IndexAr(n))];
        end

        %%%% Spectial efficiency calculation 
        
        % A. Optimal solution
        Rn_Opt = W_Opt'*W_Opt;
        SpectralEfficiency_Optimal(i,s) = abs(log2(det(eye(Ns)+(SNR(s)/Ns)*(Rn_Opt\(W_Opt'*H*F_Opt*F_Opt'*H'*W_Opt)))));
        
        % B. Hybrid solution       
        Rn_Hybrid = W_BB'*W_RF'*W_RF*W_BB;        
        SpectralEfficiency_Hybrid(i,s) = abs(log2(det(eye(Ns)+(SNR(s)/Ns)*(Rn_Hybrid\(W_BB'*W_RF'*H*F_RF*F_BB*F_BB'*F_RF'*H'*W_RF*W_BB)))));
        
        % C. Beam steering solution
        Rn_BS = W_BS'*W_BS;
        SpectralEfficiency_BeamSteering(i,s) = abs(log2(det(eye(Ns)+(SNR(s)/Ns)*(Rn_BS\(W_BS'*H*F_BS*F_BS'*H'*W_BS)))));

    end
    
end

SpectralEfficiency_Optimal = mean(SpectralEfficiency_Optimal,1); % 1 x length(SNR)
SpectralEfficiency_Hybrid = mean(SpectralEfficiency_Hybrid,1);
SpectralEfficiency_BeamSteering = mean(SpectralEfficiency_BeamSteering,1); 
end

%%
% This is the script to 
%                       1. perform hybrid precoding algorithms for
%                       millimeter wave MIMO systems [1];
%                       2. calculate the spectral efficiency obtained by
%                       optimal precoding, hybrid precoding and beam steering for comparison.
% Basic assumptions:
%                       1. The millimeter wave channel is modelled as a narrowband clustered channnel.  
%                       2. The hybrid precoding algorithm is based on the spatially sparse nature of
%                          millimeter wave propagation.
%                       3. This is the code solely for plotting Fig.3 and 4 of single data stream (Ns = 1)
%                          case in the paper [1] 
%                           
% 
% [1] O. E. Ayach et al., ?Spatially sparse precoding in millimeter wave MIMO systems,? Mar. 2013.
%
%% Clear workplace

clear variables;
close all;
clc;

%% System parameters

Nt = [64 256]; % Number of transmit antennas
Nr = [16 64]; % Number of receive antennas

Ns = 1; % Number of data streams (修正：改为标量)
NumRF = [4 6]; % Number of RF chains for precoding and combining 

NumCluster = 8; % Number of clusters
NumRay = 10; % Number of rays per cluster

AS = 7.5; % Angular spread of 7.5 degree

SNR_dB = -40:5:0; % Range of SNR in dB

ITER = 10; % Number of channel generations (为了快速测试，减少迭代次数)

%% Fig.3 
fprintf('Running simulation for 64x16 MIMO system...\n');
[Fig3_SE_Optimal,Fig3_SE_Hybrid,Fig3_SE_BeamSteering] = SpatiallySparsePrecoding(Nt(1),Nr(1),Ns,NumRF(1),NumCluster,NumRay,AS,SNR_dB,ITER);

figure(); 
l1 = plot(SNR_dB,Fig3_SE_Optimal,'-s','Color',[0 0.5 0],'LineWidth',2.0,'MarkerSize',8.0);hold on;
l2 = plot(SNR_dB,Fig3_SE_Hybrid,'-o','Color',[0 0.45 0.74],'LineWidth',2.0,'MarkerSize',8.0);hold on;
l3 = plot(SNR_dB,Fig3_SE_BeamSteering,'-d','Color',[0.85 0.33 0.1],'LineWidth',2.0,'MarkerSize',8.0);hold on;
grid on;
legend([l1 l2 l3],'Optimal unconstrained precoding','Hybrid precoding and combining','Beam steering','Location','northwest');
title('Fig.3. 64 x 16 mmWave system with 4 RF chains for sparse precoding and MMSE combining (Ns = 1)');
xlabel('SNR (dB)'); ylabel('Spectral efficiency (bits/s/Hz)');

%% Fig.4
fprintf('Running simulation for 256x64 MIMO system...\n');
[Fig4_SE_Optimal,Fig4_SE_Hybrid,Fig4_SE_BeamSteering] = SpatiallySparsePrecoding(Nt(2),Nr(2),Ns,NumRF(2),NumCluster,NumRay,AS,SNR_dB,ITER);

figure(); 
l1 = plot(SNR_dB,Fig4_SE_Optimal,'-s','Color',[0 0.5 0],'LineWidth',2.0,'MarkerSize',8.0);hold on;
l2 = plot(SNR_dB,Fig4_SE_Hybrid,'-o','Color',[0 0.45 0.74],'LineWidth',2.0,'MarkerSize',8.0);hold on;
l3 = plot(SNR_dB,Fig4_SE_BeamSteering,'-d','Color',[0.85 0.33 0.1],'LineWidth',2.0,'MarkerSize',8.0);hold on;
grid on;
legend([l1 l2 l3],'Optimal unconstrained precoding','Hybrid precoding and combining','Beam steering','Location','northwest');
title('Fig.4. 256 x 64 mmWave system with 6 RF chains for sparse precoding and MMSE combining (Ns = 1)');
xlabel('SNR (dB)'); ylabel('Spectral efficiency (bits/s/Hz)');

fprintf('Simulation completed successfully!\n');