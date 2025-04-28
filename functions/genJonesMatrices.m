function [Hi,p] = genJonesMatrices(p)

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function that calculates the Jones matrix for each fiber segment
%
% Authors: 
% Original code by S. Guerrier, C. Dorize & E. Awwad - 2022
% Modified version by E. Awwad - 2024 elie.awwad@telecom-paris.fr
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 3) Jones matrices
% Dual-polarization transmission and reflection
resolution = p.fibre.spatialRes/p.fibre.polCorrL; 
if resolution >= 1 % polarization correlation length exceeded, random SOP in new segment
    rng(p.seeds.seed_rotPol);
    p.fibre.rotPolarTab = single(asin(sqrt(rand(p.fibre.nbSegments,1)))); %random rotation of phase from 0 to segment index i, [-pi, pi]
    rng(p.seeds.seed_beta);
    p.fibre.betaTab = single(2*pi*rand(p.fibre.nbSegments,1)-pi); %for phase retardance matrix (U = R*D, D = diag(exp(+-jbeta/2)) ), from 0 to segment index i

    rng(p.seeds.seed_gamma);
    p.fibre.gammaTab = single(2*pi*rand(p.fibre.nbSegments,1)-pi); %add a 3rd random parameter to cover the whole Poincaré sphere
    fprintf('\n *** Lsegment >= polarization beat length : sharp evolution with fiber segment to polarization beat length ratio r = %d *** \n', resolution);
else 
    %data = data_loading('car_filt_high_reduced.mat');
    rng(p.seeds.seed_rotPol); % polarization correlation length not exceeded, Wiener process
    p.fibre.evolThetaTab = [single(asin(sqrt(rand()))); single(asin(sqrt(rand(p.fibre.nbSegments,1)))*resolution)];
    rng(p.seeds.seed_beta);
    p.fibre.evolBetaTab =  single(pi*resolution*(1+rand(p.fibre.nbSegments, 1)*0.0001));%single(rand()*2*pi-pi) % single((rand(p.fibre.nbSegments,1)*2*pi-pi) %single(rand(p.fibre.nbSegments,1)*2*pi-pi)*0.00001 %[single(rand()*2*pi-pi); single((rand(p.fibre.nbSegments-1,1)*2*pi-pi)*resolution)];
    rng(p.seeds.seed_gamma);
    p.fibre.evolGammaTab = [single(rand()*2*pi-pi); single((rand(p.fibre.nbSegments-1,1)*2*pi-pi)*resolution)];
    p.fibre.rotPolarTab = cumsum(p.fibre.evolThetaTab);
    p.fibre.betaTab = unwrap(cumsum(p.fibre.evolBetaTab)); 
    p.fibre.gammaTab = unwrap(cumsum(p.fibre.evolGammaTab));
    fprintf('\n *** Lsegment < polarization beat length : smooth evolution with fiber segment to polarization beat length ratio r = %d *** \n', resolution); 
end

% Polar rotation at fiber entrance
rng(p.seeds.seed_theta); 
theta = 2*pi*rand()-pi;
p.fibre.R0 = [cos(theta) sin(theta) ; sin(theta) -cos(theta)];

% Create matrices
M_i = [sqrt(1-p.fibre.AlphaPolRay), -sqrt(p.fibre.AlphaPolRay); sqrt(p.fibre.AlphaPolRay), sqrt(1-p.fibre.AlphaPolRay)]; %mirror effect on polar in each fiber segment i. M is identity if alpha null

% Generate Jones matrices of each segment
U_i = single(NaN(2,2*p.fibre.nbSegments)); %Output array storing the 2*2 Jones matrix of each segment
for n=1:p.fibre.nbSegments 
    expjBeta = exp(1j*(p.fibre.evolBetaTab(n))); 
    sinRot = sin(p.fibre.evolThetaTab(n)); cosRot = cos(p.fibre.evolThetaTab(n));%sinRot = sin(0); cosRot = cos(0);%sinRot = sin(p.fibre.evolThetaTab(n)); cosRot = cos(p.fibre.evolThetaTab(n));%sinRot = sin(0); cosRot = cos(0);%sinRot = sin(p.fibre.evolThetaTab(n)); cosRot = cos(p.fibre.evolThetaTab(n));%sinRot = sin(0); cosRot = cos(0);%sinRot = sin(p.fibre.evolThetaTab(n)); cosRot = cos(p.fibre.evolThetaTab(n)); %sinRot = sin(20); cosRot = cos(20);%
    U_i(:,1+2*(n-1):2*n) = [cosRot sinRot ; -sinRot cosRot] * [expjBeta 0 ; 0 conj(expjBeta)]*[cosRot -sinRot ; sinRot cosRot];%Calculate matrix U_forw for the current fiber segment
end

%Generate the matrices from start to segment index for each segment
p.fibre.Hi_forward = single(NaN(2,2*p.fibre.nbSegments)); %Output array storing the 2*2 Jones matrix of each segment
%Hi_forward is the matrix from the start to the current segment : Hi_forward = U1*U2*...*Ui
p.fibre.Hi_forward(:,1:2) = U_i(:,1:2); %Initial matrix is the first segment matrix
for n=2:p.fibre.nbSegments
    p.fibre.Hi_forward(:,1+2*(n-1):2*n) = U_i(:,1+2*(n-1):2*n)*p.fibre.Hi_forward(:,1+2*(n-2):2*(n-1)); %Calculate matrix from start to segment index n
end

%Generate the round-trip matrices for each segment : Hi = Hi_forward.' * M_i * Hi_forward
p.fibre.Hi = single(NaN(2,2*p.fibre.nbSegments)); %Output array storing the 2*2 Jones matrix of each segment
if p.fibre.polar
    for n=1:p.fibre.nbSegments        
        p.fibre.Hi(:,1+2*(n-1):2*n) = p.fibre.Ei(n)*p.fibre.Ai(n)*p.fibre.Hi_forward(:,1+2*(n-1):2*n).'* p.fibre.Hi_forward(:,1+2*(n-1):2*n); %p.fibre.Ei(n)*p.fibre.Ai(n)*%Calculate round-trip matrix for the current fiber segment
    end
else
    for n=1:p.fibre.nbSegments
        p.fibre.Hi(:,1+2*(n-1):2*n) =p.fibre.Ai(n)*[1 0 ; 0 1 ]; % p.fibre.Ai(n) *p.fibre.Ei(n)* *%Theoretical without polar rotation & PDL
    end
end
Hi = p.fibre.Hi; %Output of the function

