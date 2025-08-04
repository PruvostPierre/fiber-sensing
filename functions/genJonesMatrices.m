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
    rng(p.seeds.seed_rotPol); % polarization correlation length not exceeded, Wiener process
    evolThetaTab = [single(asin(sqrt(rand()))); single(asin(sqrt(rand(p.fibre.nbSegments,1)))*resolution)];
    rng(p.seeds.seed_beta);
    evolBetaTab = [single(rand()*2*pi-pi); single((rand(p.fibre.nbSegments-1,1)*2*pi-pi)*resolution)];
    rng(p.seeds.seed_gamma);
    evolGammaTab = [single(rand()*2*pi-pi); single((rand(p.fibre.nbSegments-1,1)*2*pi-pi)*resolution)];
    p.fibre.rotPolarTab = cumsum(evolThetaTab);
    p.fibre.betaTab = unwrap(cumsum(evolBetaTab)); 
    p.fibre.gammaTab = unwrap(cumsum(evolGammaTab));
    fprintf('\n *** Lsegment < polarization beat length : smooth evolution with fiber segment to polarization beat length ratio r = %d *** \n', resolution); 
end

% Polar rotation at fiber entrance
rng(p.seeds.seed_theta); 
theta = 2*pi*rand()-pi;
p.fibre.R0 = [cos(theta) sin(theta) ; sin(theta) -cos(theta)];

% Create matrices
M_i = [sqrt(1-p.fibre.AlphaPolRay), -sqrt(p.fibre.AlphaPolRay); sqrt(p.fibre.AlphaPolRay), sqrt(1-p.fibre.AlphaPolRay)]; %mirror effect on polar in each fiber segment i. M is identity if alpha null

Hi = single(NaN(2,2*p.fibre.nbSegments)); %Output array storing the 2*2 Jones matrix of each segment
if p.fibre.polar
    for n=1:p.fibre.nbSegments        
        expjGamma = exp(1j*p.fibre.gammaTab(n)); expjBeta = exp(1j*p.fibre.betaTab(n)); 
        sinRot = sin(p.fibre.rotPolarTab(n)); cosRot = cos(p.fibre.rotPolarTab(n)); 
        
        Uf = [expjGamma 0 ; 0 conj(expjGamma)] * [cosRot sinRot ; -sinRot cosRot] * [expjBeta 0 ; 0 conj(expjBeta)];%Calculate matrix U_forw for the current fiber segment
        Hi(:,1+2*(n-1):2*n) =  p.fibre.Ai(n)*p.fibre.Ei(n) * Uf.' * M_i * Uf * p.fibre.R0;%Ei is reflected amplitude and phase of segment n
    end
else
    for n=1:p.fibre.nbSegments
        Hi(:,1+2*(n-1):2*n) =  p.fibre.Ai(n)*p.fibre.Ei(n) *[1 0 ; 0 1 ] ;%Theoretical without polar rotation & PDL
    end
end
