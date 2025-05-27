function p = rxSensingCorrelation(p)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Subblock partitioning for correlation and Jones matrix extractions 
% Subblock partitioning is applied if length of captured or simulated data
% is longer than a given block length limit 'brutBlockLen'
%
% Authors: 
% Original code by S. Guerrier, C. Dorize & E. Awwad - 2022
% Modified version by E. Awwad - 2024 elie.awwad@telecom-paris.fr
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global Erx;

tic;
%Correlation process per multiple small blocks. The block size must be specified (p.corrPerBlockTwin, in seconds, with p.corrPerBlockTwin>>length(gCode))
brutBlockLen = pow2(floor(log2(p.rx.corrPerBlockTwin*p.tx.fSymb+length(p.rx.gCode)-1)));
%Test if blocklen is accurately computed (case short t_probing)
if p.rx.ErxLen < brutBlockLen % If acquired or simulated backscaterring is too short (as if T_probing << 1s)
    fprintf('\n\n blem: ErxLen is %d shorter than default blocklen %d. Note that gCode is %d \n', p.rx.ErxLen, brutBlockLen, length(p.rx.gCode));
    brutBlockLen = floor(p.tx.nbCodes*length(p.rx.gCode)/p.tx.subcarriers); 
    fprintf('\nBlock length set to %d = nbCodes*%d/subcarriers \n', brutBlockLen, length(p.rx.gCode));
elseif brutBlockLen<5*length(p.rx.gCode) %Check brutBlockLen>>hLen % Not tested yet - EA 27/05/25
    fprintf('\nBlock length %d is too short for correlation per block procedure. Must be greater than hLen=%d. Quit\n', brutBlockLen, length(Erx));
%     brutBlockLen = (p.tx.nbCodes-1)*length(gCode)-1;
    brutBlockLen = 8*length(p.rx.gCode); 
%     fprintf('\nBlock length set to %d = nbCodes*%d \n', brutBlockLen, length(gCode));
%     %return;
end

p = fcorrAndGetJones(p.rx.gCode, brutBlockLen,p); %Processing per block
fprintf('\nCode correlation AND Jones Matrix extraction processes completed within %.1f seconds. %d blocks of %d samples (>%d samples per code)\n', toc, floor(length(Erx)/(brutBlockLen-length(Erx)+1)), brutBlockLen-length(Erx)+1, length(Erx));

