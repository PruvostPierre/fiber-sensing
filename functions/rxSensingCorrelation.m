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
brutBlockLen = pow2(floor(log2(p.rx.corrPerBlockTwin*p.tx.fSymb+length(p.tx.gCode)-1)));
%Test if blocklen is accurately computed (case short t_probing)
if p.rx.ErxLen < brutBlockLen %MODEL case (mostly) : if fiber+nbcodes too small (as if T_probing << 1s)
    fprintf('\n\n blem: ErxLen is %d shorter than default blocklen %d. Note that gCode is %d \n', p.rx.ErxLen, brutBlockLen, length(p.tx.gCode));
    brutBlockLen = floor(p.tx.nbCodes*length(p.tx.gCode)/p.tx.subcarriers); 
    fprintf('\nBlock length set to %d = nbCodes*%d/subcarriers \n', brutBlockLen, length(p.tx.gCode));
elseif brutBlockLen<5*length(p.tx.gCode) %Check brutBlockLen>>hLen %FIXME
    fprintf('\nBlock length %d is too short for correlation per block procedure. Must be greater than hLen=%d. Quit\n', brutBlockLen, length(Erx));
%     brutBlockLen = (p.tx.nbCodes-1)*length(gCode)-1;
    brutBlockLen = 8*length(p.tx.gCode); 
%     fprintf('\nBlock length set to %d = nbCodes*%d \n', brutBlockLen, length(gCode));
%     %return;
end

p = fcorrAndGetJones(p.tx.gCode, brutBlockLen,p); %Processing per block
fprintf('\nCode correlation AND Jones Matrix extraction processes completed within %.1f seconds. %d blocks of %d samples (>%d samples per code)\n', toc, floor(length(Erx)/(brutBlockLen-length(Erx)+1)), brutBlockLen-length(Erx)+1, length(Erx));

