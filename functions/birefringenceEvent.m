function HiEvent = birefringenceEvent(p, segIdx, birefringenceChange)

resolution = p.fibre.spatialRes/p.fibre.polCorrL; 
%function to generate the Jones Matrices for the codes at and after the birefringence event
p.fibre.evolBetaTab(segIdx) = p.fibre.evolBetaTab(segIdx)+birefringenceChange; % add a sharp change in beta at segment segIdx
% Generate Jones matrices of each segment
U_i = single(NaN(2,2*p.fibre.nbSegments)); %Output array storing the 2*2 Jones matrix of each segment
for n=1:p.fibre.nbSegments 
    expjBeta = exp(1j*(p.fibre.evolBetaTab(n))); 
    sinRot = sin(0); cosRot = cos(0);%sinRot = sin(p.fibre.evolThetaTab(n)); cosRot = cos(p.fibre.evolThetaTab(n));%sinRot = sin(0); cosRot = cos(0);%sinRot = sin(p.fibre.evolThetaTab(n)); cosRot = cos(p.fibre.evolThetaTab(n)); %sinRot = sin(0); cosRot = cos(0);%
    U_i(:,1+2*(n-1):2*n) = [cosRot sinRot ; -sinRot cosRot] * [expjBeta 0 ; 0 conj(expjBeta)]*[cosRot -sinRot ; sinRot cosRot];%Calculate matrix U_forw for the current fiber segment
end

%Generate the matrices from start to segment index for each segment
Hi_forward = single(NaN(2,2*p.fibre.nbSegments)); %Output array storing the 2*2 Jones matrix of each segment
%Hi_forward is the matrix from the start to the current segment : Hi_forward = Ui*Ui-1*...*U1
Hi_forward(:,1:2) = U_i(:,1:2); %Initial matrix is the first segment matrix
for n=2:p.fibre.nbSegments
    Hi_forward(:,1+2*(n-1):2*n) = U_i(:,1+2*(n-1):2*n)*Hi_forward(:,1+2*(n-2):2*(n-1)); %Calculate matrix from start to segment index n
end

%Generate the round-trip matrices for each segment : Hi = Hi_forward.' * M_i * Hi_forward
HiEvent = single(NaN(2,2*p.fibre.nbSegments)); %Output array storing the 2*2 Jones matrix of each segment
if p.fibre.polar
    for n=1:p.fibre.nbSegments        
        HiEvent(:,1+2*(n-1):2*n) =p.fibre.Ai(n)*Hi_forward(:,1+2*(n-1):2*n).'* Hi_forward(:,1+2*(n-1):2*n); %Calculate round-trip matrix for the current fiber segment
    end
else
    for n=1:p.fibre.nbSegments
        HiEvent(:,1+2*(n-1):2*n) =p.fibre.Ai(n)*[1 0 ; 0 1 ]; % p.fibre.Ai(n) *p.fibre.Ei(n)* *%Theoretical without polar rotation & PDL
    end
end


