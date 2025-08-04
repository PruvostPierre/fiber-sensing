function [G1,G2] = GenGolay(basis,order,G1,G2)

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Generic function to generate Golay complementary sequences
% basis: the basis of the Golay seed: 0, 1, 10, 26 or 20 (later case: complex)
% Seed 0 means that we use sequences G1 & G2 entered in the fct arguments.
% order: the number of recursions to come up with the desired length
% (final length per sequence: basis x 2^order)
%
% Authors: 
% Original code by C. Dorize - 2018
% Modified version by E. Awwad - 2024 elie.awwad@telecom-paris.fr
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% Other seed and recursion leading to complementary sequences, to be tested
% Seed G=([1;-1]); %G1;G2
% Recursion G = ([G(1,:) G(2,:);-G(1,:) G(2,:)]);

if basis == 1 || nargin<2
    G1 = (1);
    G2 = (1);
elseif basis == 10
    G1 = ([-1 1 1 -1 1 -1 1 1 1 -1]);
    G2 = ([-1 1 1 1 1 1 1 -1 -1 1]);
elseif basis == 26
    G1 = ([1 -1 1 1 -1 -1 1 -1 -1 -1 -1 1 -1 1 -1 -1 -1 -1 1 1 -1 -1 -1 1 -1 1]);
    G2 = ([-1 1 -1 -1 1 1 -1 1 1 1 1 -1 -1 -1 -1 -1 -1 -1 1 1 -1 -1 -1 1 -1 1]);
elseif basis == 20
    G1 = [1 1i 1i 1 -1 1 1i -1 -1 1 1 -1 -1 1i 1 -1 1 1i 1i 1];
    G2 = [1 1i 1i 1 -1 1 -1i 1 1 1 -1 -1 -1 1i -1 1 -1 -1i -1i -1];
elseif basis == 0 && nargin == 4
    %G1 and G2 read from the function arguments
else
    sprintf('Wrong choice or number of arguments in fct GenGolayGeneric(basis,order,G1,G2). Return');
    return
end    
    
%order = order-log2(length(G1));
for n=1 : order
    tmp1 = G1;
    tmp2 = G2;
    G1 = ([tmp1 tmp2]);
    G2 = ([tmp1 -tmp2]);
end
