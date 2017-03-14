function rs = twofirst(phi,theta,P)

A = [1;2];
C = [2;1];
D = 1;

% phi = [1 2 3 4];
% theta = 0.5;
% P = 0.5;

% Component vector
% 1  2  3  4 
% A  B  C  D

S = [-1 1 3 0; -2 0 0 1]; % p x N
chi = abs(sign(S)); % 

Sb = [1; 2]; 

% only for first order reactions
ci = phi/sum(phi) * P * (1/(1 + D*theta));

calpha = [ci(1) ; ci(1)];

ra =  (A.*exp(C*(1-(1/(D.*theta+1)))) .* calpha ./ Sb);
ralpha = (chi .* S)' * ra;
halpha = [1 ; 5];

rs = struct('chi', chi,'Sb', Sb, 'ralpha', ralpha, 'halpha', halpha,'D',D)
    
    
