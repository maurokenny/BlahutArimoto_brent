clear p;

% alpha stable parameters (to be chosen)
% X,ALPHA,BETA,GAM,DELTA;
alpha = 2;
beta = 0;
gam = sqrt(0.5); %sqrt(0.5) to ~N(0,1)
delta = 0;
% creating p(y|x)
% range for x and y (input and output)
bound_inf = -20; %must be negative (to be chosen)
bound_sup = 20;  %must be positive (to be chosen)
step = 1;        % step in quantization (to be chosen)

% E[|X|^rho] < P
interv_P = 1:30;  %choose Power values
rho = 2; 
P = 1;   % power P
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
interv = bound_inf:step:bound_sup; % X = SaS input
x = interv;
% p(y|x) = SaS(y-x), in other words, a shift
% so we can use a toeplitz function to generate p(y|x)
l = interv+abs(bound_inf);
lin = stblpdf(l,alpha,beta,gam,delta);  %~ SalphaS
if beta == 0
    p = toeplitz(lin);
else
   c = -(interv+abs(bound_inf));
   col = stblpdf(c,alpha,beta,gam,delta);  %~ SalphaS   
   p = toeplitz(col,lin); 
end
    
% c = -(interv+abs(bound_inf));
% l = interv+abs(bound_inf);
% C is capacity
% r is the optimal input
% [C r] = BlahutArimotoConstraint(p,x,P,rho) % testing the Blahut algorithm

% Plot
optimal = zeros(numel(interv_P),1);
array = zeros(numel(interv_P),1);
for P=interv_P;
   array(P)   = BlahutArimoto_brent(p,x,P,rho);
   optimal(P) = 1/2*log2(1+P);
end
figure,plot(array),hold on, plot(optimal,'*')