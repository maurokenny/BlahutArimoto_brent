function [C0 r] = BlahutArimoto_brent(p,x,P,rho)
if ~isempty(find(p < 0))
 disp('Error: some entry in the input matrix is negative')
 C0 = 0; return;
end

% Check that the input matrix p does not have zero column
column_sum = sum(p);
if ~isempty(find(column_sum == 0))
 disp('Error: there is a zero column in the input matrix');
 C0 = 0; return;
end
% Check that the input matrix p does not have zero row
row_sum = sum(p,2);
if ~isempty(find(row_sum == 0))
 disp('Error: there is a zero row in the input matrix');
 C0= 0; return;
else
 p = diag(sum(p,2))^(-1) * p; % Make sure that the row sums are 1
end
[m n] = size(p);

% (2) Initialize capacity C0, C1.
C0 = 1; C1 = 0;
% (1) Choose any r(x) such that 0 < r(x) < 1 and
% sum_x r(x) = 1.
r = ones(1,m)/m; % initial distribution for channel input
error_tolerance = 1e-4;
% (2) Initialize capacity C0, C1.
q = zeros(m,n);
while (abs(C0-C1) > error_tolerance)
    C1 = C0;
    for j = 1:n
        q(:,j) = r'.*p(:,j);
        q(:,j) = q(:,j)/sum(q(:,j));
    end
    %(3)
    C0 = 0;
    for i = 1:m
        for j = 1:n
            if r(i) > 0 && q(i,j) > 0
                C0 = C0 + r(i)*p(i,j)*log(q(i,j)/r(i));
            end
        end
    end
    C0 = C0/log(2); % Capacity in bits

    pro = zeros(m,1);
    for j = 1:m
        pro(j) = prod(q(j,:).^p(j,:));
    end
    temp = @(B1) sum( exp(B1.*abs(x).^rho).*(1-abs(x).^rho./P).*pro');
    B0 = fzero(temp,[-1 -0.00001])

    r = exp(B0.*abs(x).^rho).*pro'/sum(exp(B0.*abs(x).^rho).*pro');
    C0-C1
end