tic; 

global J N MAX T_i dT T_f

E_avg = [];
M_sq_avg = [];

J = 1;
N = 3;
MAX = 1e5; 
T_i = 0.5;
dT = 1e-0;
T_f = 10; 

T = T_i;

while T <= T_f
    
    S = ones(1,N);
    
    % Initialize E and dE. 
    E = [];
    E_sum = 0;
    for n = 1:(N - 1)
        E_sum = E_sum + S(n) * S(n + 1);
    end
    E_sum = E_sum + S(N) * S(1);
    E_f = -J * E_sum; 
    
    M_sq = [];
    
    flip = 0; 
    count = 0;
    
    for i = 1:MAX % Run 100,000 times.
        
        [S,dE,flip] = flip_sign(S,T);
        if flip == 1
            E_f = E_f + dE; 
        end

        if count == (10 * N) % 10*N steps.
            count = 0;
            
            E = [E; E_f];
            M_sq = [M_sq; sum(S)^2];
        end
        count = count + 1;
    end
    
    E_avg = [E_avg; mean(E)];
    M_sq_avg = [M_sq_avg; mean(M_sq)];
    
    E = []; M_sq = [];
    
    T = T + dT;
end

T = toc; 
disp(T); 

% plot_avg(E_avg,M_sq_avg);

% Data extrapolation. 
data(:,1) = [T_i:dT:T_f]; 
data(:,2) = E_avg;

if N == 2
    data(:,3) = -N * J * tanh(2 * 1./[T_i:dT:T_f] * J); 
else
    data(:,3) = -N * J * (tanh(1./[T_i:dT:T_f] * J) + tanh(1./[T_i:dT:T_f] * J).^(N - 1)) ./ (1 + tanh(1./[T_i:dT:T_f] * J).^N); 
end

data(:,4) = M_sq_avg; 
data(:,5) = N * (1 + tanh(1./[T_i:dT:T_f] * J)) ./ (1 - tanh(1./[T_i:dT:T_f] * J)); 
disp(data);  