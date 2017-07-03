function plot_avg(E_avg,M_sq_avg)

global J N T_i dT T_f
syms t

Tpts = [T_i:dT:T_f]; 
figure
plot(Tpts,E_avg); 
hold all

if N == 2
    E_exp = -N * J * tanh(2 * 1/t * J); 
else
    E_exp = -N * J * (tanh(1/t * J) + tanh(1/t * J)^(N - 1)) / (1 + tanh(1/t * J)^N); 
end

plot(Tpts,subs(E_exp,Tpts));
axis([0.5 10 -5 0]);
xlabel('T'); 
ylabel('E_avg'); 
title('E_avg vs. T'); 

figure
plot(Tpts,M_sq_avg);
hold all

% beta = 1./Tpts;
% Z = 2 * exp(-beta .* -3) + 6 * exp(-beta .* 1); 
% M_sq_exp = 2 * exp(-beta .* -3) ./ Z * 9 + 6 * exp(-beta .* 1) ./ Z * 1; 

M_sq_exp = N * (1 + tanh(1./[T_i:dT:T_f] * J)) ./ (1 - tanh(1./[T_i:dT:T_f] * J)); 

plot(Tpts,M_sq_exp); 
axis([0.5 10 0 1e3]);
xlabel('T'); 
ylabel('M^2'); 
title('M^2 vs. T'); 