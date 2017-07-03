function [S,dE,flip] = flip_sign(S,T)

global J N T_i dT T_f

r = randi([1,N]);

if r == 1
    dE = 2 * J * S(r) * (S(N) + S(r + 1));
end
if r == N
    dE = 2 * J * S(r) * (S(r - 1) + S(1));
end
if r ~= 1 && r ~= N
    dE = 2 * J * S(r) * (S(r - 1) + S(r + 1));
end 

if dE < 0
    S(r) = -1 * S(r);
    flip = 1; 
else
    r2 = rand(1);
    beta = 1/T;
    
    if r2 < exp(-dE * beta)
        S(r) = -1 * S(r);
        flip = 1; 
    else
        flip = 0;
    end
end