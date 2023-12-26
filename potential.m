function [Energy, Force, K] = potential(x)
global epsilon; global sigma; global N;global rcut;
global Num; global Neighborlist; global dim; global H; global H_mod;


Energy = 0;
Force = zeros(dim * N, 1);
K = zeros(dim * N, dim * N);
disp("calling potential")

neighborlist_generate(x);

for alpha = 1:N
    for beta = 1:Num(alpha)
        gamma = Neighborlist(alpha,beta,1); %simply index of neighboring atom 
        n = zeros(dim,1);
        for i=1:dim
            n(i) = Neighborlist(alpha,beta,1+i); %translations needs to reach gamma
        end
        r_vec = x(dim*(alpha-1)+1:dim*alpha)-x(dim*(gamma-1)+1:dim*gamma); %finds normaal distance between 
        r_vec = r_vec - H_mod * n; %correction required so as to take into consideration the periodic image position of the neighboring atoms
        r = norm(r_vec);

        if (r < rcut)
            nu = 4 * epsilon * ((sigma / r)^12 - (sigma / r)^6);
            Energy = Energy + nu / 2;
            nu_d = 24 * epsilon / r * (-2 * (sigma / r)^12 + (sigma / r)^6);
            nu_dd = 24*(epsilon/r^2) * ( 26 *(sigma/r)^12 - 7*(sigma/r)^6);
            Force(dim * (alpha - 1) + 1:dim * alpha) = Force(dim * (alpha - 1) + 1:dim * alpha) - 0.5 * nu_d * r_vec / r;
            Force(dim * (gamma - 1) + 1:dim * gamma) = Force(dim * (gamma - 1) + 1:dim * gamma) + 0.5 * nu_d * r_vec / r;
            A=0.5 * nu_d * (eye(dim) - r_vec * r_vec' / r^2) +(0.5*nu_dd*(r_vec * r_vec' / r^2));
            K(dim * (alpha - 1) + 1:dim * alpha, dim * (alpha - 1) + 1:dim * alpha) = K(dim * (alpha - 1) + 1:dim * alpha, dim * (alpha - 1) + 1:dim * alpha) + A;
            K(dim * (gamma - 1) + 1:dim * gamma, dim * (gamma - 1) + 1:dim * gamma) = K(dim * (gamma - 1) + 1:dim * gamma, dim * (gamma - 1) + 1:dim * gamma) + A;
            K(dim * (alpha - 1) + 1:dim * alpha, dim * (gamma - 1) + 1:dim * gamma) = K(dim * (alpha - 1) + 1:dim * alpha, dim * (gamma - 1) + 1:dim * gamma) - A;
            K(dim * (gamma - 1) + 1:dim * gamma, dim * (alpha - 1) + 1:dim * alpha) = K(dim * (gamma - 1) + 1:dim * gamma, dim * (alpha - 1) + 1:dim * alpha) - A';

        end
    end
end
for i = 1:dim
    W(i) = i;
end
for i=1: length(W)
    for j=1:dim*N
        if W(i)==j
            K(W(i),j)=1;
            Force(W(i))=0;
        else
            K(W(i),j)=0;
            K(j, W(i))=0;
        end
    end
end
end