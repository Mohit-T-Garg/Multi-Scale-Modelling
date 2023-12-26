% Initialize constants and parameters
clc
clear all
close all

global epsilon; global sigma; global N;
global rcut; global eps_ver; global MAX_NEIGHBOR;
global x_last; global Num; global Neighborlist;
global dim; global H; global H_mod

epsilon = 0.00088;  % eV
sigma = 2.56;       % angstrom
TOL = 1e-6;
rcut = 2^(1/6.0) * sigma * 4;
eps_ver = 0.2;
MAX_NEIGHBOR = 10;

% Input spatial dimension and number of basis atoms
dim = input('Enter the spatial dimension in which deformation is occurring: ');
N = input('Enter the number of basis atoms present in the unit cell/simulation box: ');

% Initialize lattice vector matrix H
H = eye(dim);

% for i = 1:dim
%     msg = strcat('Enter the components of lattice vector number ', int2str(i), ': ');
%     disp(msg);
%     for j = 1:dim
%         msg = strcat('Component ', int2str(j), ' of lattice vector number ', int2str(i), ': ');
%         H(j, i) = input(msg);
%     end
% end

%Initialize F mat
% rix
F = [sqrt(0.9) sqrt(0.1) 0; sqrt(0.01) sqrt(0.99) 0; 0 0 1];
H_mod = F * H;


% Initialize arrays and variables
n = 0;
x = zeros(dim * N, 1);
Force = zeros(dim * N, 1);

for i = 1:dim * N
    x(i) = rand();
end

% Set the first 'dim' elements of x to 0
for i = 1:dim
    x(i) = 0;
end

x_last = x - 2 * eps_ver * rcut * ones(dim * N, 1);

% Initialize energy, force, and stiffness matrix
[Energy, f, K] = potential(x);
iter = 0;
f_prev = f;
residual = norm(f);

% Main optimization loop
while (residual > TOL)
    for i = 1:N
        if dim == 3
            plot3(x(dim * (i - 1) + 1), x(dim * (i - 1) + 2), x(dim * (i - 1) + 3), '*');
        elseif dim == 2
            plot(x(dim * (i - 1) + 1), x(dim * (i - 1) + 2), '*');
        elseif dim == 1
            plot(x(dim * (i - 1) + 1), '*');
        end
        hold on
    end

    d = pinv(K) * f;
    x = x + d;

    for i = 1:N
        temp = x(dim * (i - 1) + 1:dim * i);
        temp1 = H_mod \ temp;
        temp_ceil = ceil(temp1);
        for j = 1:dim
            temp_ceil(j) = temp_ceil(j) - 1;
        end
        x(dim * (i - 1) + 1:dim * i) = temp - H_mod * temp_ceil;
    end

    f_prev = f;
    [Energy, f, K] = potential(x);
    residual = norm(f);
    iter = iter + 1;
    fprintf('Iteration: %d, Residual: %.6f, Energy: %.6f\n', iter, residual, Energy);
end

% Visualization of crystal deformation
figure
hold on

% for i = -5:5
%     for j = -5:5
%         for k = 1:N
%             temp = x(dim * (k - 1) + 1:dim * k);
%             temp = temp + H * [i; j];
%             if k == 1
%                 plot(temp(1), temp(2), '*b');
%             elseif k == 2
%                 plot(temp(1), temp(2), '*k');
%             elseif k == 3
%                 plot(temp(1), temp(2), '*g');
%             else
%                 plot(temp(1), temp(2), '*g');
%             end
%         end
%     end
% end

% Calculate stress and stiffness
[stress, stiffness] = crystal_elasticity(x, epsilon, sigma, F);

sigmaV = stress * F';

% Display the stress and stiffness matrices
disp('Stress Matrix:');
disp(stress);

disp('Stiffness Tensor:');
disp(stiffness);

disp('Sigma Tensor:');
disp(sigmaV);
