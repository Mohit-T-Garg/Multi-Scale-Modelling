function [stress, stiffness] = crystal_elasticity(x, epsilon,sigma, F)
  global epsilon; global sigma; global N;global rcut;
  global Num; global Neighborlist; global dim; global H; global H_mod
  neighborlist_generate(x);
  stress = zeros(dim, dim);
  stiffness = zeros(dim, dim, dim, dim);
  Vo = det(H); % unit cell volume

  A = zeros(dim * N, dim, dim);
  for alpha = 1:N
    for beta = 1:Num(alpha)
      gamma = Neighborlist(alpha,beta,1);
      n = reshape(Neighborlist(alpha,beta,2:end),[],1); % relative index between alpha and gamma atoms
      r_vec = x(dim*(alpha-1)+1:dim*alpha)-x(dim*(gamma-1)+1:dim*gamma);
      r_vec = r_vec - H_mod * n;
      r = norm(r_vec);

      nu_d = 24*epsilon/r*(-2*(sigma/r)^12+(sigma/r)^6); % first derivative of nu
      nu_dd = 24*epsilon/r^2 * (26*(sigma/r)^12 - 7*(sigma/r)^6); % second derivative of nu
      
  
      stress = stress - nu_d * (r_vec/r) * (H * n)';

      tsrpdt1 = (r_vec/r) * (H * n)';
      tsrpdt2 = reshape(tsrpdt1,1,1,dim,dim);
      term1 = (nu_dd - nu_d/r) * tsrpdt1 .* tsrpdt2;
      
      t2 = reshape(r_vec/r, 1, 1, dim);
      tem1 = (-nu_dd + nu_d/r) * tsrpdt1 .* t2;

      delta_ik = reshape(eye(dim), dim,1,dim,1);
      Hn_j = reshape(H * n, 1,dim,1,1);
      Hn_l = reshape(H * n, 1,1,1,dim);
      term2 = nu_d/r * delta_ik .* Hn_j .* Hn_l;

      tem2 = - nu_d/r * delta_ik .* Hn_j ;

      newterm = tem1 + tem2;

      A(dim * (alpha - 1) + 1:dim * alpha, :, :) = A(dim * (alpha - 1) + 1:dim * alpha, :, :) + 0.5 * newterm;
      A(dim * (gamma - 1) + 1:dim * gamma, :, :) = A(dim * (gamma - 1) + 1:dim * gamma, :, :) - 0.5 * newterm;
      
      stiffness = stiffness + term1 + term2;
    end
  end  
  stress = stress / (2*Vo);
  stiffness = stiffness / (2*Vo);

[Energy, f, K] = potential(x);

for i = 1:dim
    W(i) = i;
end
for i=1: length(W)
    for j=1:dim*N
        if W(i)==j
            A(W(i),j)=1;
        else
            A(W(i),j)=0;
            A(j, W(i))=0;
        end
    end
end

A_reshaped = reshape(A, dim * N, []);

% Compute the product A * K
AK_product = (A_reshaped' * pinv(K)) * A_reshaped;

result_matrix = reshape(AK_product, dim, dim, dim, dim);
stiffness = stiffness - result_matrix;

end
