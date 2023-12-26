[stress, stiffness] = crystal_elasticity(x, epsilon, sigma, F);

sigmaV = stress * F';

% Display the stress and stiffness matrices
disp('Stress Matrix:');
disp(stress);

disp('Formation Matrix:');
disp(F);

disp('Sigma Tensor:');
disp(sigmaV);

% Assuming stiffness_tensor is a 3x3x3x3 stiffness tensor

% Get the size of the stiffness tensor
tensor_size = size(stiffness);

% Display the stiffness tensor using loops
disp('Stiffness Tensor:');
for i = 1:tensor_size(1)
    for j = 1:tensor_size(2)
        for k = 1:tensor_size(3)
            for l = 1:tensor_size(4)
                fprintf('C(%d, %d, %d, %d) = %f\n', i, j, k, l, stiffness(i, j, k, l));
            end
        end
    end
end
