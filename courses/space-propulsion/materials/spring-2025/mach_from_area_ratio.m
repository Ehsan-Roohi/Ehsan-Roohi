% Define constants
gamma = 1.4;  % Specific heat ratio for air
% Input the area ratio (A_e / A*)
A_Astar = input('Enter the value of the area ratio (A_e / A*): ');
% Define the area ratio equation as a function of Mach number
area_ratio_func = @(Me) (1/Me) * ((2/(gamma + 1)) * (1 + (gamma - 1)/2 * Me^2))^((gamma + 1)/(2*(gamma - 1))) - A_Astar;
% Define the bisection method parameters
tol = 1e-6;  % Tolerance for convergence
max_iter = 1000;  % Maximum number of iterations
lower_bound = 1;  % Start search from Mach number 1 (supersonic region)
upper_bound = 10;  % Upper bound for Mach number (this is a guess)
% Bisection method loop
for iter = 1:max_iter
    Me_mid = (lower_bound + upper_bound) / 2;  % Midpoint
    f_mid = area_ratio_func(Me_mid);  % Evaluate the function at midpoint
    % Check if solution is found within tolerance
    if abs(f_mid) < tol
        break;
    elseif area_ratio_func(lower_bound) * f_mid < 0
        upper_bound = Me_mid;  % The solution is in the lower half
    else
        lower_bound = Me_mid;  % The solution is in the upper half
    end
end