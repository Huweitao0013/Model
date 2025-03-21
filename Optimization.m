% Actual output values
actual_values = [
    13213.442, 13933.297, 18145.845, 17301.998,	18370.146, 22532.400, 23217.565, 23517.097, 19563.263, 19887.315; % I = 40
    7772.532, 7806.804, 7865.874, 8716.409, 9740.254, 10683.295, 9003.96, 8316.149, 7449.454, 9339.298; % I = 9
    4957.479, 4931.845, 5083.98, 5414.998, 5883.663, 5792.828, 5758.834, 6142.794, 6117.159, 5548.186; % I = 6
];
currents = [40, 9, 6]; % I

% Define the objective function for optimization
function error = multi_objective(params, actual_values, currents)
    S = params(1); % Optimization variable S
    DLi = params(2); % Optimization variable DLi
    Dcl = params(3); % Optimization variable Dcl
    t_max_values = params(4:6); % t_max corresponds to I=40, 9, 6
    errors = zeros(1, length(actual_values));

    % Calculate the error for each current condition separately
    for idx = 1:3 
        I = currents(idx);
        t_max = t_max_values(idx);
        CLi = 100 * ones(1, 10);
        Ccl = 100 * ones(1, 10);
        WESC = 0.8 * 23684 * ones(1, 10);
        u = 0 * ones(1, 10);
        k = 1e-11;
        time_step = 1;
    
        for t = 1:t_max
            S_eff = S * (0.5 .^ 1.5) * ones(1, 9);
            RS = 0.000025 ./ S_eff;
            K_eff = 96485.33 * 96485.33 / 8.314 / 298 * (CLi(1:9) * DLi * (0.5 .^ 1.5) + Ccl(1:9) * Dcl * (0.5 .^ 1.5));
            RL = 0.000025 ./ K_eff;
            j0 = 96485.33 * k * (WESC .^ 0.5) .* ((23684 - WESC) .^ 0.5) .* (CLi .^ 0.5);
            Rct = 8.314 * 298 ./ (10^6 * j0 * 96485.33 * 0.000025);
    
            % Create Rmatrix
            Rmatrix = zeros(9, 9);
            for n = 1:9
                for m = 1:9
                    if n == m
                        Rmatrix(n, m) = RS(n) + Rct(n) + Rct(n + 1) + RL(n);
                    elseif m == n + 1
                        if n + 1 <= 9
                            Rmatrix(n, m) = -Rct(n + 1);
                        end
                    elseif m == n - 1
                        if n - 1 >= 1
                            Rmatrix(n, m) = -Rct(n);
                        end
                    else
                        Rmatrix(n, m) = 0;
                    end
                end
            end
    
            % Create Vvector (9x1)
            Vvector = zeros(9, 1);
            for n = 1:9
                if n == 1
                    Vvector(n) = u(n) + I * Rct(n) + I * RL(n) - u(n + 1);
                else
                    Vvector(n) = u(n) + I * RL(n) - u(n + 1);
                end
            end
    
            % Calculate Ivector
            Ivector = Rmatrix \ Vvector;
    
            % Calculate is,n, il,n, Saj,m
            is = Ivector';
            il = I - Ivector';
            Saj = zeros(1, 10);
            Saj(1) = I - Ivector(1);
            for m = 2:9
                Saj(m) = Ivector(m-1) - Ivector(m);
            end
            Saj(10) = Ivector(9);
    
            % Create dWESC,m/dt gradCLi,m and gradCcl,m
            dWESC_dt = -Saj / (96485.33 * 0.000025 * 0.5);
            gradCLi = zeros(1, 10);
            for m = 1:9
                gradCLi(m) = (CLi(m+1) - CLi(m)) / 0.000025;
            end
            gradCLi(10) = (100 - CLi(10)) / 0.000025;
    
            gradCcl = zeros(1, 10);
            for m = 1:9
                gradCcl(m) = (Ccl(m+1) - Ccl(m)) / 0.000025;
            end
            gradCcl(10) = (100 - Ccl(10)) / 0.000025;
    
            % Calculate gradfL,m
            gradfL = zeros(1, 10);
            for m = 1:9
                gradfL(m) = -RL(m) .* il(m);
            end
            gradfL(10) = -RL(9) .* il(9);
    
            % Calculate dCLi,m/dt and dCcl,m/dt
            % Update dCLi_dt and dCcl_dt
            % Set dRLi_dt as JLi/h/0.5 and dRcl_dt as Jcl/h/0.5
            dCLi_dt = zeros(1, 10);
            dRLi_dt = zeros(1, 10);
            for n = 1:10
                for m = 1:10
                    dRLi_dt(m) = (-DLi * gradCLi(m) - 96485.33 / (8.314 * 298) * DLi * CLi(m) * gradfL(m)) / 0.000025 / 0.5;
                end
                if n == 1
                    dCLi_dt(n) = (-dWESC_dt(n) - dRLi_dt(n));
                else
                    dCLi_dt(n) = (-dWESC_dt(n) - dRLi_dt(n) + dRLi_dt(n-1));
                end
            end
    
            dCcl_dt = zeros(1, 10);
            dRcl_dt = zeros(1, 10);
            for n = 1:10
                for m = 1:10
                    dRcl_dt(m) = (-Dcl * gradCcl(m) - 96485.33 / (8.314 * 298) * Dcl * Ccl(m) * -gradfL(m)) / 0.000025 / 0.5;
                end
                if n == 1
                    dCcl_dt(n) = (- dRcl_dt(n));
                else
                    dCcl_dt(n) = (- dRcl_dt(n) + dRcl_dt(n-1));
                end
            end
    
            WESC = WESC + dWESC_dt .* time_step;
            CLi = CLi + dCLi_dt .* time_step;
            Ccl = Ccl + dCcl_dt .* time_step;
            u = 0.9 - 8.314 * 298 / 96485.33 * log(WESC ./ (18947.2 - WESC));
        end

         % Error calculation (with corrected indexing)
        actual = actual_values(idx, :); % Actual values
        predicted = WESC + CLi;    % Simulated values
        errors(idx) = max(abs(predicted - actual));
    end

    % Return the maximum value of all errors to ensure local minimization
    error = sum(errors);
end

% Initial guess values
initial_guess = [0.0182, 1.9538e-12, 7.6080e-11, 25000, 1000, 50000];

% Invoke the optimizer by passing parameters through an anonymous function
optimized_params = fminsearch(@(params) multi_objective(params, actual_values, currents), initial_guess);

% Output the optimization results
fprintf('S: %.4f\n', optimized_params(1));
fprintf('DLi: %.4e\n', optimized_params(2));
fprintf('Dcl: %.4e\n', optimized_params(3));
fprintf('t_max (I=40): %.0f\n', optimized_params(4));
fprintf('t_max (I=9): %.0f\n', optimized_params(5));
fprintf('t_max (I=6): %.0f\n', optimized_params(6));