% Initialize parameters
S = 0.0309;
CLi = 100 * ones(1, 10);
Ccl = 100 * ones(1, 10);
WESC = 0.8 * 23684 * ones(1, 10);
t_max = 3604;
DLi = 3.9713e-12;
Dcl = 2.9996e-10;
time_step = 1;
u = 0 * ones(1, 10);
k = 1e-11;
I = 40;

% Introduce logic mask to record fixed WESC
fixed_WESC = false(1, 10);

for t = 1:1:t_max
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
    
    % Calculate dWESC,m/dt gradCLi,m and gradCcl,m
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
            dRLi_dt(m) = (-DLi * (0.5 .^ 1.5) * gradCLi(m) - 96485.33 / (8.314 * 298) * DLi * (0.5 .^ 1.5) * CLi(m) * gradfL(m)) / 0.000025 / 0.5;
        end
        if fixed_WESC(n)
            dWESC_dt(n) = 0; % Focus dWESC_dt as 0
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
            dRcl_dt(m) = (-Dcl * (0.5 .^ 1.5) * gradCcl(m) - 96485.33 / (8.314 * 298) * Dcl * (0.5 .^ 1.5) * Ccl(m) * -gradfL(m)) / 0.000025 / 0.5;
        end
        if n == 1
            dCcl_dt(n) = (-dRcl_dt(n));
        else
            dCcl_dt(n) = (-dRcl_dt(n) + dRcl_dt(n-1));
        end
    end

    % Update WESC
    WESC = WESC + dWESC_dt .* time_step;

    % If WESC is less than or equal to 0, set it to 1e-8 and fix dWESC_dt to 0
    fixed_indices = (WESC <= 0);
    WESC(fixed_indices) = 1e-8;
    dWESC_dt(fixed_indices) = 0;
    fixed_WESC(fixed_indices) = true;

    % For the fixed WESC and dWESC_dt, retain their values unchanged
    WESC(fixed_WESC) = 1e-8;
    dWESC_dt(fixed_WESC) = 0;

    % Update CLi 和 Ccl
    CLi = CLi + dCLi_dt .* time_step;
    Ccl = Ccl + dCcl_dt .* time_step;

    % Update u
    u = 0.9 - 8.314 * 298 / 96485.33 * log(WESC ./ (18947.2 - WESC));
end

% Output results
disp('At t = max:');
disp('WESC,m:');
disp(WESC);
disp('CLi,m:');
disp(CLi);
disp('WESC,m + CLi,m:');
disp(WESC + CLi);
  
