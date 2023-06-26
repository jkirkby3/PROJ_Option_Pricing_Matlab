function params = getModel( model, case_id )
% Helper for getting the model used in testing
params = {};

if nargin < 2
    case_id = 1;
end

if model == 1 %BSM (Black Scholes Merton)
    if case_id == 1
        params.sigmaBSM = 0.45;    %CHOOSE
    elseif case_id == 2
        params.sigmaBSM = 0.1;    %CHOOSE
    elseif case_id == 3
        params.sigmaBSM = 0.2;    %CHOOSE    
    elseif case_id == 4
        params.sigmaBSM = 0.3;    %CHOOSE  
    elseif case_id == 5
        params.sigmaBSM = 0.4;    %CHOOSE  
    elseif case_id == 6
        params.sigmaBSM = 0.5;    %CHOOSE 
    end
    params.sigma = params.sigmaBSM;
    params.jumpModel = 0;
    
elseif model == 2 %CGMY
    params.C  = 0.03; 
    params.G  = 1; 
    params.MM = 10; 
    params.Y  = 1.3;
    
elseif model == 3 %NIG
    params.alpha = 6.1882;
    params.beta  = -3.8941;
    params.delta = 0.1622;
    
elseif model == 4 %MJD (Merton Jump Diffusion)
    params.sigma  = 0.12;
    params.lam    = 0.4;
    params.muj    = -0.12;
    params.sigmaj = 0.18;
    
    params.jumpModel = 1;
    
elseif model == 5 %Kou Double Expo
    params.sigma = 0.15;
    params.lam   = 3;
    params.p_up  = 0.2;
    params.eta1  = 25;
    params.eta2  = 10;
    
    params.jumpModel = 2;
    
    if case_id == 1
        % params.lam   = 3;
    elseif case_id == 2
        params.lam   = 1.5;
    elseif case_id == 3
        params.lam   = 2;
    elseif case_id == 4
        params.lam   = 2.5;
    elseif case_id == 5
        params.lam   = 3;
    elseif case_id == 6
        params.lam   = 3.5;
    end

elseif model == 8 % Variance Gamma 
    params.sigma = 0.2; 
    params.nu = 0.85;  
    params.theta = 0;    
    
elseif model == 9 % Bilateral Gamma 
    if case_id == 1
        params.alpha_p = 1.18; 
        params.lam_p = 10.57;  
        params.alpha_m = 1.44; 
        params.lam_m = 5.57;
    elseif case_id == 2
        params.alpha_p = 1.18; 
        params.lam_p = 4.50;  
        params.alpha_m = 1.44; 
        params.lam_m = 3.75;
    elseif case_id == 3
        params.alpha_p = 1.18; 
        params.lam_p = 10.57;  
        params.alpha_m = 1.18; 
        params.lam_m = 5.57;
    elseif case_id == 4
        params.alpha_p = 1.18; 
        params.lam_p = 7.68;  
        params.alpha_m = 1.18; 
        params.lam_m = 7.68;
    end
end

end

