function gf = gen_func(b, dt, death_prob)
% Calculates the generating function
    if b == 0
        gf = 1;
    else
        es = exp(b*dt*(1:length(death_prob)));
        gf = death_prob*es.';
    end
end
