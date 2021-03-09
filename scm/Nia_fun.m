function nia = Nia_fun(state)
    
    offset = 0.5;
    nia = 600 * ones(size(state));
    affected = state > offset;
    nia(affected) = -300 * exp(-3*(state(affected) - offset)) + 600;
end