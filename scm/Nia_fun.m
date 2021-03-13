function nia = Nia_fun(state)
    
    offset = 3;
    nia = 600 * ones(size(state));
    affected = state > offset;
    nia(affected) = -500 * exp(-1*(state(affected) - offset)) + 600;
    nia = max(nia, 300);
end