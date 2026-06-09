function val = g_func(x,y)

tol = 1e-12;

if abs(y) < tol
    val = sin(x);

elseif abs(x) < tol
    val = sin(y);

else
    val = 0;
end

% val = 0; %test

end