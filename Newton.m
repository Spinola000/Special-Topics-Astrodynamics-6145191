function solution = Newton(func, x0, tol, maxit)
% Function to apply the Newton method to find the zeros of a function.
% func - Function which the zeros are to be found.
% x0 - Initial guess.
% tol - tolerance applied to the functio i.e how close will it get to zero.
% maxit - Maximum number of iterations allowed.
%
    syms x
    f_sym = func(x);               
    f_prime = diff(f_sym);         
    f2 = matlabFunction(f_prime);  
    f1 = matlabFunction(f_sym);    
    
    x_curr = x0;
    
    for i = 1:maxit
        x_next = x_curr - f1(x_curr) / f2(x_curr);
        
        if abs(x_next - x_curr) < tol
            break;
        end
        
        x_curr = x_next;
    end
    
    solution = x_curr;
end