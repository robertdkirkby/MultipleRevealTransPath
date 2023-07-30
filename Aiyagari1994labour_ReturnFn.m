function F=Aiyagari1994labour_ReturnFn(h,aprime, a, z,r,w,mu_c,mu_l,psi,tau)
% The return function is essentially the combination of the utility
% function and the constraints.
% Is a minor extension of Aiyagari (1994) which endogenizes labor supply
% (h) and adds a labor income tax (tau).

F=-Inf;
c=(1-tau)*w*h*z+(1+r)*a-aprime; % Budget Constraint

if c>0
    if mu_c==1
        F=log(c)+log(1-h);
    else
        F=(c^(1-mu_c))/(1-mu_c) + psi*((1-h)^(1-mu_l))/(1-mu_l);
    end
end

end