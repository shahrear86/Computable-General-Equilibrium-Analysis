% Md. Shahrear Zaman
% shahrear.zaman1971@gmail.com
% student.eco86@gmail.com
% 09:59 P.M.
% 11/08/2020 
% 11:36 A.M.
% 12/08/2020 % I have developed a Computable General Equilibrium model for a Power Generator 
% 12:44 P.M. 
% 14/01/2021 % I have tried to adjust statistical uncertainty by the calculated t-value with the Computable General Equilibrium Model
% 08:25 P.M. 
% 22/01/2021 

clear all
ncol = size(make,2);
Y = make(:,ncol);
k = columns(make);
[T  c] = size(Y) ;
N=T;
[count, value] = runlength (Y);

PP = [ ]
for i =1:T
      PP1 = (rows(Y(find(Y == Y(i))))/rows(Y))
      PP =[PP PP1]
end

mu  = sum(PP*Y);
s   = sum(Y)/N;
sigma2 = sum(((Y-mu).^2)'*PP');
sigma_hat2 = sum((Y- mean(Y)).^2)/(T - k);

w = 10  % weight on the random number
% making 1000 replication of the t value from the sample
% randn(1) generates a normally distributed random number, you can also check with 
% a gamma distributed random number by using randg(1) 
t = [ ];
 for itr = 1 : 1000
 tvelu =  (((sum(Y)+w*randn(1)) /rows(Y)) - mu) / sqrt(sigma2/T);
 t = [t  [tvelu]];
 end
 
tdim = prctile(t,[30,60,90]);
tdim
 
% Check t-values from the tdim, where three percentiles have taken from the replicated 1000 t-values  
% Adjust t-value with the coefficients in the equations to solve the system of the equations
% As fsolve has some limitations, I have adjusted t-value with the coefficients manually 
% by following a simplified form 
 
function y = f (p)
y = zeros (6, 1);
y(1) =  -0.37182*(0*p(1)+0*p(2)+0*p(3)+0.833333333*p(4)+0*p(5)+0*p(6)+0.166666667*p(7) - 0.083333333*p(1)-	0.166666667*p(2)-	0*p(3)-	0.083333333*p(4)-	0.833333333*p(5)-	0*p(6)-	0.166666667*p(7));
y(2) =  -0.37182*(0*p(1)+0*p(2)+1*p(3)+0*p(4)+0*p(5)+0*p(6)+0*p(7)-0*p(1)- 0*p(2)-1*p(3)-0*p(4)-0*p(5)-0*p(6)-0*p(7));
y(3) =  -0.37182*(0*p(1)+0*p(2)+0*p(3)+0.9375*p(4)+0*p(5)+0*p(6)+0.0625*p(7)-0*p(1)-0*p(2)-0*p(3)-0.6875*p(4)-	0*p(5)-	0*p(6)-	0.0625*p(7));
y(4) =  -0.37182*(0*p(1)+0*p(2)+0*p(3)+0*p(4)+0.941176471*p(5)+0*p(6)+0.058823529*p(7)-0*p(1)-0*p(2)- 0*p(3)- 0*p(4)-1.764705882*p(5)-0*p(6)-0.058823529*p(7));
y(5) =  -0.37182*(0.138888889*p(1)+0.138888889*p(2)+0*p(3)+0*p(4)+0*p(5)+0.694444444*p(6)+0.027777778*p(7)-0.208333333*p(1)-0.208333333*p(2)-0.416666667*p(3)-0.027777778*p(4)-0*p(5)-0.069444444*p(6)-0.069444444*p(7));
y(6) =  -0.37182*(0.095238095*p(1)+0.111111111*p(2)+0.476190476*p(3)+0*p(4)+0.26984127*p(5)+0*p(6)+0.047619048*p(7)-0*p(1)-0*p(2)-0*p(3)-0.174603175*p(4)-0*p(5)-0.714285714*p(6)-0*p(7));
endfunction

[x, fval, info] = fsolve (@f, [1; 2; 3; 4; 5; 6; 7])

p=x

Geninc =72 % Power Generator Box is total using and total making 72 unit value of the power
% You can check the indirect utility by increasing or decreasing the value of the Geninc , Geninc = 100 for example  
roh = 0.2  % degree of the elasticity of the substitutability
% You can check the indirect utility by increasing or decreasing the value of the roh , roh = 0.9 for example 
r = roh / (roh -1 ) % Elasticity of the Substitutions

n = rows(x)
Xn = []
for  i = 1:n
     Xn = (((p(i))^(r-1))*Geninc)/(sum((p).^(r)))
end

% Calculating the Indirect Utility from the utility function in the form of the Constant Elasticity of the Substitutions 
v = (((sum((p).^(r))) / ((sum((p).^(r))).^(roh))).^(1/roh))*Geninc
*************************************************************************************************************************
% *Additional Comment : In case of the misspecification of the model and for the other kinds of mistakes 
%                       do not hesitate to inform me.  

roh1 = [ 0.2 ,0.9 ,0.1 ,0.3 , 0.4 , 0.5, 0.8]'  % degree of the elasticity of the substitutability
% You can check the indirect utility by increasing or decreasing the value of the roh , roh = 0.9 for example 
r1 = roh1 ./ (roh1 -1 ) % Elasticities of the Substitutions

n = rows(x)
Xn = []
for  i = 1:n
     Xn1 = (((p(i)).^(r1-1)).*Geninc)./(sum((p).^(r1)))
end

% Calculating the Indirect Utility from the utility function in the form of
% the Constant Elasticities of the Substitutions 
v1 = (((sum((p).^(r1))) ./ ((sum((p).^(r1))).^(roh1))).^(1./roh1)).*Geninc

v2 = ((((p).^(r1)) ./ (((p).^(r1)).^(roh1))).^(1./roh1)).*Geninc


p0 = ones(7,1);

v02 = ((((p0).^(r1)) ./ (((p0).^(r1)).^(roh1))).^(1./roh1)).*Geninc


CV + Geninc  = (v02./((((p).^(r1)) ./ (((p).^(r1)).^(roh1))).^(1./roh1)));
CV = (v02./((((p).^(r1)) ./ (((p).^(r1)).^(roh1))).^(1./roh1))) - Geninc

% Comment: According to the result of this Power Generator Model, we need to find out the prime source of 
%          the extra power ( very high large unit) and need to improve the environment of the system, we need 
%          to reduce the system loss ( high unit) and need to build up very very strong power supply 
%          for the transit three and make the transit three strong enough to handle the power to carry the process.       




    
     