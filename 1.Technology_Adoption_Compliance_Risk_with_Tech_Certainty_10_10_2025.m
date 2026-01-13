% Date: Janury 12, 2026

% This set of MATLAB scripts calculates the investment thresholds, along with 
% optimal emissions and violation levels from the simulations reported in 
% Arguedas, C., Peinado, F., and Zofío, J.L. (2026) "Incentives for Green 
% Technology Adoption and Compliance under Risk Aversion and Technological Uncertainty",
% Environmental and Resource Economics, 89(6), art. nº 6.
% 
% Optimal emissions are calculated by solving the condition presented in 
% equation (2) of the article, minimizing firms' expected disutility in terms of 
% compliance and non-compliance (see equation (1)). 
% The scripts have been run in version R2024b of MATLAB and use the function
% 'vpasolve' from the Symbolic Math Toolbox to determine optimal emissions.
%
% This first script focuses on the simulations in subsection 4.2 of the paper 
% with technological certainty.
%
% Main outputs of the script: Degree of Relative Risk Aversion (Ar),
% Investment threshold (I), Optimal emissions with the old technology,
% Optimal emissions with the new technology, Violation level with the old technology
% Violation level with the new technology.

% Follow steps [1], [2], ... to ensure a correct replication of the results.

% Notation:
% pi=monitoring probability.
% i=fixed investment cost of installing the cleaner/new technology.
% e0=actual emissions with old technology.
% e1=actual emissions with new technology.
% e0c=optimal actual emissions with old technology.
% e1c=optimal actual emissions with new technology.
% r0=declared (reported) emissions with old technology.
% r1=declared (reported) emissions with new technologies.
% tau=tax on declared emissions.
% ff=fixed part of the fine.
% m=sanction multiplier.
% ef=exponent of the sanction function.
% te0=efficiency of the old technology (less efficient the larger this
% value).
% te1=efficiency of the new technology (less efficient the
% larger this value).
% gamma=gamma parameter of the hyperbolic disutility function.
% mu=mu parameter of the hyperbolic disutility function.
% nu=nu parameter of the hyperbolic disutility function.

clear all;
global pi rho i tau e0c e1c ff ef te0 te1 gamma mu nu m
syms e0 e1 r0 r1 

% [1] Assign values to the model parameters. Appendix A.2 in the article shows how the choice of 
% different values of gamma, mu and nu, results in specific formulations of the hyperbolic 
% disutility function: Power, Logarithmic (Bernoulli), Negative exponential, or Quadratic).

pi=0.5
tau=20
ff=0 %Use ff=39.9999999 as ff=40
m=1
ef=2
te0=10*10000
te1=5*10000
gamma=2
mu=1
nu=0

% [2] Defines the abatement costs functions for the old technology (0) and 
% the new technology (1), as well as their respective first derivatives (0d, 1d)

c0=te0/e0-1; % Maximum emissions with te0 are te0 & Maximum cost when e0=0 is Inf
c0d=-te0/(e0^2);
c1=te1/e1-1; % Maximum emissions with te0 are te1 & Maximum cost when e1=0 is Inf
c1d=-te1/(e1^2);

% Calculates the optimal emissions depending on the technology (0 or 1)

ope0eqn=c0d+tau==0;
ope0 = vpasolve(ope0eqn, e0, [0 te0]);
e0c=ope0;
ope1eqn=c1d+tau==0;
ope1 = vpasolve(ope1eqn, e1, [0 te1]);
e1c=ope1;

% [3] Assign abatement costs functions like in [2] substituting e0 and e1 
% by e0c and e1c, respectively. Assign sanctioning functions with e0c and e1c 
% keeping them equivalent to the case without technology uncertainty to allow comparability.
% Notice that f0d corresponds to the first derivative of f0 evaluated at 
% the violation level, that is, at e0c-r0 (which is not the same as de first derivative 
% of f0 evaluated at r0). It must be updated manually, and the same applies
% to f1d.

c0=te0/e0c-1;
c1=te1/e1c-1;

f0=ff*m*(e0c-r0)+m*(e0c-r0)^ef;
f0d=m*ff+m*ef*e0c-m*ef*r0;
f1=ff*m*(e1c-r1)+m*(e1c-r1)^ef;
f1d=m*ff+m*ef*e1c-m*ef*r1;

% [4] Backwards loop for several i values. Do not directly specify i values
% in i, instead, set counter=maximum value of i to be studied.

for counter=829:-1:0; %828 threshold when pi=0.5 tau=20 ff=0 m=1 ef=2 te0=10*10000 te1=5*10000 gamma=2 mu=1 nu=0

i=counter;

% Old Technology

% [5] Disutility function (see Appendix A.2). Notice that dN0d and dS0d correspond 
% to the first derivatives with respect to N0=(ci0+tau*ri0) and S0=(ci0+tau*ri0+fi0)
% of dN0 and dS0, respectively. They must be updated manually at this step.

dN0 = (gamma-1)*(mu*(c0+tau*r0)/(gamma-1)-nu)^gamma;
dN0d = gamma*mu*(mu*(c0+tau*r0)/(gamma-1)-nu)^(gamma-1);
dS0 = (gamma-1)*(mu*(c0+tau*r0+f0)/(gamma-1)-nu)^gamma;
dS0d = gamma*mu*(mu*(c0+tau*r0+f0)/(gamma-1)-nu)^(gamma-1);

% Proposition 1: Equation to implicitly obtain declared emissions

eqn = (dS0d*pi*f0d)/((1-pi)*dN0d+pi*dS0d) == tau;

% [6] The function 'vpasolve' solves the equation numerically. With the range
% [0 e0c] we consider only real solutions between 0 (which constitutes the
% lower bound because firms won't declare negative emissions) and the
% optimal emissions level with the old technology (which constitutes the upper
% bound because firms won't declare more emissions than actual emissions).
% Be careful to keep Fo the same as the sanctioning function in [3]
% substituting r0 by OT (the amount of declared emissions with the old technology)

OT = vpasolve(eqn, r0,[0 e0c]); % Amount of declared emissions with the old technology.

Vot=e0c-OT; % Violation level with the old technology

Cot=c0; % Abatement costs with the old technology

Tot=OT*tau; %Total taxes paid for declared emissions with the old technology

Fot=m*ff*(e0c-OT)+m*(e0c-OT)^ef; % Fine for the violation level with the old technology, Parametrizations 1 and 3: Fo=ff*(ei0c-O)+(ei0c-O)^2; Parametrization 2: Fo=ff*5*(ei0c-O)+5*(ei0c-O)^2;

% Shows total costs to get an idea about how large they are numerically (the independent variable)

CostsOT=Cot+Tot+Fot

% [7] Disutility assuming the old technology. Make sure that the disutility
% functions multiplying (1-pi) and pi are the same as the ones written
% in [5] but substituting c0, tau*r0 and f0 by Cot, Tot and Fot respectively.

Dot=(1-pi)*((gamma-1)*(mu*(Cot+Tot)/(gamma-1)-nu)^gamma)+pi*((gamma-1)*(mu*(Cot+Tot+Fot)/(gamma-1)-nu)^gamma);

% New Technology

% [8] Follow the same reasoning as in [5] when setting these disutility
% functions and keep them like the ones in [5] adding i to the total costs
% and being careful to work with 1h functions

dN1 = (gamma-1)*(mu*(c1+tau*r1+i)/(gamma-1)-nu)^gamma;
dN1d = gamma*mu*(mu*(c1+tau*r1+i)/(gamma-1)-nu)^(gamma-1);
dS1 = (gamma-1)*(mu*(c1+tau*r1+f1+i)/(gamma-1)-nu)^gamma;
dS1d = gamma*mu*(mu*(c1+tau*r1+f1+i)/(gamma-1)-nu)^(gamma-1);

% Proposition 1: Equation to implicitly obtain declared emissions

eqn = (dS1d*pi*f1d)/((1-pi)*dN1d+pi*dS1d) == tau;

% [9] The function 'vpasolve' solves the equation numerically. With the range
% [0 e1c] we consider only real solutions between 0 (which constitutes the
% lower bound because firms won't declare negative emissions) and the
% optimal emissions level with the new tech (which constitutes the upper
% bound because firms won't declare more emissions than actual
% emissions). Be careful to keep Fnt the same as the sanctioning function in [3]
% substituting r1 by NT (the amount of declared emissions with the new technology)

NT = vpasolve(eqn, r1,[0 e1c]); % Amount of declared emissions with the new technology

Vnt=e1c-NT; % Violation level with the least efficient new technology

Cnt=c1; % Abatement costs with the new technology

Tnt=NT*tau; % Total taxes paid for declared emissions with the new technology

Int=i; % Fixed investment cost in the new technology

Fnt=m*ff*(e1c-NT)+m*(e1c-NT)^ef; % Fine for the violation level with the new technology

% Shows total costs with the new technology to get an idea about how large
% they are numerically (the independent variable).

Costsn=Cnt+Tnt+Int+Fnt;

% [10] Disutility assuming the new technology. Make sure that the disutility
% functions multiplying (1-pi) and pi are the same as the ones written
% in [8] but substituting c1, tau*r1, i and f1 by Cnt, Tnt, Int and Fnt, respectively

Dnt=(1-pi)*((gamma-1)*(mu*(Cnt+Tnt+Int)/(gamma-1)-nu)^gamma)+pi*((gamma-1)*(mu*(Cnt+Tnt+Int+Fnt)/(gamma-1)-nu)^gamma)

% Now we check if the disutility with the old technology is larger than the
% disutility with the new technology. When this is the case, the firms invests in
% the new technology. This indicates the program that an approximate I threshold 
% has been found and therefore the loop ends.

if Dot-Dnt>0
    break
end
    
end

Ar=-mu*Costsn/(mu*Costsn/(gamma-1)-nu); % Calculates the degree of Relative Risk Aversion (Ar) using eq. (A2.6)
fprintf("Degree of Relative Risk Aversion (Ar) is = %s\n", Ar);
X=['Indifferent if the investment cost of the new technology (I) is aproximately between [',num2str(i),',',num2str(i+1),')'];
disp(X)
fprintf('Optimal emissions with the old technology = %s\n',ope0);
fprintf('Optimal emissions with the new technology = %s\n',ope1);
fprintf('Violation level with the old technology = %s\n',Vot);
fprintf('Violation level with the new technology = %s\n',Vnt);

