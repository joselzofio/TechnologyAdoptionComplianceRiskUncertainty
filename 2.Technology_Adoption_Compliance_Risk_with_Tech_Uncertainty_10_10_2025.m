% Date: January 12, 2025

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
% This second script focuses on the simulations in subsection 4.2 of the paper 
% with technological uncertainty.
%
% Main outputs of the script: Degree of Relative Risk Aversion (Ar),
% Investment threshold (I), Optimal emissions with the old technology,
% Expected optimal emissions with the new technology, 
% Optimal emissions with the least efficient new technology,
% Optimal emissions with the most efficient new technology,
% Violation level with the old technology,
% Expected violation level with the new technology,
% Violation level with the least efficient new technology,
% Violation level with the most efficient new technology.

% Follow steps [1], [2], ... to ensure a correct replication of the results.

% Notation:
% pi=monitoring probability.
% i=fixed investment cost of installing the cleaner/new technology.
% e0=actual emissions with old technology.
% e1h=actual emissions with the least efficient new technology.
% e1l=actual emissions with the most efficient new technology.
% e0c=optimal actual emissions with old technology.
% e1hc=optimal actual emissions with the least efficient new technology.
% e1lc=optimal actual emissions with the most efficient new technology.
% r0=declared (reported) emissions with old technology.
% r1=declared (reported) emissions with new technologies.
% tau=tax on declared emissions.
% ff=fixed part of the fine.
% m=sanction multiplier.
% ef=exponent of the sanction function.
% te0=efficiency of the old technology (less efficient the larger this
% value).
% te1h=efficiency of the least efficient new technology (less efficient the
% larger this value).
% te1l=efficiency of the most efficient new technology (less efficient the
% larger this value).
% alpha=likelihood of least efficient abatement new technology.
% gamma=gamma parameter of the hyperbolic disutility function.
% mu=mu parameter of the hyperbolic disutility function.
% nu=nu parameter of the hyperbolic disutility function.


clear all;
global pi i tau e0c e1hc e1lc alpha ff m ef te0 te1h te1l gamma mu nu
syms e0 e1h e1l r0 r1


% [1] Assign values to the model parameters. Appendix A.2 in the article shows how the choice of 
% different values of gamma, mu and nu, results in specific formulations of the hyperbolic 
% disutility function: Power, Logarithmic (Bernoulli), Negative exponential, or Quadratic).

pi=0.5
tau=20
ff=0 %Use ff=39.9999999 as ff=40
m=1
ef=2
te0=10*10000*1.0
te1h=7.5*10000*1.0
te1l=2.5*10000*1.0
alpha=0.5
gamma=2
mu=1
nu=0

% [2] Defines the abatement costs functions for the old technology (0) and the
% new technology -least efficient (1h) and most efficient (1l)-, as well as
% their respective first derivatives (0d, 1hd, 1ld)

c0=te0/e0-1; % Maximum emissions with te0 are te0 & Maximum cost when e0=0 is Inf
c0d=-te0/(e0^2);
c1h=te1h/e1h-1; % Maximum emissions with te0 are te1h & Maximum cost when e1h=0 is Inf
c1hd=-te1h/(e1h^2);
c1l=te1l/e1l-1; % Maximum emissions with te0 are te1hl & Maximum cost when e1l=0 is Inf
c1ld=-te1l/(e1l^2);

% Calculates the optimal emissions depending on the technology (0, 1h, or 1l)

ope0eqn=c0d+tau==0;
ope0 = vpasolve(ope0eqn, e0, [0 te0]);
e0c=ope0;
ope1heqn=c1hd+tau==0;
ope1h = vpasolve(ope1heqn, e1h, [0 te1h]);
e1hc=ope1h;
ope1leqn=c1ld+tau==0;
ope1l = vpasolve(ope1leqn, e1l, [0 te1l]);
e1lc=ope1l;

% [3] Assign abatement costs fucntions like in [2] substituting e0, e1h and e1l 
% by e0c, e1hc and e1lc respectively. Assign sanctioning functions with e0c, e1hc 
% and e1lc keeping them equivalent to the case without technology uncertainty to allow comparability.
% Notice that f0d corresponds to the first derivative of f0 evaluated at 
% the violation level, this is, at e0c-r0 (which is not the same as de first derivative 
% of f0 evaluated at r0). It must be updated manually, and the same applies to 
% fi1hd and fi1ld

c0=te0/e0c-1;
c1h=te1h/e1hc-1;
c1l=te1l/e1lc-1;

f0=ff*m*(e0c-r0)+m*(e0c-r0)^ef;
f0d=ff*m+m*ef*e0c-m*ef*r0;
f1h=ff*m*(e1hc-r1)+m*(e1hc-r1)^ef;
f1hd=ff*m+m*ef*e1hc-m*ef*r1;
f1l=ff*m*(e1lc-r1)+m*(e1lc-r1)^ef;
f1ld=ff*m+m*ef*e1lc-m*ef*r1;

% [4] Backwards loop for several i values. Do not directly specify i values
% in i, instead, set counter=maximum value of i to be studied.

for counter=846:-1:0; % 845 threshold when pi=0.5 tau=20 ff=0 m=1 ef=2 te0=10*10000 te1h=7.5*10000 te1l=2.5*10000 gamma=2 mu=1 nu=0

i=counter;

% Old Technology

% [5] Disutility function (see Appendix A.2). Notice that dN0d and dS0d correspond to the first
% derivatives with respect to N0=(ci0+tau*ri0) and S0=(ci0+tau*ri0+fi0)
% of dN0 and dS0 respectively, and they must be updated manually at this step.

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
% substituting r0 by OT (the amount of declared emissions with the old
% technology).

OT = vpasolve(eqn, r0,[0 e0c]); %Amount of declared emissions with the old technology

Vot=e0c-OT; %Violation level with the old technology

Cot=c0; %Abatement costs with the old technology

Tot=OT*tau; %Total taxes paid for declared emissions with the old technology

Fot=ff*m*(e0c-OT)+m*(e0c-OT)^ef; %Sanction for the violation level with the old technology

% Shows total costs to get an idea about how large they are numerically
% (the independent variable).

CostsOT=Cot+Tot+Fot

% [7] Disutility assuming the old technology. Make sure that the disutility
% functions multiplying (1-pi) and pi are the same as the ones written
% in [5] but substituting c0, tau*r0 and f0 by Cot, Tot and Fot respectively.

Dot=(1-pi)*((gamma-1)*(mu*(Cot+Tot)/(gamma-1)-nu)^gamma)+pi*((gamma-1)*(mu*(Cot+Tot+Fot)/(gamma-1)-nu)^gamma)

% Least Efficient New Technology (higher abatement cost than with the alternative New Technology)

% [8] Follow the same reasoning as in [5] when setting these disutility
% functions and keep them like the ones in [5] adding i to the total costs
% and being careful to work with 1h functions

dN1h = (gamma-1)*(mu*(c1h+tau*r1+i)/(gamma-1)-nu)^gamma;
dN1hd = gamma*mu*(mu*(c1h+tau*r1+i)/(gamma-1)-nu)^(gamma-1);
dS1h = (gamma-1)*(mu*(c1h+tau*r1+f1h+i)/(gamma-1)-nu)^gamma;
dS1hd = gamma*mu*(mu*(c1h+tau*r1+f1h+i)/(gamma-1)-nu)^(gamma-1);

% Proposition 1: Equation to implicitly obtain declared emissions.

eqn = (dS1hd*pi*f1hd)/((1-pi)*dN1hd+pi*dS1hd) == tau;

% [9] The function 'vpasolve' solves the equation numerically. With the range
% [0 e1hc] we consider only real solutions between 0 (which constitutes the
% lower bound because firms won't declare negative emissions) and the
% optimal emissions level with the least efficient new tech (which constitutes the upper
% bound because firms won't declare more emissions than actual
% emissions). Be careful to keep Fnth the same as the sanctioning function in [3]
% substituting r1 by NTh (the amount of declared emissions with the least efficient new technology)

NTh = vpasolve(eqn, r1,[0 e1hc]); %Amount of declared emissions with the least efficient new technology

Vnth=e1hc-NTh; %Violation level with the least efficient new technology

Cnth=c1h; %Abatement costs with the least efficient new technology

Tnth=NTh*tau; %Total taxes paid for declared emissions with the least efficient new technology

Inth=i; %Fixed investment cost in the new technology

Fnth=ff*m*(e1hc-NTh)+m*(e1hc-NTh)^ef; % Sanction for the violation level with the least efficient new technology

% Shows total costs to get an idea about how large they are numerically (the independent variable)

CostsNTh=Cnth+Tnth+Inth+Fnth;

% [10] Disutility assuming the least efficient new technology. Make sure that the disutility
% functions multiplying (1-pi) and pi are the same as the ones written
% in [8] but substituting c1h, tau*r1, i and 
% f1h by Cnth, Tnth, Inth and Fnth, respectively
Dnth=(1-pi)*((gamma-1)*(mu*(Cnth+Tnth+Inth)/(gamma-1)-nu)^gamma)+pi*((gamma-1)*(mu*(Cnth+Tnth+Inth+Fnth)/(gamma-1)-nu)^gamma)

% Most Efficient New Technology (lower abatement cost than with the alternative New Technology)

% [11] Follow the same reasoning as in [5] when setting these disutility
% functions and keep them like the ones in [5] adding i to the total costs
% and being careful to work with 1l functions

dN1l = (gamma-1)*(mu*(c1l+tau*r1+i)/(gamma-1)-nu)^gamma;
dN1ld = gamma*mu*(mu*(c1l+tau*r1+i)/(gamma-1)-nu)^(gamma-1);
dS1l = (gamma-1)*(mu*(c1l+tau*r1+f1l+i)/(gamma-1)-nu)^gamma;
dS1ld = gamma*mu*(mu*(c1l+tau*r1+f1l+i)/(gamma-1)-nu)^(gamma-1);

% Proposition 1: Equation to implicitly obtain declared emissions

eqn = (dS1ld*pi*f1ld)/((1-pi)*dN1ld+pi*dS1ld) == tau;

% [12] The function 'vpasolve' solves the equation numerically. With the range
% [0 e1lc] we consider only real solutions between 0 (which consitutes the
% lower bound because firms won't declare negative emissions) and the
% optimal emissions level with the most efficient new tech (which constitutes the upper
% bound because firms won't declaring more emissions than actual
% emissions). Be careful to keep Fntl the same as the sanctioning function in [3]
% substituting r1 by NTl (the amount of declared emissions with the most efficient new technology).

NTl = vpasolve(eqn, r1,[0 e1lc]); %Amount of declared emissions with the most efficient new technology

Vntl=e1lc-NTl; %Violation level with the most efficient new technology

Cntl=c1l; %Abatement costs with the most efficient new technology

Tntl=NTl*tau; %Total taxes paid for declared emissions with the most efficient new technology

Intl=i; %Fixed investment cost in the new technology

Fntl=ff*m*(e1lc-NTl)+m*(e1lc-NTl)^ef; %Sanction for the violation level with the most efficient new technology

% Shows total costs to get an idea about how large they are numerically (the independent variable)

CostsNTl=Cntl+Tntl+Intl+Fntl;

% [13] Disutility assuming the most efficient new technology. Make sure that the disutility
% functions multiplying (1-pi) and pi are the same as the ones written
% in [11] but substituting c1l, tau*r1, i and 
% f1l by Cntl, Tntl, Intl and Fntl, respectively.

Dntl=(1-pi)*((gamma-1)*(mu*(Cntl+Tntl+Intl)/(gamma-1)-nu)^gamma)+pi*((gamma-1)*(mu*(Cntl+Tntl+Intl+Fntl)/(gamma-1)-nu)^gamma)

% Now we calculate the expected disutility with the possible new technologies

EDnt=alpha*Dnth+(1-alpha)*Dntl;

% Now we check if the disutility with the old technology is larger than the
% expected disutility with the new technology. When this is the case, the firms invest
% in the new technology. This indicates the program that an approximate I threshold 
% has been found and therefore the loop ends.

if Dot-EDnt>0
    break
end
    
end
Ar=-mu*CostsOT/(mu*CostsOT/(gamma-1)-nu); % Calculates the degree of Relative Risk Aversion (Ar) using eq (A2.6)
fprintf("Degree of Relative Risk Aversion (Ar) is = %s\n", Ar);
X=['Indifferent if the installing cost of the new technology (I) is aproximately between [',num2str(i),',',num2str(i+1),')'];
disp(X)
fprintf('Optimal emissions with the old technology=%s\n',ope0);
eope1=alpha*ope1h+(1-alpha)*ope1l;
fprintf('Expected optimal emissions with the new technology=%s\n',eope1);
fprintf('Optimal emissions with the least efficient new technology=%s\n',ope1h);
fprintf('Optimal emissions with the most efficient new technology=%s\n',ope1l);
fprintf('Violation level with the old technology=%s\n',Vot);
eVnt=alpha*Vnth+(1-alpha)*Vntl;
fprintf('Expected violation level with the new technology=%s\n',eVnt);
fprintf('Violation level with the least efficient new technology=%s\n',Vnth);
fprintf('Violation level with the most efficient new technology=%s\n',Vntl);