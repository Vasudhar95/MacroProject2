%endogeneous variables
var z c n k y U_c;

%predetermined variables
predetermined_variables k;

%exogeneous shocks
varexo e_x;

%parameters
parameters alpha beta delta sigma zss sigmaz rhoz eta;

%parameters: steady state variables declared as parameters within Dynare
parameters css yss nss kss U_css;

%------
%Calibration
%------

%declare global parameters
global sigma beta delta alpha zss;

%------
%parameterization
%------

sigma = 2; %risk aversion
beta = 0.99; %Discount Factor
delta = 0.025; %depreciation rate
alpha = 0.36; %capital share in prouction function
sigmaz = 0.007; %TFP shock
rhoz = 0.95; %autocorrelation of TFP
zss = 1; %Steady state TFP

%----
%steady state and calibrated parameters
%------

feed = [
    11; %k
    0.60; %c
    0.5 %y
    ];


options_fsolve_silent = optimset('display', 'iter', 'maxfunevals', 100000, 'maxiter', 100000, 'tolfun', 10e-25, 'tolx, 10e-25);

[outc,fval,eflag]=fsolve(@RBC_calibration2,feed,options_fsolve_silent);

kss = out(1);
css = out(2);
yss = out(3);

% you can define additional equations and vairables (as many as you want)
investss = delta*kss;
yss = css + investss;
U_css = css^(-sigma)


%%%% Start model

model;

%%%%% insert model here %%%%%%

% Euler equation
1      =       beta*((U_css*exp(U_c(+1)))/(U_css*exp(U_c)))*(z(+1)*alpha*(kss*exp(k(+1)))^(alpha - 1) + 1 - delta_;

% Resource constraint
yss*exp(y)     =  css*exp(c) + kss*exp(k(+1)) - (1 - delta)*kss*exp(k);

% production function
yss*exp(y)   =   z*(kss*exp(k))^(alpha);

% Law of motion for capital
kss*exp(k(+1)) = (1 - delta)*kss*exp(k) + investss*dexp(invest);

% Definition of marginal utility of consumption
U_css*exp(U_c) = (css*exp(c))^(-sigma);

% End model

end;

% 5. Computation and Simulation

initval;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% State initial values : note that if the model is written with
% (xss*exp(x)) for variable x, this means that x is already expressed in
% logs , so the initial values for variables in this form are 0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

z          -1;
k          -0;
c          -0;
invest     -0;
y          -0;
U_c        -0;

end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

shocks;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% State shocks in the model and specify size of shocks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

var e_z-sigmaz^2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Check residuals : these should all be exactly zero if the model is
% written correctly

resid(1);

%Choose algerithm that Dynare uses to find the steady state

stady(solve_algo-0);

check;

%Fix an intial seef dor shocks so that every time the model is simulated 
%the same shock series is used to compute second moments

set_dynare_seed(100000);

%Perform a stochastic simulation using a first-order approximation
%for 2100 periods and using an HP filter with smoothing parameter 1600
%on the simulated series to compute simulated second moments




%%%%%%%%%%%%%%%%%%%

%%%%% modification of functional form for #3
% kss_chi = ((chi/2)*((k_t+1/k_t)-1)^2)*(k_t), where chi ? 0

