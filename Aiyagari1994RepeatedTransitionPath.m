% Example of computing the repeated reveal of a general eqm transition path for the model of Aiyagari (1994).
%
% These codes set up and solve an extended Aiyagari (1994) model with endogenous labor 
% and adding labor income tax, for a given parametrization. 
% They then show how to solve for the general equilibrium transition path in reposonse 
% to a series of 'surprise' changes in the path of taxes.
% So in period 1 everyone wakes up expecting that the tax rate will
% increase in 10 periods time and then remain at that new higher level forever
% (the first transitition path reveal), but then when we get to period five 
% everyone 'wakes up' again and discovers that actually the taxes increased 
% in period five and will remain at this new higher rate forever (the
% second transition path reveal)
%
% To solve we essentially just first solve for the first transition path
% reveal. Then we take where the agent distribution was in period five of
% this and use that as the starting point to solve for the second
% transition path reveal.
%
% One clarifying comment: when solving transition paths it is important that the
% economy ends up in a stationary equilbrium, but where the economy starts
% is irrelevant; all we need to know about the start of the transitition
% path is what the agent distribution is (there is no need for this
% starting point of a transition to be anything to do with any kind of
% equilibrium concept). If you look at the inputs to the command for
% solving the transition path you will see that the agent distribution is
% only input that is related to how the economy was prior to the transition
% path being revealed.
%
% If you wanted to do three or more reveals, you can just easily extend
% this concept: solve first reveal, find agent distribution for start of
% second reveal, solve second reveal, find agent distribution for start of
% third reveal, solve third reveal, etc.
%
% Note: I am lazy solving the stationary initial and final equilbria. We
% know an equilibrium relationship that will hold between r and w, and we
% could use this so that we just find r to clear the capital market general
% eqm condition (and not need to solve w to clear the labor market general
% eqm condition). This would be faster, but this problem is simple enough
% anyway, so I just take the lazy option.
%
% Note: Would be better with n_d=101, n_k=2^10, and transpathoptions.tolerance=10^(-5)
% But this is just an illustration/example, so used something a little
% smaller so it solves a bit easier/faster.

%% Set some basic variables

n_d=51;
n_k=2^9;%2^9;
n_z=7; %21;

%Parameters
Params.beta=0.96; %
Params.alpha=0.36;
Params.delta=0.08; % depreciation rate
Params.mu_c=2;% utility, curvature of consumption
Params.sigma=0.2;
Params.rho=0.6;
Params.psi=10;
Params.mu_l=2; % utility, curvature of leisure

% Tax before and after reform (note that the different transition paths
% will just be about when the tax rate changes, to keep things simple I use
% the same final tax rate for both, but this is not important to anything
% about how the codes function)
Params.tau_initial=0.1; % labor income tax
Params.tau_final=0.12;

%% Since this example is intended to show working of transtion paths, make it verbose (print output)
transpathoptions.verbose=1;

%% Set up the exogenous shock process
% Create markov process for the exogenous labor efficiency units, z.
[z_grid,pi_z]=discretizeAR1_FarmerToda(0,Params.rho,sqrt((1-Params.rho^2)*Params.sigma^2),n_z);
% Note: sigma is standard deviations of z, input needs to be standard deviation of the innovations
% Because (log of) z is AR(1), the variance of the innovations is (1-rho^2)*sigma^2
z_grid=exp(z_grid);

%% Grids

% Grid for fraction of time worked
h_grid=linspace(0,1,n_d)';

% In the absence of idiosyncratic risk, the steady state equilibrium is given by
r_ss=1/Params.beta-1;
K_ss=((r_ss+Params.delta)/Params.alpha)^(1/(Params.alpha-1)); %The steady state capital in the absence of aggregate uncertainty.

% Set grid for asset holdings
k_grid=10*K_ss*(linspace(0,1,n_k).^3)'; % linspace ^3 puts more points near zero, where the curvature of value and policy functions is higher and where model spends more time

% Bring model into the notational conventions used by the toolkit
d_grid=h_grid;
a_grid=k_grid;
% pi_z;
% z_grid

n_a=n_k;
% n_z

%%
% Create functions to be evaluated
FnsToEvaluate.K = @(h,aprime,a,z) a; % aggregate assets (which is this periods state)
FnsToEvaluate.L = @(h,aprime,a,z) h*z; % labor supply (in labor efficiency units)
FnsToEvaluate.TaxRevenue = @(h,aprime,a,z,w,tau) tau*w*h*z; % aggregate assets (which is this periods state)

% Now define the functions for the General Equilibrium conditions
    % Should be written as LHS of general eqm eqn minus RHS, so that the closer the value given by the function is to 
    % zero, the closer the general eqm condition is to holding.
GeneralEqmEqns.LaborMarket = @(w,K,L,alpha) w-(1-alpha)*(K^alpha)*(L^(-alpha)); % wage equals marginal product of (effective) labor
GeneralEqmEqns.CapitalMarket = @(r,K,L,alpha,delta) r-(alpha*(K^(alpha-1))*(L^(1-alpha))-delta); % interest rate equals marginal product of capital
GeneralEqmEqns.GovBudgetConstraint = @(G,TaxRevenue) G-TaxRevenue; % Government runs a balanced budget
% Inputs can be any parameter, price, or aggregate of the FnsToEvaluate

%%
DiscountFactorParamNames={'beta'};

ReturnFn=@(h,aprime, a, z,r,w,mu_c,mu_l,psi,tau)... 
    Aiyagari1994labour_ReturnFn(h,aprime, a, z,r,w,mu_c,mu_l,psi,tau);

%%

% Use the toolkit to find the equilibrium price index
GEPriceParamNames={'w','r','G'};
% Set initial value for interest rates (Aiyagari proves that with idiosyncratic
% uncertainty, the eqm interest rate is limited above by it's steady state value
% without idiosyncratic uncertainty, that is that r<r_ss).
Params.w=1.18;
Params.r=0.038; % interest rate
Params.G=0.024; % Government spending (just burnt, unrealistic but keeps model simple), this is just an initial guess, it will be determined in general equilibrium

% Note, for the stationary equilibrium we could actually skip the w and
% LaborMarket and just use an expression that ties r and w in stationary
% equilibrium (as the basic codes Aiyagari (1994)). I don't do this just as
% it won't work in the transition path and so it is conceptually easier to
% follow (if a slight waste computationally) if I just use all three even
% for the stationary general equilbrium.

%% Compute the initial general equilibrium
Params.tau=Params.tau_initial;

% Solve for the stationary general equilbirium
vfoptions=struct(); % Use default options for solving the value function (and policy fn)
simoptions=struct(); % Use default options for solving for stationary distribution
heteroagentoptions.verbose=1; % verbose means that you want it to give you feedback on what is going on

fprintf('Calculating price vector corresponding to the stationary general eqm \n')
[p_eqm_init,~,GeneralEqmCondn_init]=HeteroAgentStationaryEqm_Case1(n_d, n_a, n_z, 0, pi_z, d_grid, a_grid, z_grid, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, [], [], [], GEPriceParamNames,heteroagentoptions, simoptions, vfoptions);

p_eqm_init % The equilibrium values of the GE prices
% Note: GeneralEqmCondn_init will be essentially zero, it is the value of the general equilibrium equation

% For the transition path we will need the initial agents distribution
Params.w=p_eqm_init.w;
Params.r=p_eqm_init.r;
Params.G=p_eqm_init.G;
[~,Policy_init]=ValueFnIter_Case1(n_d,n_a,n_z,d_grid,a_grid,z_grid, pi_z, ReturnFn, Params, DiscountFactorParamNames, [], vfoptions);

StationaryDist_init=StationaryDist_Case1(Policy_init,n_d,n_a,n_z,pi_z);

% Following line is just a check
AggVars_init=EvalFnOnAgentDist_AggVars_Case1(StationaryDist_init, Policy_init, FnsToEvaluate, Params, [], n_d, n_a, n_z, d_grid, a_grid, z_grid);

%% Compute the final general equilbrium
Params.tau=Params.tau_final;

% Note: if the change in parameters affected pi_z this would need to be recalculated here.

disp('Calculating price vector corresponding to the final stationary eqm')
[p_eqm_final,~,GeneralEqmCondn_final]=HeteroAgentStationaryEqm_Case1(n_d, n_a, n_z, 0, pi_z, d_grid, a_grid, z_grid, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, [], [], [], GEPriceParamNames,heteroagentoptions, simoptions, vfoptions);

p_eqm_final % The equilibrium values of the GE prices
% Note: GeneralEqmCondn_final will be essentially zero, it is the value of the general equilibrium equation

% For the transition path we will need the final value function
Params.r=p_eqm_final.r;
[V_final,Policy_final]=ValueFnIter_Case1(n_d,n_a,n_z,d_grid,a_grid,z_grid, pi_z, ReturnFn,Params, DiscountFactorParamNames,[],vfoptions);

StationaryDist_final=StationaryDist_Case1(Policy_final,n_d,n_a,n_z,pi_z);
AggVars_final=EvalFnOnAgentDist_AggVars_Case1(StationaryDist_final, Policy_final, FnsToEvaluate, Params, [], n_d, n_a, n_z, d_grid, a_grid, z_grid);

% surf(k_grid*ones(1,n_s),ones(n_a,1)*s_grid',V_final)

% Just double-check that the initial and final equilibria seem to be accurate
GeneralEqmCondn_init
GeneralEqmCondn_final

%% Compute the transition path for the 'first reveal'
% For this we need the following extra objects: PricePathOld, PriceParamNames, ParamPath, ParamPathNames, T, V_final, StationaryDist_init
% (already calculated V_final & StationaryDist_init above)

% Number of time periods to allow for the transition (if you set T too low
% it will cause problems, too high just means run-time will be longer).
T=150

% We want to look at a one off unanticipated change of beta. ParamPath & PathParamNames are thus given by
ParamPath1.tau=[Params.tau_initial*ones(10,1);Params.tau_final*ones(T-10,1)]; % For each parameter that changes value, ParamPath is matrix of size T-by-1
% (the way ParamPath is set is designed to allow for a series of changes in the parameters)

% We need to give an initial guess for the price path on interest rates
% (this is deliberately not a good guess, so you can see that the transition path can be found)
PricePath0.w=[linspace(p_eqm_init.w, p_eqm_final.w, floor(T/2))'; p_eqm_final.w*ones(T-floor(T/2),1)]; % For each price, PricePath0 is matrix of size T-by-1
PricePath0.r=[linspace(p_eqm_init.r, p_eqm_final.r, floor(T/2))'; p_eqm_final.r*ones(T-floor(T/2),1)]; % For each price, PricePath0 is matrix of size T-by-1
PricePath0.G=[linspace(p_eqm_init.G, p_eqm_final.G, floor(T/2))'; p_eqm_final.G*ones(T-floor(T/2),1)]; % For each price, PricePath0 is matrix of size T-by-1

% General equilibrium conditions (for the transition path)
TransPathGeneralEqmEqns.LaborMarket = @(w,K,L,alpha) w-(1-alpha)*(K^alpha)*(L^(-alpha)); % wage equals marginal product of (effective) labor
TransPathGeneralEqmEqns.CapitalMarket = @(r,K,L,alpha,delta) r-(alpha*(K^(alpha-1))*(L^(1-alpha))-delta); % interest rate equals marginal product of capital
TransPathGeneralEqmEqns.GovBudgetConstraint = @(G,TaxRevenue) G-TaxRevenue; % Government runs a balanced budget
% Note: For this model the transition path has the same general equilibrium conditions as the stationary equilibrium, but this will not always be true for more complex models.

transpathoptions.GEnewprice=3;
% Need to explain to transpathoptions how to use the GeneralEqmEqns to
% update the general eqm transition prices (in PricePath).
transpathoptions.GEnewprice3.howtoupdate={... % a row is: GEcondn, price, add, factor
    'LaborMarket','w',0,0.1;...
    'CaptialMarket','r',0,0.1;...
    'GovBudgetConstraint','G',0,0.1;...
    }; % CaptialMarket GE condition will be positive if r is too big, so subtract
% Note: the update is essentially new_price=price+factor*add*GEcondn_value-factor*(1-add)*GEcondn_value
% Notice that this adds factor*GEcondn_value when add=1 and subtracts it what add=0
% A small 'factor' will make the convergence to solution take longer, but too large a value will make it 
% unstable (fail to converge). Technically this is the damping factor in a shooting algorithm.


% Now just run the TransitionPath_Case1 command (all of the other inputs
% are things we had already had to define to be able to solve for the
% initial and final equilibria)
transpathoptions.weightscheme=1;
transpathoptions.verbose=1;

transpathoptions.tolerance=10^(-4); % default is 10^(-5).

PricePath1=TransitionPath_Case1(PricePath0, ParamPath1, T, V_final, StationaryDist_init, n_d, n_a, n_z, pi_z, d_grid,a_grid,z_grid, ReturnFn, FnsToEvaluate, TransPathGeneralEqmEqns, Params, DiscountFactorParamNames, transpathoptions);

% Before we can solve the second reveal transition path we need the agent
% distribution in period 5 (which is when the second reveal occours).
[VPath,PolicyPath]=ValueFnOnTransPath_Case1(PricePath1, ParamPath1, T, V_final, Policy_final, Params, n_d, n_a, n_z, pi_z, d_grid, a_grid,z_grid, DiscountFactorParamNames, ReturnFn, transpathoptions, vfoptions, simoptions);
AgentDistPath=AgentDistOnTransPath_Case1(StationaryDist_init, PolicyPath,n_d,n_a,n_z,pi_z,T,transpathoptions,simoptions);
% Get the agent distribution in period 5 (the last dimension of AgentDistPath indexes the transition path time periods)
AgentDist_Period5ofFirstTransition=AgentDistPath(:,:,5); % (a,z,t)

%% Compute the transition path for the second reveal
% We just need to redefine the ParamPath and the initial agent
% distribution, everything else is the same.

ParamPath2.tau=Params.tau_final*ones(T,1); % For each parameter that changes value, ParamPath is matrix of size T-by-1
% Note: this is period 5 to T+5 under the timing of the first reveal transition path

PricePath2=TransitionPath_Case1(PricePath0, ParamPath2, T, V_final, AgentDist_Period5ofFirstTransition, n_d, n_a, n_z, pi_z, d_grid,a_grid,z_grid, ReturnFn, FnsToEvaluate, TransPathGeneralEqmEqns, Params, DiscountFactorParamNames, transpathoptions);

% Comment: This second reveal transition includes 5 extra periods beyond the end of the
% previous one. This is irrelevant (as with or without the 5 'extra' periods
% we have already converged into the final stationary equilibrium as long as 
% T is large enough to begin), and we could delete the last five periods from 
% the second reveal transition path if we wanted.

%% Plot the interest rate and wage over the entirety of both reveals
subplot(2,1,1); plot([PricePath1.r(1:4);PricePath2.r(5:end)]);
title('Interest rate')
subplot(2,1,2); plot([PricePath1.w(1:4);PricePath2.w(5:end)]);
title('Wage')


