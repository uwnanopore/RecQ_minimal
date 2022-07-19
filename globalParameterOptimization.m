%% load the data and set number of monte carlo trials 

load('forceDependentStatsStructure')
% loads the variable 'finalStats' which contains all data and all
% analysis products. Any variable with 'd' in front of it is an error bar,
% e.g. finalStats.velocity = measurements of velocity
%      finalStats.dvelocity = error on velocity
nMonte = 100;

%% Initial conditions

% 
guessVals1 = [300 2 200 10 180 100 100 0.01 0.01 0.005 0.005];
guessVals2 = [1 1 1 1 1 1 1 0.001 0.001 0.001 0.001];
guessVals3 = [10 10 10 10 10 10 10 0.005 0.005 0.005 0.005];
guessVals4 = [1000 1000 1000 1000 1000 1000 1000 0.03 0.03 0.03 0.03];
guessVals5 = [200 10 400 2 100 200 200 0.01 0.01 0.005 0.005];
guessVals6 = [100 100 100 100 100 100 100 0.02 0.02 0.01 0.01];
guessVals7 = [1000 5 1000 5 500 100 500 0.03 0.02 0.01 0.01];
guessVals8 = guessVals7/2;
guessVals9 = [500 50 500 50 500 500 500 0.001 0.001 0.001 0.001];
guessVals10 = [500 50 300 50 500 500 500 0.03 0.03 0.03 0.03];
guessVals11 = [200 50 200 50 200 200 200 0.001 0.03 0.001 0.03];
guessVals12 = [300 50 500 50 300 1000 500 0.03 0.001 0.03 0.001];
guessVals13 = [100 50 100 50 500 100 300 0.02 0.02 0.02 0.02];
guessVals14 = [100 0.1 100 0.1 100 100 100 0.015 0.0015 0.0015 0.0015];
guessVals15 = [50 10 500 50 1000 500 1000 0.015 0.0015 0.0015 0.0015];

guessVals = {guessVals15}; % choose a subset of initial conditions to run
%%  Monte-carlo jitter the data
finalStats.monte.tauDep_ff = cell(1,nMonte); %set up cell array of monte carlo trials
finalStats.monte.tauDep_fb = cell(1,nMonte);
finalStats.monte.tauDep_bf = cell(1,nMonte);
finalStats.monte.tauIndep_ff = cell(1,nMonte);
finalStats.monte.tauIndep_bf = cell(1,nMonte);
finalStats.monte.tauIndep_fb = cell(1,nMonte);
finalStats.monte.pbackDep = cell(1,nMonte);
finalStats.monte.pbackIndep = cell(1,nMonte);
finalStats.monte.params = zeros(length(guessVals),nMonte,11); %% setup output parameters
finalStats.monte.paramsWithVelocityConstraint = zeros(length(guessVals),nMonte,11); %% set output parameters with data constrained to total velocity
finalStats.monte.likelihood = zeros(length(guessVals),nMonte); 

for ii = 1:nMonte
    finalStats.monte.tauDep_ff{ii} = finalStats.tauDep_ff + randn(1,12).*finalStats.dtauDep_ff;
    finalStats.monte.tauDep_bf{ii} =  finalStats.tauDep_bf + randn(1,6).*finalStats.dtauDep_bf;
    finalStats.monte.tauDep_fb{ii} =  finalStats.tauDep_fb + randn(1,6).*finalStats.dtauDep_fb;
    finalStats.monte.tauIndep_ff{ii} =  finalStats.tauIndep_ff + randn(1,12).*finalStats.dtauIndep_ff;
    finalStats.monte.tauIndep_bf{ii} =  finalStats.tauIndep_bf + randn(1,6).*finalStats.dtauIndep_bf;
    finalStats.monte.tauIndep_fb{ii} =  finalStats.tauIndep_fb + randn(1,6).*finalStats.dtauIndep_fb;
    finalStats.monte.pbackDep{ii} =  finalStats.pbackDep + randn(1,12).*finalStats.dpbackDep;
    finalStats.monte.pbackIndep{ii} =  finalStats.pbackIndep + randn(1,12).*finalStats.dpbackIndep; 
end
finalStats.README = ' alpha (km1) beta(km2) gamma (k1) delta(k2)';
finalStats.identity = {'k1 = x1', 'km1 = x2', 'k2 = x3', 'km2 = x4',...
    'kd = x5', 'kh = x6', 'kp = x7', 'alpha = x8', 'beta = x9',...
    'gamma = x10', 'delta = x11'}'; % variable names 
% in the final version of the paper these change
% expressions are defined below. km1 and k1 refer to ADP-bound
% conformational changes and km2 and k2 refer to ATP-bound conformational
% changes
%% Set up anonymous functions and cost function
% the first entry in each anonymous function is the applied voltage
velocityTotal = @(v,k1,km1,k2,km2,kd,kh,kp,a,b,g,d) ...
   ( ((k2.*exp(g.*v) + kh + km2.*exp(-b.*v)))./(k2.*exp(g.*v).*kh) + ...
    ((k1.*exp(d.*v) + kd + km1.*exp(-a.*v)))./(k1.*exp(d.*v).*kd) + 1./kp).^-1; 

p_D = @(v,km1,kd,a) km1./(km1 + kd.*exp(a.*v)); % _D = ATP-dependent step
tff_D = @(v,km1,kd,a,k2,d) (kd + k2.*exp(d.*v) + km1.*exp(-a.*v))./(k2.*exp(d.*v).*(kd+km1.*exp(-a.*v))); 
tbf_D = @(v,km1,kd,a) 1./(kd + km1.*exp(-a.*v));
tfb_D = @(v,k2,d) 1./(k2.*exp(d.*v));

p_I = @(v,km2,kh,b) km2./(km2 + kh.*exp(b.*v)); % _I = ATP-independent step
tff_I = @(v,km2,kh,b,k1,g,kp) 1./(kh + km2.*exp(-b.*v)) + 1./kp +  1./(k1.*exp(g.*v));
tbf_I = @(v,km2,kh,b) 1./(kh + km2.*exp(-b.*v));
tfb_I = @(v,k1,g) 1./(k1.*exp(g.*v));
 
lowerz = [0 0 0 0 0 0 0 0 0 0 0]; % Set lower and upper bounds for constrained fit
upperz = [Inf Inf Inf Inf Inf Inf Inf 0.025 0.025 0.025 0.025];

opts = optimset('maxfunevals',10^6,'maxiter',10^6); % options for fminsearch
% opts = optimset('PlotFcns',@optimplotfval);
for kk = 1:length(guessVals)
    for jj = 1:nMonte
    disp(jj)

    chiP_D = @(x)  sum((finalStats.monte.pbackDep{jj} - p_D(finalStats.voltage,x(2),x(5),x(8))).^2./(finalStats.dpbackDep).^2);
    chiff_D = @(x) sum((finalStats.monte.tauDep_ff{jj} - tff_D(finalStats.voltage,x(2),x(5),x(8),x(3),x(11))).^2./(finalStats.dtauDep_ff).^2);
    chibf_D = @(x) sum((finalStats.monte.tauDep_bf{jj} - tbf_D(finalStats.voltage(1:6),x(2),x(5),x(8))).^2./(finalStats.dtauDep_bf).^2);
    chifb_D = @(x) sum((finalStats.monte.tauDep_fb{jj} - tfb_D(finalStats.voltage(1:6),x(3),x(11))).^2./(finalStats.dtauDep_fb).^2); 
    
    chiVelocity = @(x) sum( ((finalStats.velocity-velocityTotal(finalStats.voltage,x(1),x(2),x(3),x(4),x(5),x(6),x(7),x(8),x(9),x(10),x(11)))./finalStats.dvelocity).^2);

    chiP_I = @(x) sum((finalStats.monte.pbackIndep{jj} - p_I(finalStats.voltage,x(4),x(6),x(9))).^2./(finalStats.dpbackIndep).^2);
    chiff_I = @(x) sum((finalStats.monte.tauIndep_ff{jj} - tff_I(finalStats.voltage,x(4),x(6),x(9),x(1),x(10),x(7))).^2./(finalStats.dtauIndep_ff).^2);
    chibf_I = @(x) sum((finalStats.monte.tauIndep_bf{jj} - tbf_I(finalStats.voltage(1:6),x(4),x(6),x(9))).^2./(finalStats.dtauIndep_bf).^2);
    chifb_I = @(x) sum((finalStats.monte.tauIndep_fb{jj} - tfb_I(finalStats.voltage(1:6),x(1),x(10))).^2./(finalStats.dtauIndep_fb*0+mean(finalStats.dtauIndep_fb).^2));

    chiTotal = @(x) chiP_D(x) + chiff_D(x) + chibf_D(x) + chifb_D(x)...
                  + chiP_I(x) + chiff_I(x) + chibf_I(x) + chifb_I(x);
    chiTotalWithVelocityConstraint = @(x) chiP_D(x) + chiff_D(x) + chibf_D(x) + chifb_D(x) + ...
                                            chiP_I(x) + chiff_I(x) + chibf_I(x) + chifb_I(x) + ...
                                            chiVelocity(x); % fminsearch only works in a single vector output x which is why we need to write the 

    finalStats.monte.params(kk,jj,:) = ...
        fminsearchbnd(chiTotal,guessVals{kk},lowerz,upperz,opts); % perform parameter optimization without using the velocity
    
    [finalStats.monte.paramsWithVelocityConstraint(kk,jj,:),finalStats.monte.likelihood(kk,jj),exitFlag] = ...
        fminsearchbnd(chiTotalWithVelocityConstraint,guessVals{kk},lowerz,upperz,opts); % perform parameter optimization using the velocity 
    % - this ultimately stabilizes the solutions 
%         disp(exitFlag)
    end
end




