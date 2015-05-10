clear *; close all; clc

global best_yet initials
best_yet = 0;
consts = get_consts();

% States: x=[y,z,th,psi,dy,dz,dth,dpsi,m]'
fuel = consts.m_nofuel+consts.max.m_fuel;
initials = [10 150 0.1 0 0 0 0 0 fuel;
            50 1500 0.2 0 0 0 0 0 fuel;
            100 25 0 0 0 0 0 0 fuel;
            100 500 pi/2 0 0 0 0 0 fuel;
            10 1500 179*pi/180 0 0 0 0 0 fuel];


opts.MaxFunEvals  = 50000;
opts.TolX = 1e-6;
opts.StopFitness = 1e-6;
opts.CMA.active = 1;
opts.Restarts = 0;
opts.LogPlot = 'off';
opts.StopOnStagnation = 'on';
opts = cmaes('defaults', opts);

sigma = [.01 .01 30 .01 .01 .01 3 .01 .01 .2 .6]';

x0 = [0.00204073753336390; -0.000908019845082666;94.6797088671666;-0.000565059278701018;-0.000481162874801398;0.00119017653895017;8.89500774975685;0.00379848720593460;-0.000741772844496262;0.529571440663460;1.85372872712185]';

[XMIN,FMIN,COUNTEVAL,STOPFLAG,OUT,BESTEVER] = cmaes('costfun',x0,sigma,opts);
