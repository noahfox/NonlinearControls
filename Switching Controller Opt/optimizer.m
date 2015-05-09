clear *; close all; clc

global best_yet initials
best_yet = 0;
consts = get_consts();

% States: x=[y,z,th,psi,dy,dz,dth,dpsi,m]'
fuel = consts.m_nofuel+0.7*consts.max.m_fuel;
initials = [10 150 0 0 0 0 0 0 fuel; % 0
    10 150 pi/2 0 0 0 0 0 fuel; % 1
    10 25 0 0 0 0 0 0 fuel; % 2
    10 25 pi/4 0 0 0 0 0 fuel; % 3
    100 25 pi/4 0 0 0 0 0 fuel; % 4
    100 450 pi 0 0 0 0 0 fuel; % 5
    10 450 pi 0 0 0 0 0 fuel; % 6
    100 1500 pi 0 0 0 0 0 fuel; % 7
    100 1500 pi/2 0 0 0 0 0 fuel; % 8
    100 700 pi/2 0 0 0 0 0 fuel]; % 9


opts.MaxFunEvals  = 50000;
opts.TolX = 1e-6;
opts.StopFitness = 1e-6;
opts.CMA.active = 1;
opts.Restarts = 0;
opts.LogPlot = 'off';
opts.StopOnStagnation = 'on';
opts = cmaes('defaults', opts);

% sigma = [.001 .001 10 .001 .001 .001 1 .001 0.001 0.1 0.1]';
sigma = 0.01;

% x0 = [0.001 0.001 100 0.001 0.001 0.001 10 0.001 0.001 .5 2.0... 
%       0.1 100 10 0.000001 0.001 0.001 .001 0.001 0.001 .1 10];
% x0 = [ -0.7544 1.1544 0 0 -5.5728 8.5728 0 0 -0.1575];
x0 = [-0.684995189331813,1.15059724593106,-0.102287336164053,0.0697255362865406,-5.61834282543692,8.46906736340270,-0.0865564201780617,-0.00516796088993788,-0.150346777975641];

[XMIN,FMIN,COUNTEVAL,STOPFLAG,OUT,BESTEVER] = cmaes('costfun',x0,sigma,opts);
