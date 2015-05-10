% Function to setup and perform any one-time computations
% Input parameters
%   x - state of the rocket, x=[y,z,th,psi,dy,dz,dth,dpsi,m]^T
%   consts - structure that contains various system constants
% Output parameters
%   ctrl  -  any student defined control parameters
function ctrl = student_setup(x0, consts, qr)
global done
done = 0;

% qr = [0.00204073753336390; -0.000908019845082666;94.6797088671666;-0.000565059278701018;-0.000481162874801398;0.00119017653895017;8.89500774975685;0.00379848720593460;-0.000741772844496262;0.529571440663460;1.85372872712185]';

[~,P] = lqrmaker(qr);

K = [-1.04716829144030,1.54078107984502,-0.0749218299462610,0.0563990208890100,-10.5290588084430,13.4833642227773,-0.0944399610180675,-0.0444385584797427,-0.141424923092204;
    0.156183172683227,0.0123302410978552,-69.3705867184020,123.748316667828,1.93832562254661,0.0383528368305634,-103.803019040439,68.2406097022938,-0.00163293034371243];

syms x1 x2 x3 x4 x5 x6 x7 x8 x9 real

f_vec = [   x5
    x6
    x7
    x8
    0
    -consts.g
    0
    0
    0] ;

g_vec = [0,    0 ;
    0,    0 ;
    0,    0 ;
    0,    0 ;
    -consts.gamma*sin(x4+x3)/x9,    0 ;
    consts.gamma*cos(x4+x3)/x9,    0 ;
    -consts.L*consts.gamma*sin(x4)/consts.J,    0 ;
    0, 1/consts.JT ;
    -1,    0] ;

x = [x1 x2 x3 x4 x5 x6 x7 x8 x9]';
V = x'*P*x;
LfV(x1,x2,x3,x4,x5,x6,x7,x8,x9) = jacobian(V,x)*f_vec;
LgV(x1,x2,x3,x4,x5,x6,x7,x8,x9) = jacobian(V,x)*g_vec;
ctrl.LfVfun = matlabFunction(LfV);
ctrl.LgVfun = matlabFunction(LgV);
ctrl.K = K;
ctrl.p = x0;
end

function [K,P] = lqrmaker(qr)
consts = get_consts();
global target
syms y z th ps dy dz dth dpsi m ft T real

f_vec = [   dy
            dz
            dth
            dpsi
            0
            -consts.g
            0
            0
            0] ;

g_vec = [0,    0 ;
        0,    0 ;
        0,    0 ;
        0,    0 ;
        -consts.gamma*sin(ps+th)/m,    0 ;
        consts.gamma*cos(ps+th)/m,    0 ;
        -consts.L*consts.gamma*sin(ps)/consts.J,    0 ;
        0, 1/consts.JT ;
        -1,    0] ;

x = [y z th ps dy dz dth dpsi m];
u = [ft; T];

sys = f_vec + g_vec*u;

ft_0 = (consts.m_nofuel*consts.g)/consts.gamma;

target = [0 consts.L 0 0 0 0 0 0 consts.max.m_fuel];
eq = [0 consts.L  0 0 0 0 0 0 consts.max.m_fuel 1.7168 0];
eq_vars = [y z th ps dy dz dth dpsi m ft T];

A_sym = jacobian(sys,x);
B_sym = jacobian(sys,u);

A = eval(subs(A_sym,eq_vars,eq));
B = eval(subs(B_sym,eq_vars,eq));

Q = diag(abs(qr(1:9)));
R = diag(abs(qr(10:11)));
[K,P] = lqr(A,B,Q,R);
end