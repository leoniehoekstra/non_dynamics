function out = colpitts
out{1} = @init;
out{2} = @fun_eval;
out{3} = @jacobian;
out{4} = @jacobianp;
out{5} = @hessians;
out{6} = @hessiansp;
out{7} = @der3;
out{8} = [];
out{9} = [];

% --------------------------------------------------------------------------
function dydt = fun_eval(t, kmrgd, par_g, par_Q)
e = exp(-kmrgd(2));
scale = 2*par_g/par_Q;
dydt = zeros(3,1);
dydt(1) = scale*(1 - e + kmrgd(3));
dydt(2) = scale*kmrgd(3);
dydt(3) = -(par_Q/(4*par_g))*(kmrgd(1) + kmrgd(2)) - (1/par_Q)*kmrgd(3);

% --------------------------------------------------------------------------
function [tspan, y0, options] = init
handles = feval(@colpitts);
y0 = [0, 0, 0];
options = odeset('Jacobian', handles(3), 'JacobianP', handles(4), ...
    'Hessians', handles(5), 'HessiansP', handles(6), 'Der3', handles(7));
tspan = [0 10];

% --------------------------------------------------------------------------
function jac = jacobian(t, kmrgd, par_g, par_Q)
e = exp(-kmrgd(2));
scale = 2*par_g/par_Q;
jac = zeros(3,3);
jac(1,1) = 0;
jac(1,2) = scale*e;
jac(1,3) = scale;
jac(2,1) = 0;
jac(2,2) = 0;
jac(2,3) = scale;
jac(3,1) = -(par_Q/(4*par_g));
jac(3,2) = -(par_Q/(4*par_g));
jac(3,3) = -1/par_Q;

% --------------------------------------------------------------------------
function jacp = jacobianp(t, kmrgd, par_g, par_Q)
e = exp(-kmrgd(2));
common = 1 - e + kmrgd(3);
jacp = zeros(3,2);
jacp(1,1) = (2/par_Q)*common;
jacp(1,2) = -(2*par_g/par_Q^2)*common;
jacp(2,1) = (2/par_Q)*kmrgd(3);
jacp(2,2) = -(2*par_g/par_Q^2)*kmrgd(3);
jacp(3,1) = (par_Q/(4*par_g^2))*(kmrgd(1) + kmrgd(2));
jacp(3,2) = -(1/(4*par_g))*(kmrgd(1) + kmrgd(2)) + (1/par_Q^2)*kmrgd(3);

% --------------------------------------------------------------------------
function hess = hessians(t, kmrgd, par_g, par_Q)
e = exp(-kmrgd(2));
scale = 2*par_g/par_Q;
hess = zeros(3,3,3);
hess1 = zeros(3);
hess1(2,2) = -scale*e;
hess(:,:,1) = hess1;
hess(:,:,2) = zeros(3);
hess(:,:,3) = zeros(3);

% --------------------------------------------------------------------------
function hessp = hessiansp(t, kmrgd, par_g, par_Q)
e = exp(-kmrgd(2));
hessp = zeros(3,3,2);
hessp1 = zeros(3);
hessp1(2,2) = (2/par_Q)*e;
hessp2 = zeros(3);
hessp2(2,2) = -(2*par_g/par_Q^2)*e;
hessp(:,:,1) = hessp1;
hessp(:,:,2) = hessp2;

% --------------------------------------------------------------------------
function tens3 = der3(t, kmrgd, par_g, par_Q)
e = exp(-kmrgd(2));
scale = 2*par_g/par_Q;
tens3 = zeros(3,3,3,3);
tens3(1,2,2,2) = scale*e;
