function output = Property_Calculator

% this function will calculate and display the material properties from the
% orthotropic stiffness matrix.

%input = '2x2_0.5_01_AS';
%input = '7x7_0.5_01_AS';
%input = 'Hex_0.5_01_AS';
%input = 'HexInter_0.5_01_AS';
%input = 'HexInter_0.5_02_AS';
%input = 'SquareInter_0.5_01_AS';
%input = 'SquareInter_0.5_01_AS_(Test)';

input = '2x2_0.5_01_TS';

[C,~,~] = MOC(input);

S = inv(C);

x0 = zeros(1,12);

x0(1) = C(1,1) - 2*C(1,2)^2/(C(2,2) + C(2,3)); %E11
x0(9) = C(1,2)/(C(2,2) + C(2,3)); %v12
x0(6) = C(6,6); %g12
    K23 = 1/2*(C(2,2) + C(2,3));
x0(4) = 1/2*(C(2,2) - C(2,3)); %g23
x0(2) = 4*x0(4)*K23/(K23 + (1 + 4*K23*x0(9)^2/x0(1))*x0(4)); %E22
x0(7) = x0(2)/(2*x0(4)) - 1; %v23
x0(3) = x0(2); %E33
x0(5) = x0(6); %G31
x0(8) = x0(9); %v31
x0(10) = x0(7); %v32
x0(11) = x0(8); %v13
x0(12) = x0(9); %v21

options = optimset('PlotFcns',@optimplotfval,'MaxFunEvals',100000,'TolFun',1e-10);
x = fminsearch(@(x) myfun(x,S),x0,options);


output.E11 = x(1);
output.E22 = x(2);
output.E33 = x(3);
output.G23 = x(4);
output.G31 = x(5);
output.G12 = x(6);
output.v23 = x(7);
output.v31 = x(8);
output.v12 = x(9);
output.v32 = x(10);
output.v13 = x(11);
output.v21 = x(12);
end

function F = myfun(x,S)

E11 = x(1);
E22 = x(2);
E33 = x(3);
G23 = x(4);
G31 = x(5);
G12 = x(6);
v23 = x(7);
v31 = x(8);
v12 = x(9);
v32 = x(10);
v13 = x(11);
v21 = x(12);

s = zeros(6,6);

s(1,1) = [1/E11];
s(1,2) = [-v21/E22];
s(2,1) = [-v12/E11];
s(2,2) = [1/E22];
s(1,3) = [-v31/E33];
s(3,1) = [-v13/E11];
s(2,3) = [-v32/E33];
s(3,2) = [-v23/E22];
s(3,3) = [1/E33];
s(4,4) = [1/G23];
s(5,5) = [1/G31];
s(6,6) = [1/G12];



f(1)  = (S(1,1)-s(1,1));
f(2)  = (S(1,2)-s(1,2));
f(3)  = (S(2,1)-s(2,1));
f(4)  = (S(2,2)-s(2,2));
f(5)  = (S(1,3)-s(1,3));
f(6)  = (S(3,1)-s(3,1));
f(7)  = (S(2,3)-s(2,3));
f(8)  = (S(3,2)-s(3,2));
f(9)  = (S(3,3)-s(3,3));
f(10) = (S(4,4)-s(4,4));
f(11) = (S(5,5)-s(5,5));
f(12) = (S(6,6)-s(6,6));

f = abs(f);
%f(2:6) = f(2:6)*1e1;
%f(7:end) = f(7:end)*1e5;
F = mean(f);
end