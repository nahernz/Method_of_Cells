function output = CCM(input)

% output = CCM(input)
%
% CCM computes the composite properties by the composite cylinder model 
% from the component level properties.
% 
% INPUT is a struct with the following fields
%
% Fiber
%   :E23f
%   :E12f
%   :v23f
%   :v12f
%   :G23f
%   :G12f
% Matrix
%   :Em
%   :vm
% Volume Fraction
%   :Vf
%
% OUTPUT is a struct with the following fields
%   :E11
%   :E22
%   :v12
%   :v23
%   :G12
%   :G23
%   :K23

% sample input

input.E11f = 251000; % MPa
input.E22f =  40400; % MPa
input.v23f =  0.256; 
input.v12f =  0.321; 
input.G23f =  16080; % MPa
input.G12f =  30700; % MPa

input.Em  =   3310; % MPa
input.vm  =  0.318;

input.Vf  =    0.5;

% build the fiber compliance matrix
S11 = 1/input.E11f;
S22 = 1/input.E22f;
S12 = -input.v12f/input.E11f;
S23 = -input.v23f/input.E22f;
S44 = 1/input.G23f;
S66 = 1/input.G12f;

S = [S11 S12 S12  0   0   0  ;
     S12 S22 S23  0   0   0  ;
     S12 S23 S22  0   0   0  ;
      0   0   0  S44  0   0  ;
      0   0   0   0  S66  0  ;
      0   0   0   0   0  S66];

C = inv(S);

% find the bulk moduli
input.K23f = 1/2*(C(2,2)+C(2,3));
input.Km  = input.Em/(3*(1-2*input.vm));
input.Gm  = input.Em/(2*(1+input.vm));

fprintf('\ninput =\n\n')
disp(input)

E11f  = input.E11f;
E22f  = input.E22f;
v23f = input.v23f;
v12f = input.v12f; 
G23f = input.G23f;
G12f = input.G12f;
K23f = input.K23f;
Em   = input.Em;
vm   = input.vm;
K23m = input.Km;
Gm   = input.Gm;
Vf   = input.Vf;


output.E11 = E11f*Vf + Em*(1-Vf) + (4*Vf*(1-Vf)*(v12f-vm)^2*Gm)/...
    (((1-Vf)*Gm/K23f) + (Vf*Gm/K23m) + 1);

output.v12 = v12f*Vf + vm*(1-Vf) + (Vf*(1-Vf)*(v12f-vm)*(Gm/K23m-Gm/K23f))/...
    (((1-Vf)*Gm/K23f) + (Vf*Gm/K23m) + 1);

output.K23 = K23m + Vf/(1/(K23f-K23m)+(1-Vf)/(K23m+Gm));

output.G12 = Gm*((G12f*(1+Vf) + Gm*(1-Vf))/(G12f*(1-Vf) + Gm*(1+Vf)));

output.G23 = Gm*(1 + (Vf/(Gm/(G23f-Gm) + (K23m+2*Gm)*(1-Vf)/(2*(K23m+Gm)))));

output.E22 = (4*output.G23*output.K23)/...
    (output.K23 + (1 + (4*output.K23*output.v12^2/output.E11))*output.G23);

output.v23 = output.E22/(2*output.G23) - 1;

Cout = [output.E11+4*output.v12^2*output.K23, 2*output.v12*output.K23, 2*output.v12*output.K23, 0, 0, 0;
        2*output.v12*output.K23, output.K23 + output.G23, output.K23 - output.G23, 0, 0, 0;
        2*output.v12*output.K23, output.K23 - output.G23, output.K23 + output.G23, 0, 0, 0;
        0, 0, 0, output.G23, 0, 0;
        0, 0, 0, 0, output.G23, 0;
        0, 0, 0, 0, 0, output.G23]
