function Build_MOC_Input

% function Build_MOC_Input
% This function is used to build MOC input files for the MOC code written
% by Michael Kaplan and Rehan Nawaz with the assistance of Professor
% Anthony Waas
% 
% by Michael Kaplan

% create the filespace

filename = 'FiberCentered_0.5_01.moci';
fid = fopen(filename,'w+');

% create the input header
fprintf(fid,'%% MOC Input: %s\n',filename);
time = datestr(clock);
fprintf(fid,'%% Created: %s\n%%\n',time);

%--------------------------------------------------------------------------
%                      FIBER CONSITUIVE PROPERTIES
%--------------------------------------------------------------------------
fprintf(fid,'*CONSTITUENTS\n');

nmats = 2;  % number of materials present
fprintf(fid,'NMATS=%i\n',nmats);

% MATERIAL 1 (fiber)
mat = 1;
cmod = 1;
%   1 = general elastic
%   2 = isotropic elastic
%   add other models later

if cmod == 1
    E11 = 251e3;   % MPa
    E22 = 40.4e3;    % MPa
    v12 = .321;     
    G12 = 30.7e3;    % MPa
    v23 = 0.256;
    G23 = 1/2*(E22/(v23+1));  % MPa
    
    
    % write input
    fprintf(fid,'MAT=%i,CMOD=%i\n',mat,cmod);
    fprintf(fid,'EA=%8.3e\n',E11);
    fprintf(fid,'ET=%8.3e\n',E22);
    fprintf(fid,'NUA=%0.3f\n',v23);
    fprintf(fid,'NUT=%0.3f\n',v12);
    fprintf(fid,'GA=%8.3e\n',G23);
    fprintf(fid,'GT=%8.3e\n%%\n',G12);

elseif cmod == 2
    E = 3.31e3;   % MPa
    v = 0.318;    

    
    % write input
    fprintf(fid,'MAT=%i,CMOD=%i\n',mat,cmod);
    fprintf(fid,'E=%8.3e\n',E);
    fprintf(fid,'NU=%0.3f\n%%\n',v);
else
    disp('add more material models')
end

% MATERIAL 2 (matrix)
mat = 2;
cmod = 2;
%   1 = general elastic
%   2 = isotropic elastic
%   add other models later

if cmod == 1
    E11 = 251e3;   % MPa
    E22 = 40.4e3;    % MPa
    v12 = .321;     
    G12 = 30.7e3;    % MPa
    v23 = 0.256;
    G23 = 1/2*(E22/(v23+1));  % MPa
    
    
    % write input
    fprintf(fid,'MAT=%i,CMOD=%i\n',mat,cmod);
    fprintf(fid,'EA=%8.3e\n',E11);
    fprintf(fid,'ET=%8.3e\n',E22);
    fprintf(fid,'NUA=%0.3f\n',v23);
    fprintf(fid,'NUT=%0.3f\n',v12);
    fprintf(fid,'GA=%8.3e\n',G23);
    fprintf(fid,'GT=%8.3e\n%%\n',G12);

elseif cmod == 2
    E = 3.31e3;   % MPa
    v = 0.318;     

    
    % write input
    fprintf(fid,'MAT=%i,CMOD=%i\n',mat,cmod);
    fprintf(fid,'E=%8.3e\n',E);
    fprintf(fid,'NU=%0.3f\n%%\n',v);
else
    disp('add more material models')
end

%--------------------------------------------------------------------------
%                            CELL ARCHITECTURE 
%--------------------------------------------------------------------------

fprintf(fid,'*CELL\n');

amod = 2;
%   1 = 4 cell square
%   2 = fiber centered square
%   3 = hex packed rectangle 
%   4 = user defined
%   add other models

if amod == 1
    Vf = 0.5;  % fiber volume fraction
    Df = 5e-3;   % fiber diameter (mm)
    
    fprintf(fid,'AMOD=%i\n',amod);
    fprintf(fid,'VF=%0.3f\n',Vf);
    fprintf(fid,'DF=%8.3e\n%%\n',Df);
    
elseif amod == 2
    Vf = 0.5;  % fiber volume fraction
    Df = 5e-3;   % fiber diameter (m)
    
    fprintf(fid,'AMOD=%i\n',amod);
    fprintf(fid,'VF=%0.3f\n',Vf);
    fprintf(fid,'DF=%8.3e\n%%\n',Df);
    
elseif amod == 3
    Vf = 0.5;  % fiber volume fraction
    Df = 5e-3;   % fiber diameter (m)
    
    fprintf(fid,'AMOD=%i\n',amod);
    fprintf(fid,'VF=%0.3f\n',Vf);
    fprintf(fid,'DF=%8.3e\n%%\n',Df);
    
elseif amod == 4
    DIM = [2,2]; %[H,L]
    
    H = [1,1];  % fiber volume fraction
    L = [1,1];   % fiber diameter (m)
    SM = [1,2,2,2];
    
    fprintf(fid,'AMOD=%i\n',amod);
    
    fprintf(fid,'DIM=');
    for i = 1:size(DIM,2)
        fprintf(fid,'%i,',DIM(i));
    end
    fprintf(fid,'\n');
    
    fprintf(fid,'H=');
    for i = 1:size(H,2)
        fprintf(fid,'%8.3e,',H(i));
    end
    fprintf(fid,'\n');
    
    fprintf(fid,'L=');
    for i = 1:size(L,2)
        fprintf(fid,'%8.3e,',L(i));
    end
    fprintf(fid,'\n');
    
    fprintf(fid,'SM=');
    for i = 1:size(SM,2)
        fprintf(fid,'%i,',SM(i));
    end
    fprintf(fid,'\n%%\n');    
    
end
% NOTE 
% currently all architectures have the same input, but in the future, we
% may have the ability or need to add more inputs for the architectures

%--------------------------------------------------------------------------
%                           LOADING CONDITIONS
%--------------------------------------------------------------------------

fprintf(fid,'*LOADING\n');

lmod = 1;
%   1 = axial strain
%   2 = tangential strain
%   3 = shear strain
%   add other cases

if lmod == 1 % axial strain
    loads = [1e-3]; 
    Nloads = size(loads,2);
    
    fprintf(fid,'LMOD=%i\n',lmod);
    fprintf(fid,'NL=%i\n',Nloads);
    fprintf(fid,'L=');
    for i = 1:Nloads
        fprintf(fid,'%0.5f,',loads(i));
    end
    fprintf(fid,'\n');
    
elseif lmod == 2 % tangential strain
    loads = [1e-3,2e-3];
    Nloads = size(loads,2);
    
    fprintf(fid,'LMOD=%i\n',lmod);
    fprintf(fid,'NL=%i\n',Nloads);
    fprintf(fid,'L=');
    for i = 1:Nloads
        fprintf(fid,'%0.5f,',loads(i));
    end
    fprintf(fid,'\n');
    
elseif lmod == 3 % shear strain
    loads = [1e-3,2e-3];
    Nloads = size(loads,2);
    
    fprintf(fid,'LMOD=%i\n',lmod);
    fprintf(fid,'NL=%i\n',Nloads);
    fprintf(fid,'L=');
    for i = 1:Nloads
        fprintf(fid,'%0.5f,',loads(i));
    end
    fprintf(fid,'\n');
end

% end the input file

fprintf(fid,'%%\n*END');
fclose(fid);



