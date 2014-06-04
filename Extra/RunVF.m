function RunVF
close all

N = 10;% number of volume fractions to calculate


PLOT_RUC = true;
PLOT_XY  = false;

LoadCase = 2;
%  1 = Axial Tension
%  2 = Transverse Tension
%  3 = Axial Shear

% setup the filespace
mkdir('INPUTS');
mkdir('OUTPUTS');
sleg = cell(N,1);

Vf = linspace(.6239,.7160,N);

H = [0,0.0006145,0.0006145,0.002458,0.0006145,0.0006145,0];
L = H; % square symmetric element

SM = [2,2,2,2,2,2,2;
      2,2,2,1,2,2,2;
      2,2,1,1,1,2,2;
      2,1,1,1,1,1,2;
      2,2,1,1,1,2,2;
      2,2,2,1,2,2,2;
      2,2,2,2,2,2,2];  

fiber  = 1;
matrix = 2;

% initialize the volumes
Vfiber  = 0;
Vmatrix = 0;

% calculate the inner square length
dx = sum(L);
dy = sum(H);

% loop over elements
for x = 1:size(L,2);
    for y = 1:size(H,2);
        if     SM(y,x) == fiber % calculate the volume of the fiber
            Vfiber = Vfiber + L(x)*H(y);
        elseif SM(y,x) == matrix % calculate the volume of the surrounding matrix
            Vmatrix = Vmatrix + L(x)*H(y);
        end
    end
end
Vfiber/(Vfiber+Vmatrix);

for i = 1:size(Vf,2)
    
    % set files
    vf = Vf(i);
    
    if     LoadCase == 1; sLoad = 'AT';
    elseif LoadCase == 2; sLoad = 'TT';
    elseif LoadCase == 3; sLoad = 'AS';
    end;
    
    input  = sprintf('INPUTS/Example_8D_%s_vf(%0.2f)',sLoad, vf);
    output = sprintf('OUTPUTS/Example_8D_%s_vf(%0.2f)',sLoad, vf);
    
    % Calculate the needed extra matrix volume
    syms Vm_Needed
    Vm_Needed = solve(vf == Vfiber/(Vfiber+Vmatrix+Vm_Needed),Vm_Needed);
    Vm_Needed = eval(Vm_Needed);
    
    % calculate the needed height and lengths of the surrounding matrix
    syms T
    T = solve(Vm_Needed==4*T^2+4*T*dx,T);
    T = max(eval(T)); % the thickness of the surrounding membrane
    
    % Complete the H and L vectors
    H(1) = T;
    H(end) = T;
    L(1) = T;
    L(end) = T;
    
    if PLOT_RUC % plot the geometry
        x0 = 0;
        figure(i)
        hold on
        % loop over the elements
        for x = 1:size(L,2);
            y0 = 0;
            for y = 1:size(H,2);
                if     SM(y,x) == fiber;  color = 'b';
                elseif SM(y,x) == matrix; color = 'y';
                end
                X = [x0, x0, x0+L(x), x0+L(x)];
                Y = [y0, y0+H(y), y0+H(y), y0];
                fill(X,Y,color)
                y0 = y0+H(y);
            end
            x0 = x0+L(x);
        end
        axis equal
        axis square
        box off
        set(gca,'color','none')
        title(['RUC for V_F = ',num2str(vf,4)],'FontSize',12,'FontWeight','bold')
        xlabel('X (mm)','FontSize',12,'FontWeight','bold')
        ylabel('Y (mm)','FontSize',12,'FontWeight','bold')
    end
   
end
end



function InputBuilder(input,LoadCase,vf,H,L)

fID = fopen([input,'.mac'],'w+'); % open the input file that we are going to write

    if     LoadCase == 1; sLoad = 'AT';
    elseif LoadCase == 2; sLoad = 'TT';
    elseif LoadCase == 3; sLoad = 'AS';
    end;

fprintf(fID,'MAC/GMC Example 8D (volume fraction = %5.2f)\n',vf);

fprintf(fID,'*CONSTITUENTS\n');
fprintf(fID,'  NMATS=2\n');
fprintf(fID,'# -- IM7 fiber\n');
fprintf(fID,'  M=1 CMOD=6 MATID=U MATDB=1 &\n   EL=275.E3,16.E3,0.31,0.31,105.E3,0.,0.\n');
fprintf(fID,'# --  Epoxy matrix\n');
fprintf(fID,'  M=2 CMOD=6 MATID=U MATDB=1 &\n   EL=3.5E3,3.5E3,0.35,0.35,1.29E3,0.,0.\n');

fprintf(fID,'*RUC\n');
fprintf(fID,'   MOD=2 ARCHID=99\n');
fprintf(fID,'   NB=7 NG=7\n');

% write the H and L vectors
fprintf(fID,'   H=');
fprintf(fID,'%0.7f,',H([1:end-1]));
fprintf(fID,'%0.7f\n',H(end));

fprintf(fID,'   L=');
fprintf(fID,'%0.7f,',L([1:end-1]));
fprintf(fID,'%0.7f\n',L(end));

fprintf(fID,'   SM=2,2,2,2,2,2,2\n');
fprintf(fID,'   SM=2,2,2,1,2,2,2\n');
fprintf(fID,'   SM=2,2,1,1,1,2,2\n');
fprintf(fID,'   SM=2,1,1,1,1,1,2\n');
fprintf(fID,'   SM=2,2,1,1,1,2,2\n');
fprintf(fID,'   SM=2,2,2,1,2,2,2\n');
fprintf(fID,'   SM=2,2,2,2,2,2,2\n');

% set the load case
fprintf(fID,'*MECH\n');
if LoadCase == 1
    fprintf(fID,'  LOP=1\n');
    fprintf(fID,'  NPT=2 TI=0.,1. MAG=0.,0.02 MODE=1\n');
    fprintf(fID,'  NPT=2 TI=0.,1. MAG=0.,-0.015 MODE=1\n');
elseif LoadCase == 2
    fprintf(fID,'  LOP=2\n');
    fprintf(fID,'  NPT=2 TI=0.,1. MAG=0.,0.05 MODE=1\n');
    fprintf(fID,'  NPT=2 TI=0.,1. MAG=0.,-0.04 MODE=1\n');
elseif LoadCase == 3
    fprintf(fID,'  LOP=6\n');
    fprintf(fID,'  NPT=2 TI=0.,1. MAG=0.,0.1 MODE=1\n');
end

fprintf(fID,'*SOLVER\n');
fprintf(fID,'  METHOD=1 NPT=2 TI=0.,1. STP=0.004 ITMAX=500 ERR=1.E-5\n');
fprintf(fID,'  NLEG=1 NINTEG=1\n');

fprintf(fID,'*XYPLOT\n');
fprintf(fID,'  FREQ=1\n');
fprintf(fID,'  MACRO=1\n');

name = sprintf('Example_8D_%s_vf(%0.2f)',sLoad, vf);
if LoadCase == 1
    fprintf(fID,'   name=%s X=1 Y=7\n',name);
elseif LoadCase == 2
    fprintf(fID,'   name=%s X=2 Y=8\n',name);
elseif LoadCase == 3
    fprintf(fID,'   name=%s X=6 Y=12\n',name);
end
fprintf(fID,'  MICRO=0\n');

fprintf(fID,'*FAILURE_SUBCELL\n');
fprintf(fID,'  DMAX=0.99999\n');
fprintf(fID,'  NMAT=2\n');
fprintf(fID,'  MAT=1 NCRIT=1\n');
fprintf(fID,'   CRIT=1 X11=4600. COMPR=DIF XC11=2100. ACTION=1\n');
fprintf(fID,'  MAT=2 NCRIT=1\n');
fprintf(fID,'   CRIT=7 X11=0.0125 X22=0.0125 X33=0.0125 X23=0.0103 X13=0.0103 X12=0.0103 COMPR=DIF &\n');
fprintf(fID,'       XC11=0.0307 XC22=0.0307 XC33=0.0307 ACTION=4\n');
fprintf(fID,'    CMOD=4 PDMOD=3\n');
fprintf(fID,'     K11=2.5,4.5 K22=2.5,4.5 K33=2.5,4.5\n');
fprintf(fID,'     B11=1. B22=1. B33=1. B44=1. B55=1. B66=1. &\n');
fprintf(fID,'        B42=0.5 B43=0.5 B51=0.5 B53=0.5 B61=0.5 &\n');
fprintf(fID,'        B62=0.5\n');
fprintf(fID,'     MMLAW=1 GIC=1.8 GIIC=5.2 GIIIC=5.2\n');
fprintf(fID,'     LEN=0.\n');
fprintf(fID,'    CMODC=4 PDMODC=3\n');
fprintf(fID,'     KC11=5.0,0.9 KC22=5.0,0.9 KC33=5.0,0.9\n');
fprintf(fID,'	   BC11=1. BC22=1. BC33=1.\n');
fprintf(fID,'     MMLAWC=4 WSC=5.55E-6\n');

fprintf(fID,'*PRINT\n');
fprintf(fID,'  NPL=6\n');

fprintf(fID,'*END');

fclose(fID);

end
