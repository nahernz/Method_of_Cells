function [C_star,s_avg,e_avg] = MOC(input)

close all

%input = '2x2_0.5_01_AS';
%input = '7x7_0.5_01_AS';
%input = 'Hex_0.5_01_AS';
input = 'HexInter_0.5_01_AS';
%input = 'SquareInter_0.5_01_AS';
%input = 'SquareInter_0.5_01_AS_(Test)';

%input = '2x2_0.5_01_TS';

inputfile = ['Inputs/',input,'.moci'];

[mat,arch,load,matlab] = MOC_read(inputfile);


% Materials
nmat = size(mat,1);
Cn = cell(nmat,1);

% Build the C matrix depending on the material
for i = 1:nmat
    cmod = mat{i}.cmod;
    if cmod == 1
        E11 = mat{i}.E11;
        E22 = mat{i}.E22;
        v23 = mat{i}.V23;
        v12 = mat{i}.V12;
        G23 = mat{i}.G23;
        G12 = mat{i}.G12;
        
        
        del = E22*v12^2+E11*(v23-1);
        Cn{i} = [E11^2*(v23-1)/del,-E11*E22*v12/del,-E11*E22*v12/del,0,0,0;
                -E11*E22*v12/del, E22*(E22*v12^2-E11)/((1+v23)*del),-E22*(E22*v12^2+E11*v23)/((1+v23)*del),0,0,0; 
                -E11*E22*v12/del,-E22*(E22*v12^2+E11*v23)/((1+v23)*del), E22*(E22*v12^2-E11)/((1+v23)*del),0,0,0;
                0,0,0,2*G23,0,0;
                0,0,0,0,2*G12,0;
                0,0,0,0,0,2*G12];
            
    elseif cmod == 2
        E = mat{i}.E;
        v = mat{i}.V;
        
        c1 = E*v/((1+v)*(1-2*v));
        c2 = E/(2*(1+v));
        Cn{i} = [c1+2*c2,c1,c1,0,0,0;
                 c1,c1+2*c2,c1,0,0,0;
                 c1,c1,c1+2*c2,0,0,0;
                 0,0,0,c2,0,0;
                 0,0,0,0,c2,0;
                 0,0,0,0,0,c2];
    end
end
    

% Cell Architecture
%amod = arch.amod;
   
    L = arch.l;
    H = arch.h;

    
    SM = arch.sm;
   

Ng = size(L,2);
l=zeros(1,Ng);
l(1,1)=L(1,1)/2;
for i=2:Ng
    l(1,i)=sum(L(1,1:i-1))+L(1,i)/2;
end

Nb = size(H,2);
h=zeros(1,Nb);
h(1,1)=H(1,1)/2;
for i=2:Nb
    h(1,i)=sum(H(1,1:i-1))+H(1,i)/2;
end


%GIVENS
ex = 0;
ey = 0;
ez = 0;
exy = 0;
eyz = 0;
exz = 0;

lmod = load.lmod;
nloads = load.nl;
loads = load.l;

% create loop for each load, assign strains for each load, compute outputs

for nl = 1:nloads
    if lmod == 1
        ex = loads(nl);
        
    elseif lmod == 2
        ey = loads(nl);
        
    elseif lmod == 3
        eyz = loads(nl);
        
    elseif lmod == 4
        error('myApp:argChk', 'ADD FUNCTIONALITY FOR GENERAL LOADING')
        
    end
    
    eglobal = [ex;ey;ez;eyz;exz;exy];
 
    %CREATE MATRICES

    A = zeros(Ng*Nb*6,Ng*Nb*6);
    K = zeros(Ng*Nb*6,6);
    count = 1;
    %Ag and J
   
    %e33
    for b = 1:Nb
        for g = 1:Ng
            A(count,(g-1)*6+(b-1)*6*Ng+3)=L(1,g);
        end
        K(count,3) = sum(L);
        count = count + 1;
    end
    %e22
    for g = 1:Ng
        for b = 1:Nb
            A(count,(g-1)*6+(b-1)*6*Ng+2)=H(1,b);
        end
        K(count,2) = sum(H);
        count = count + 1;
    end
    %e11
    for i = 1:Nb*Ng
        A(count,(i-1)*6+1)=1;
        K(count,1) = 1;
        count = count + 1;
    end
    %e13
    for b = 1:Nb
        for g = 1:Ng
            A(count,(g-1)*6+(b-1)*6*Ng+5)=L(1,g);
        end
        K(count,5) = sum(L);
        count = count + 1;
    end
    %e12
    for g = 1:Ng
        for b = 1:Nb
            A(count,(g-1)*6+(b-1)*6*Ng+6)=H(1,b);
        end
        K(count,6) = sum(H);
        count = count + 1;
    end
    %e23
    for g = 1:Ng
        for b = 1:Nb
            A(count,(g-1)*6+(b-1)*6*Ng+4)=H(1,b)*L(1,g);
        end
    end
    K(count,4) = sum(H)*sum(L);
    count = count + 1;

    %Am
    %r22
    for g = 1:Ng
        for b = 1:(Nb-1)
            
            matn = SM(b,g);
                A(count,(g-1)*6+(b-1)*6*Ng+1)=Cn{matn}(2,1);
                A(count,(g-1)*6+(b-1)*6*Ng+2)=Cn{matn}(2,2);
                A(count,(g-1)*6+(b-1)*6*Ng+3)=Cn{matn}(2,3);    
                        
            matm = SM(b+1,g);
                A(count,(g-1)*6+(b)*6*Ng+1)=-Cn{matm}(2,1);
                A(count,(g-1)*6+(b)*6*Ng+2)=-Cn{matm}(2,2);
                A(count,(g-1)*6+(b)*6*Ng+3)=-Cn{matm}(2,3);   

            count = count + 1;
        end
    end
    %r33
    for b = 1:Nb
        for g = 1:(Ng-1)
            
            matn = SM(b,g);
                A(count,(g-1)*6+(b-1)*6*Ng+1)=Cn{matn}(3,1);
                A(count,(g-1)*6+(b-1)*6*Ng+2)=Cn{matn}(3,2);
                A(count,(g-1)*6+(b-1)*6*Ng+3)=Cn{matn}(3,3);    
                        
            matm = SM(b,g+1);
                A(count,(g)*6+(b-1)*6*Ng+1)=-Cn{matm}(3,1);
                A(count,(g)*6+(b-1)*6*Ng+2)=-Cn{matm}(3,2);
                A(count,(g)*6+(b-1)*6*Ng+3)=-Cn{matm}(3,3);   

            count = count + 1;
        end
    end
    %r12
    for g = 1:Ng
        for b = 1:(Nb-1)
            
            matn = SM(b,g);
                A(count,(g-1)*6+(b-1)*6*Ng+6)=Cn{matn}(6,6);
                
            matm = SM(b+1,g);
                A(count,(g-1)*6+(b)*6*Ng+6)=-Cn{matm}(6,6);

            count = count + 1;
        end
    end
    %r13
    for b = 1:Nb
        for g = 1:(Ng-1)
                       
            matn = SM(b,g);
                A(count,(g-1)*6+(b-1)*6*Ng+5)=Cn{matn}(5,5);            

            matm = SM(b,g+1);
                A(count,(g)*6+(b-1)*6*Ng+5)=-Cn{matm}(5,5);

            count = count + 1;
        end
    end
    %r23
    for g = 1:Ng
        for b = 1:(Nb-1)
            
            matn = SM(b,g);
                A(count,(g-1)*6+(b-1)*6*Ng+4)=Cn{matn}(4,4);            

            matm = SM(b+1,g);
                A(count,(g-1)*6+(b)*6*Ng+4)=-Cn{matm}(4,4);

            count = count + 1;
        end
    end
    for b = 1:Nb
        for g = 1:(Ng-1)
            if count<=Nb*Ng*6
                
                matn = SM(b,g);
                    A(count,(g-1)*6+(b-1)*6*Ng+4)=Cn{matn}(4,4);            

                matm = SM(b,g+1);
                    A(count,(g)*6+(b-1)*6*Ng+4)=-Cn{matm}(4,4);

                count = count + 1;
            end
        end
    end

    B = A\K;
    
    esub = B*eglobal;
    C=zeros(Nb*Ng*6,Nb*Ng*6);
    for b=1:Nb
        for g=1:Ng
            i=(b-1)*Ng+g;
            matn = SM(b,g);
                C((6*i-5):(6*i),(6*i-5):(6*i))=Cn{matn};
        end
    end
    rsub = C*esub;
    
  %{ 
    for i=1:size(esub,1)
        if abs(rsub(i,1)) < 0.00001
            rsub(i,1) = 0;
        end
        if abs(esub(i,1)) < 0.00001
            esub(i,1) = 0;
        end
    end
    %}
    
    e1=zeros(Nb,Ng);
    e2=zeros(Nb,Ng);
    e3=zeros(Nb,Ng);
    e23=zeros(Nb,Ng);
    e13=zeros(Nb,Ng);
    e12=zeros(Nb,Ng);
    for b=1:Nb
        for g=1:Ng
            e1(b,g)  = esub((g-1)*6+(b-1)*6*Ng+1,1);
            e2(b,g)  = esub((g-1)*6+(b-1)*6*Ng+2,1);
            e3(b,g)  = esub((g-1)*6+(b-1)*6*Ng+3,1);
            e23(b,g) = esub((g-1)*6+(b-1)*6*Ng+4,1);
            e13(b,g) = esub((g-1)*6+(b-1)*6*Ng+5,1);
            e12(b,g) = esub((g-1)*6+(b-1)*6*Ng+6,1);
        end
    end
    r1=zeros(Nb,Ng);
    r2=zeros(Nb,Ng);
    r3=zeros(Nb,Ng);
    r23=zeros(Nb,Ng);
    r13=zeros(Nb,Ng);
    r12=zeros(Nb,Ng);
    for b=1:Nb
        for g=1:Ng
            r1(b,g)  = rsub((g-1)*6+(b-1)*6*Ng+1,1);
            r2(b,g)  = rsub((g-1)*6+(b-1)*6*Ng+2,1);
            r3(b,g)  = rsub((g-1)*6+(b-1)*6*Ng+3,1);
            r23(b,g) = rsub((g-1)*6+(b-1)*6*Ng+4,1);
            r13(b,g) = rsub((g-1)*6+(b-1)*6*Ng+5,1);
            r12(b,g) = rsub((g-1)*6+(b-1)*6*Ng+6,1);
        end
    end
%     s = [r1;r2;r3;r23;r13;r12];
%     e = [e1;e2;e3;e23;e13;e12];
    
    s_avg  = zeros(6,1);
    e_avg  = zeros(6,1);
    C_star = zeros(6,6);
    V      = zeros(1,1);
    
    for b=1:Nb
        for g=1:Ng
            sbg   = [r1(b,g);r2(b,g);r3(b,g);r23(b,g);r13(b,g);r12(b,g)];  % subcell stress
            ebg   = [e1(b,g);e2(b,g);e3(b,g);e23(b,g);e13(b,g);e12(b,g)];  % subcell strain
            Vbg   = H(b)*L(g);                                             % subcell volume
            s_avg = s_avg + sbg*Vbg;
            e_avg = e_avg + ebg*Vbg;
            V     = V + Vbg;
        end
    end
    
    i = 1;
    for b=1:Nb
        for g=1:Ng
            matn = SM(b,g);

            C_star = C_star + 1/(sum(H)*sum(L))*H(b)*L(g)*Cn{matn}*B((i-1)*6+1:(i)*6,1:6);
            i = i + 1;
            
        end
    end
    
    s_avg = s_avg/V;
    e_avg = e_avg/V;
    
    
    % plot and output matlab data
    
    StressTitles = {'Axial Stress (\sigma_1_1)','Transverse Stress (\sigma_2_2)',...
        'Transverse Stress (\sigma_3_3)', 'Axial Shear Stress (\sigma_2_3)',...
        'Transverse Shear Stress (\sigma_1_3)','Transverse Shear Stress (\sigma_1_2)'};
    
    StrainTitles = {'Axial Strain (\epsilon_1_1)','Transverse Strain (\epsilon_2_2)',...
        'Transverse Strain (\epsilon_3_3)', 'Axial Shear Strain (\epsilon_2_3)',...
        'Transverse Shear Strain (\epsilon_1_3)','Transverse Shear Strain (\epsilon_1_2)'};
    
%     StiffTitles = {'Axial Stiffness (E_1_1)','Transverse Stiffness (E_2_2)',...
%         'Transverse Stiffness (E_3_3)','Axial Shear Stiffness (G_2_3)',...
%         'Transverse Shear Stiffness (G_1_3)','Transverse Shear Stiffness (G_1_2)'};

    if isfield(matlab,'ns')
        for i = 1:matlab.ns

            if     matlab.s(i)== 1; P_var = r1;
            elseif matlab.s(i)== 2; P_var = r2;
            elseif matlab.s(i)== 3; P_var = r3;
            elseif matlab.s(i)== 4; P_var = r23;
            elseif matlab.s(i)== 5; P_var = r13;    
            elseif matlab.s(i)== 6; P_var = r12;
            end

            Cell_Plot(P_var,L,H)
            title(StressTitles{matlab.s(i)})
        end
    end
    
    if isfield(matlab,'ne')    
        for i = 1:matlab.ne

            if     matlab.e(i)== 1; P_var = e1;
            elseif matlab.e(i)== 2; P_var = e2;
            elseif matlab.e(i)== 3; P_var = e3;
            elseif matlab.e(i)== 4; P_var = e23;
            elseif matlab.e(i)== 5; P_var = e13;    
            elseif matlab.e(i)== 6; P_var = e12;
            end

            Cell_Plot(P_var,L,H)
            title(StrainTitles{matlab.e(i)})
        end
    end
    
    % Values needed to compare to CCM model:
    E11 = C_star(1,1) - 2*C_star(1,2)^2/(C_star(2,2) + C_star(2,3))
    v12 = C_star(1,2)/(C_star(2,2) + C_star(2,3))
    G12 = C_star(6,6)
    K23 = 1/2*(C_star(2,2) + C_star(2,3))
    G23 = 1/2*(C_star(2,2) - C_star(2,3))
    E22 = 4*G23*K23/(K23 + (1 + 4*K23*v12^2/E11)*G23)
    v23 = E22/(2*G23) - 1

   %E112 = s_avg(1)/e_avg(1);
   %E222 = s_avg(2)/e_avg(2)
    % output to file
    
    
end % loading
end

