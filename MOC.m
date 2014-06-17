function [C,s,e] = MOC(input)

input = 'input003';
inputfile = ['Inputs/',input,'.moci'];

[mat,arch,load] = MOC_read(inputfile);


% Materials
nmat = size(mat,1);
Cn = cell(nmat,1);

% Build the C matrix depending on the 
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
amod = arch.amod;

if ismember(amod,[1 2 3])
    error('myApp:argChk', 'ADD FUNCTIONALITY FOR PREDEFINED ARCHITECTURES')
    
elseif amod == 4
    L = arch.l;
    H = arch.h;
    
    c = 1;
    for i = 1:size(L,2)
        for j = 1:size(H,2)           
            SM(i,j) = arch.sm(c);
            c = c+1;
        end
    end
end    

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
    %SOLVE
    %}
    
    %{
    %----------------------------------------------
    %                 try it my self
    %----------------------------------------------
    % strain vector looks like:
    %
    % [e(11) e(12) e(13) e(14) e(21) ... e(NbNg)] transposed
    % loop through gammas then betas
    
    % e22 (Paley eqn. 35)
    for g = 1:Ng
        for b = 1:Nb
            spot = (g-1)*6+(b-1)*6*Ng;
            Ag(count,spot+2) = H(b);
        end
        J(count,2) = sum(H);
        count = count+1;
    end
    
    % e33 (Paley eqn. 36)
    for b = 1:Nb
        for g = 1:Ng
            spot = (g-1)*6+(b-1)*6*Ng;
            Ag(count,spot+3) = L(g);
        end
        J(count,3) = sum(L);
        count = count+1;
    end
    
    % 2*e12 (Paley eqn. 39)
    for g = 1:Ng
        for b = 1:Nb
            spot = (g-1)*6+(b-1)*6*Ng;
            Ag(count,spot+6) = H(b);
        end
        J(count,6) = sum(H);
        count = count+1;
    end
    
    % 2*e13 (Paley eqn. 40)
    for b = 1:Nb
        for g = 1:Ng
            spot = (g-1)*6+(b-1)*6*Ng;
            Ag(count,spot+5) = L(g);
        end
        J(count,5) = sum(L);
        count = count+1;
    end
    
    % 2*e23 (Paley eqn. 32)
    for b = 1:Nb
        for g = 1:Ng
            spot = (g-1)*6+(b-1)*6*Ng;
            Ag(count,spot+4) = (H(b)*L(g))/(sum(H)*sum(L));
        end        
    end
    J(count,4) = 1;
    count = count+1;
    
    % e11 (e11(bg) = e11(avg))
    for b = 1:Nb
        for g = 1:Ng
            spot = (g-1)*6+(b-1)*6*Ng;
            Ag(count,spot+1) = 1;
            J(count,1) = 1;
            count = count+1;
        end        
    end
    
    % build Am
    count = 1;
    
    %r22
    for g = 1:Ng
        for b = 1:(Nb-1)
            spot = (g-1)*6+(b-1)*6*Ng;
            matn = SM(b,g);
                Am(count,(g-1)*6+(b-1)*6*Ng+1)=Cn{matn}(2,1);
                Am(count,(g-1)*6+(b-1)*6*Ng+2)=Cn{matn}(2,2);
                Am(count,(g-1)*6+(b-1)*6*Ng+3)=Cn{matn}(2,3);    
                        
            matm = SM(b+1,g);
                Am(count,(g-1)*6+(b)*6*Ng+1)=-Cn{matm}(2,1);
                Am(count,(g-1)*6+(b)*6*Ng+2)=-Cn{matm}(2,2);
                Am(count,(g-1)*6+(b)*6*Ng+3)=-Cn{matm}(2,3);   

            count = count + 1;
        end
    end
    
    %r33
    for b = 1:Nb
        for g = 1:(Ng-1)
            
            matn = SM(b,g);
                Am(count,(g-1)*6+(b-1)*6*Ng+1)=Cn{matn}(3,1);
                Am(count,(g-1)*6+(b-1)*6*Ng+2)=Cn{matn}(3,2);
                Am(count,(g-1)*6+(b-1)*6*Ng+3)=Cn{matn}(3,3);    
                        
            matm = SM(b,g+1);
                Am(count,(g)*6+(b-1)*6*Ng+1)=-Cn{matm}(3,1);
                Am(count,(g)*6+(b-1)*6*Ng+2)=-Cn{matm}(3,2);
                Am(count,(g)*6+(b-1)*6*Ng+3)=-Cn{matm}(3,3);   

            count = count + 1;
        end
    end
    
    %r12
    for g = 1:Ng
        for b = 1:(Nb-1)
            
            matn = SM(b,g);
                Am(count,(g-1)*6+(b-1)*6*Ng+6)=Cn{matn}(6,6);
                
            matm = SM(b+1,g);
                Am(count,(g-1)*6+(b)*6*Ng+6)=-Cn{matm}(6,6);

            count = count + 1;
        end
    end
    
    %r13
    for b = 1:Nb
        for g = 1:(Ng-1)
                       
            matn = SM(b,g);
                Am(count,(g-1)*6+(b-1)*6*Ng+5)=Cn{matn}(5,5);            

            matm = SM(b,g+1);
                Am(count,(g)*6+(b-1)*6*Ng+5)=-Cn{matm}(5,5);

            count = count + 1;
        end
    end
    
    %r23
    for g = 1:Ng
        for b = 1:(Nb-1)
            
            matn = SM(b,g);
                Am(count,(g-1)*6+(b-1)*6*Ng+4)=Cn{matn}(4,4);            

            matm = SM(b+1,g);
                Am(count,(g-1)*6+(b)*6*Ng+4)=-Cn{matm}(4,4);

            count = count + 1;
        end
    end
    
    %r32
%     for b = 1:Nb
%         for g = 1:(Ng-1)
%             if count<=Nb*Ng*6
%                 
%                 matn = SM(b,g);
%                     Am(count,(g-1)*6+(b-1)*6*Ng+4)=Cn{matn}(4,4);            
% 
%                 matm = SM(b,g+1);
%                     Am(count,(g)*6+(b-1)*6*Ng+4)=-Cn{matm}(4,4);
% 
%                 count = count + 1;
%             end
%         end
%     end
    
    %SOLVE
    A = [Am;Ag]
    K = [zeros(5*(Nb*Ng)-2*(Nb+Ng)-1,6);J];
    
    %----------------------------------------------
    %}
    
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
    
   
    for i=1:size(esub,1)
        if abs(rsub(i,1)) < 0.00001
            rsub(i,1) = 0;
        end
        if abs(esub(i,1)) < 0.00001
            esub(i,1) = 0;
        end
    end
    
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
    s = [r1;r2;r3;r23;r13;r12];
    e = [e1;e2;e3;e23;e13;e12];
    
    s_avg = [mean(mean(r1));mean(mean(r2));mean(mean(r3));mean(mean(r23));mean(mean(r13));mean(mean(r12))]
    e_avg = [mean(mean(e1));mean(mean(e2));mean(mean(e3));mean(mean(e23));mean(mean(e13));mean(mean(e12))]
    
    % output to file
    
    
end % loading
end

