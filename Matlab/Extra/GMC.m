function []=GMC(H,L,grid)
 
%CONSTANTS
%material constants
E11f = 231*10^9;
E22f = 15*10^9;
v12f = .14;
G12f = 24*10^9;
G23f = 5.01*10^9;
Em = 3*10^9;
vm = .36;
v23f = (E22f/(2*G23f))-1;
delf = E22f*v12f^2+E11f*(v23f-1);
Cf = [E11f^2*(v23f-1)/delf,-E11f*E22f*v12f/delf,-E11f*E22f*v12f/delf,0,0,0;...
    -E11f*E22f*v12f/delf, E22f*(E22f*v12f^2-E11f)/((1+v23f)*delf),-E22f*(E22f*v12f^2+E11f*v23f)/((1+v23f)*delf),0,0,0; ...
    -E11f*E22f*v12f/delf,-E22f*(E22f*v12f^2+E11f*v23f)/((1+v23f)*delf), E22f*(E22f*v12f^2-E11f)/((1+v23f)*delf),0,0,0;...
    0,0,0,2*G23f,0,0;...
    0,0,0,0,2*G12f,0;...
    0,0,0,0,0,2*G12f];
c1 = Em*vm/((1+vm)*(1-2*vm));
c2 = Em/(2*(1+vm));
Cm = [c1+2*c2,c1,c1,0,0,0;...
    c1,c1+2*c2,c1,0,0,0;...
    c1,c1,c1+2*c2,0,0,0;...
    0,0,0,c2,0,0;...
    0,0,0,0,c2,0;...
    0,0,0,0,0,c2];
%lengths
    L=[1/sqrt(3),1-1/sqrt(3)];
    %L=[1/(2*sqrt(3)),1/(2*sqrt(3)),(1-1/sqrt(3))/2,(1-1/sqrt(3))/2];
    %L=[1/(4*sqrt(3)),1/(4*sqrt(3)),1/(4*sqrt(3)),1/(4*sqrt(3)),(1-1/sqrt(3))/4,(1-1/sqrt(3))/4,(1-1/sqrt(3))/4,(1-1/sqrt(3))/4];
    %L = [1/(4*sqrt(3)),1/(4*sqrt(3)),1/(4*sqrt(3)),1/(4*sqrt(3))];
    %x=sqrt(1/(6*sqrt(3)));
    %L = [1/(sqrt(3)*2)-x/2,x-1/(sqrt(3)*2),(1/sqrt(3)-x)/2,(1/sqrt(3)-x)/2,x-1/(sqrt(3)*2),1/(sqrt(3)*2)-x/2];
Ng = size(L,2);
l=zeros(1,Ng);
l(1,1)=L(1,1)/2;
for i=2:Ng
    l(1,i)=sum(L(1,1:i-1))+L(1,i)/2;
end
    H=[1/sqrt(3),1-1/sqrt(3)];
    %H=[1/(2*sqrt(3)),1/(2*sqrt(3)),(1-1/sqrt(3))/2,(1-1/sqrt(3))/2];
    %H=[1/(4*sqrt(3)),1/(4*sqrt(3)),1/(4*sqrt(3)),1/(4*sqrt(3)),(1-1/sqrt(3))/4,(1-1/sqrt(3))/4,(1-1/sqrt(3))/4,(1-1/sqrt(3))/4];
    %H = [1/6,1/6,1/6,1/6,1/6,1/6];
    %H = [x/2,1/2-x,x/2,x/2,1/2-x,x/2];
Nb = size(H,2);
h=zeros(1,Nb);
h(1,1)=H(1,1)/2;
for i=2:Nb
    h(1,i)=sum(H(1,1:i-1))+H(1,i)/2;
end
%grid (1 = fiber, 0 = matrix)
     grid = [1,0;...
         0,0];
%     grid = [1,1,0,0;...
%         1,1,0,0;...
%         0,0,0,0;...
%         0,0,0,0];
%     grid = [1,1,1,1,0,0,0,0;...
%         1,1,1,1,0,0,0,0;...
%         1,1,1,1,0,0,0,0;...
%         1,1,1,1,0,0,0,0;...
%         0,0,0,0,0,0,0,0;...
%         0,0,0,0,0,0,0,0;...
%         0,0,0,0,0,0,0,0;...
%         0,0,0,0,0,0,0,0];
%     grid = [1,0,0,1;...
%         0,0,0,0;...
%         0,1,1,0;...
%         0,1,1,0;...
%         0,0,0,0;...
%         1,0,0,1];
%     grid = [1,1,0,0,1,1;... 
%         0,0,0,0,0,0;...
%         0,1,1,1,1,0;...
%         0,1,1,1,1,0;...
%         0,0,0,0,0,0;...
%         1,1,0,0,1,1];
%     
%GIVENS
ex = .005;
ey = 0;
ez = 0;
exy = 0;
exz = 0;
eyz = 0;
eglobal = [ez;ex;ey;exy;eyz;exz];
 
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
A;
%Am
%r22
for g = 1:Ng
    for b = 1:(Nb-1)
        if grid(b,g)==1
            A(count,(g-1)*6+(b-1)*6*Ng+1)=Cf(2,1);
            A(count,(g-1)*6+(b-1)*6*Ng+2)=Cf(2,2);
            A(count,(g-1)*6+(b-1)*6*Ng+3)=Cf(2,3);           
        else
            A(count,(g-1)*6+(b-1)*6*Ng+1)=Cm(2,1);
            A(count,(g-1)*6+(b-1)*6*Ng+2)=Cm(2,2);
            A(count,(g-1)*6+(b-1)*6*Ng+3)=Cm(2,3);
        end
        if grid(b+1,g)==1
            A(count,(g-1)*6+(b)*6*Ng+1)=-Cf(2,1);
            A(count,(g-1)*6+(b)*6*Ng+2)=-Cf(2,2);
            A(count,(g-1)*6+(b)*6*Ng+3)=-Cf(2,3);
        else
            A(count,(g-1)*6+(b)*6*Ng+1)=-Cm(2,1);
            A(count,(g-1)*6+(b)*6*Ng+2)=-Cm(2,2);
            A(count,(g-1)*6+(b)*6*Ng+3)=-Cm(2,3);
        end
        count = count + 1;
    end
end
%r33
for b = 1:Nb
    for g = 1:(Ng-1)
        if grid(b,g)==1
            A(count,(g-1)*6+(b-1)*6*Ng+1)=Cf(3,1);
            A(count,(g-1)*6+(b-1)*6*Ng+2)=Cf(3,2);
            A(count,(g-1)*6+(b-1)*6*Ng+3)=Cf(3,3);    
        else
            A(count,(g-1)*6+(b-1)*6*Ng+1)=Cm(3,1);
            A(count,(g-1)*6+(b-1)*6*Ng+2)=Cm(3,2);
            A(count,(g-1)*6+(b-1)*6*Ng+3)=Cm(3,3);
        end
        if grid(b,g+1)==1
            A(count,(g)*6+(b-1)*6*Ng+1)=-Cf(3,1);
            A(count,(g)*6+(b-1)*6*Ng+2)=-Cf(3,2);
            A(count,(g)*6+(b-1)*6*Ng+3)=-Cf(3,3);   
        else
            A(count,(g)*6+(b-1)*6*Ng+1)=-Cm(3,1);
            A(count,(g)*6+(b-1)*6*Ng+2)=-Cm(3,2);
            A(count,(g)*6+(b-1)*6*Ng+3)=-Cm(3,3);
        end
        count = count + 1;
    end
end
%r12
for g = 1:Ng
    for b = 1:(Nb-1)
        if grid(b,g)==1
            A(count,(g-1)*6+(b-1)*6*Ng+6)=Cf(6,6);            
        else
            A(count,(g-1)*6+(b-1)*6*Ng+6)=Cm(6,6);
        end
        if grid(b+1,g)==1
            A(count,(g-1)*6+(b)*6*Ng+6)=-Cf(6,6);
        else
            A(count,(g-1)*6+(b)*6*Ng+6)=-Cm(6,6);
        end
        count = count + 1;
    end
end
%r13
for b = 1:Nb
    for g = 1:(Ng-1)
        if grid(b,g)==1
            A(count,(g-1)*6+(b-1)*6*Ng+5)=Cf(5,5);            
        else
            A(count,(g-1)*6+(b-1)*6*Ng+5)=Cm(5,5);
        end
        if grid(b,g+1)==1
            A(count,(g)*6+(b-1)*6*Ng+5)=-Cf(5,5);
        else
            A(count,(g)*6+(b-1)*6*Ng+5)=-Cm(5,5);
        end
        count = count + 1;
    end
end
%r23
for g = 1:Ng
    for b = 1:(Nb-1)
        if grid(b,g)==1
            A(count,(g-1)*6+(b-1)*6*Ng+4)=Cf(4,4);            
        else
            A(count,(g-1)*6+(b-1)*6*Ng+4)=Cm(4,4);
        end
        if grid(b+1,g)==1
            A(count,(g-1)*6+(b)*6*Ng+4)=-Cf(4,4);
        else
            A(count,(g-1)*6+(b)*6*Ng+4)=-Cm(4,4);
        end
        count = count + 1;
    end
end
for b = 1:Nb
    for g = 1:(Ng-1)
        if count<=Nb*Ng*6
            if grid(b,g)==1
                A(count,(g-1)*6+(b-1)*6*Ng+4)=Cf(4,4);            
            else
                A(count,(g-1)*6+(b-1)*6*Ng+4)=Cm(4,4);
            end
            if grid(b,g+1)==1
                A(count,(g)*6+(b-1)*6*Ng+4)=-Cf(4,4);
            else
                A(count,(g)*6+(b-1)*6*Ng+4)=-Cm(4,4);
            end
            count = count + 1;
        end
    end
end
%SOLVE
B = A\K;
esub = B*eglobal;
C=zeros(Nb*Ng*6,Nb*Ng*6);
for b=1:Nb
    for g=1:Ng
        i=(b-1)*Ng+g;
        if grid(b,g)==1
            C((6*i-5):(6*i),(6*i-5):(6*i))=Cf;
        else
            C((6*i-5):(6*i),(6*i-5):(6*i))=Cm;
        end
    end
end
rsub = C*esub;
for i=1:size(esub,1)
    if rsub(i,1) < 0.0000001
        rsub(i,1) = 0;
    end
    if esub(i,1) < 0.0000001
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
    s = [r1;r2;r3;r23;r13;r12]
    e = [e1;e2;e3;e23;e13;e12]
    
% figure(1)
% contour(l,h,e2);
% title('E22');
% figure(2)
% contour(l,h,e3);
% title('E33');
% figure(3)
% contour(l,h,e23);
% title('E23');
% figure(4)
% contour(l,h,r2);
% title('S22');
% figure(5)
% contour(l,h,r3);
% title('S33');
% figure(6)
% contour(l,h,r23);
% title('S23');
