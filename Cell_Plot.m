function Cell_Plot(value, L, H)
% Cell_Plot(value, L, H)
% Plot the value over the RUC. Value can be any stress or strain

X = zeros(size(L,2)+1,1);
Y = zeros(size(H,2)+1,1);

% So it is displayed correctly must flip about horizontal
H = fliplr(H);
value = flipud(value);

X(1) = 0;
for i = 1:max(size(L))
    X(i+1) = X(i)+L(i);
end
Y(1) = 0;
for j = 1:max(size(H))
    Y(j+1) = Y(j)+H(j);
end

V = zeros(size(value)+1);
V(1:end-1,1:end-1) = value;

H

figure(1)
pcolor(X,Y,V)
colorbar

end

