x = [0.125, 0.125, 0.5, 0.125, 0.125];

%v = 3*x(3)^2+4*x(4)^2

SM = [0 0 1 0 0;
      0 1 1 1 0;
      1 1 1 1 1;
      0 1 1 1 0;
      0 0 1 0 0];

vfiber = 0;
vall   = 0;
for b = 1:size(x,2)
    for g = 1:size(x,2)
        if SM(g,b) == 1;
            vfiber = vfiber + x(b)*x(g);
        end
        vall = vall + x(b)*x(g);
    end
end

sf = sqrt(1/vfiber);

x_scaled = sf*x

vfiber_scaled = 0;
vall   = 0;
for b = 1:size(x,2)
    for g = 1:size(x,2)
        if SM(g,b) == 1;
            vfiber_scaled = vfiber_scaled + x_scaled(b)*x_scaled(g);
        end
        vall = vall + x_scaled(b)*x_scaled(g);
    end
end
vfiber_scaled