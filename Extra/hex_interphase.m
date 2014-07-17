function hex_interphase(df,vf,it)
% produce a hex packed cell with certain volume fraction and diameter
% df = diameter of fiber
% vf = fiber volume fraction
% ir = normalized interphase thickness (.02 = 2% of the fiber radius)


PLOT_RUC = true;

if vf > .9
    error('FIBER VOLUME FRACTION EXCEEDS MAXIMUM')
    return
end

%x = [0.1387, 0.1387, 0.5547, 0.1387, 0.1387];
dx = [0.277350098112615, 0.138675049056307, 0.138675049056307];
x0  = [0, dx(1), dx(1)+dx(2), sum(dx)];
xfiber = sqrt(pi/4*df^2)*x0;                  % x locations of center fiber
xinter = xfiber(2:end) + it*df/2;
sf = sqrt((pi/2*df^2)/(vf*sqrt(3)));

xleft = sf*(sqrt(3)/2);
xfiber2 = xleft-fliplr(xfiber);
xinter2 = xfiber2(1:end-1)-it*df/2;

x = sort([xfiber,xinter,xfiber2,xinter2]);

yfiber2 = sf*1/2-(xfiber);
yinter2 = yfiber2(2:end)-it*df/2;
y = sort([xfiber,xinter,yfiber2,yinter2]);

for i = 1:size(x,2)-1
    xc(i) = (x(i)+x(i+1))/2;
end
for i = 1:size(y,2)-1
    yc(i) = (y(i)+y(i+1))/2;
end

sm = zeros(size(yc,2),size(xc,2));
for b = 1:size(xc,2)
    for g = 1:size(yc,2)
        % test if center fiber or interphase
        if     xc(b) < xinter(3) && yc(g) < xinter(1) % section 1
            sm(g,b) = 3; % interphase
        elseif xc(b) < xinter(2) && yc(g) < xinter(2) % section 2
            sm(g,b) = 3; % interphase
        elseif xc(b) < xinter(1) && yc(g) < xinter(3) % section 3
            sm(g,b) = 3; % interphase
        end
        
        % test if quarter fiber interphase
        if     xc(b) > xinter2(1) && yc(g) > yinter2(1) % section 1
            sm(g,b) = 3; % interphase
        elseif xc(b) > xinter2(2) && yc(g) > yinter2(2) % section 2
            sm(g,b) = 3; % interphase
        elseif xc(b) > xinter2(3) && yc(g) > yinter2(3) % section 3
            sm(g,b) = 3; % interphase
        end
        
        if xc(b) < xfiber(4) && yc(g) < xfiber(2) % section 1
            sm(g,b) = 1; % fiber
        elseif xc(b) < xfiber(3) && yc(g) < xfiber(3) % section 2
            sm(g,b) = 1; % fiber
        elseif xc(b) < xfiber(2) && yc(g) < xfiber(4) % section 3
            sm(g,b) = 1; % fiber
        end
        
        % test if quarter fiber
        if     xc(b) > xfiber2(1) && yc(g) > yfiber2(2) % section 1
            sm(g,b) = 1; % fiber
        elseif xc(b) > xfiber2(2) && yc(g) > yfiber2(3) % section 2
            sm(g,b) = 1; % fiber
        elseif xc(b) > xfiber2(3) && yc(g) > yfiber2(4) % section 3
            sm(g,b) = 1; % fiber
        end
           
        if sm(g,b) == 0; % matrix
            sm(g,b) = 2; 
        end
    end
end

SM = [rot90(sm,2),flipud(sm);fliplr(sm),sm];

for i = 1:size(x,2)-1
    l(i) = x(i+1)-x(i);
end
for i = 1:size(y,2)-1
    h(i) = y(i+1)-y(i);
end
L = [fliplr(l),l];
H = [h,fliplr(h)];

% vfiber = 0;
% vall   = 0;
% for b = 1:size(L,2)
%     for g = 1:size(H,2)
%         if SM(g,b) == 1;
%             vfiber = vfiber + L(b)*H(g);
%         end
%         vall = vall + L(b)*H(g);
%     end
% end
% VF = vfiber/vall;

if PLOT_RUC % plot the geometry
        x0 = 0;
        figure(1)
        hold on
        % loop over the elements
        for x = 1:size(L,2);
            y0 = 0;
            for y = 1:size(H,2);
                if     SM(y,x) == 1;  color = 'b';
                elseif SM(y,x) == 2; color = 'y';
                elseif SM(y,x) == 3; color = 'r';
                end
                X = [x0, x0, x0+L(x), x0+L(x)];
                Y = [y0, y0+H(y), y0+H(y), y0];
                fill(X,Y,color)
                y0 = y0+H(y);
            end
            x0 = x0+L(x);
        end
        axis equal
        axis tight
        box off
        set(gca,'color','none')
        title(['RUC for V_F = ',num2str(vf,4)],'FontSize',12,'FontWeight','bold')
        xlabel('X (mm)','FontSize',12,'FontWeight','bold')
        ylabel('Y (mm)','FontSize',12,'FontWeight','bold')
end
    
end