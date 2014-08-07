function ArchPlot(input)

close all

%input = '2x2_0.5_01_AS';
input = '7x7_0.5_01_AS';
%input = 'Hex_0.5_01_AS';
%input = 'HexInter_0.5_01_AS';
%input = 'HexInter_0.5_02_AS';
%input = 'SquareInter_0.5_01_AS';
%input = 'SquareInter_0.5_01_AS_(Test)';

%input = '2x2_0.5_01_TS';

inputfile = ['Inputs/',input,'.moci'];

[~,arch,~,~] = MOC_read(inputfile);

L  = arch.l
H  = fliplr(arch.h)
SM = flipud(arch.sm)

x0 = 0;
figure(1)
hold on
% loop over the elements
for x = 1:size(L,2);
    y0 = 0;
    for y = 1:size(H,2);
        if     SM(y,x) == 1; color = [.75 0 0];
        elseif SM(y,x) == 2; color = [.2 .2 1];
        elseif SM(y,x) == 3; color = [1 1 1];
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
%title(['RUC for ',input],'FontSize',12,'FontWeight','bold')
xlabel('X (mm)','FontSize',12,'FontWeight','bold')
ylabel('Y (mm)','FontSize',12,'FontWeight','bold')