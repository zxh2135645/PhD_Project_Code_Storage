function PlotBullsEye(SegPixPerc)


c = createBullseye([0 1 1 0; 1 2 4 45; 2 3 6 0; 3 4 6 0]);
set(c,'Color','w','LineWidth',1)

% Filling the bullseye, vector by vector
fillBullseye(flip(SegPixPerc(1:6)),3,4,60,420);
fillBullseye(flip(SegPixPerc(7:12)),2,3,60,420);
fillBullseye(flip(SegPixPerc(13:16)),1,2,45,405);
axis image;

uistack(c,'top');
set(gca,'Color',[50/255 50/255 50/255]);
%caxis([0, 100]);
colormap parula;
clr = colorbar;
%clr.Location = 'east';
clr.Color = [1 1 1];
%clr.TickLabels = {'0%', '100%'};
clr.Ticks = [0, 100];
%clr.Limits = [0, 100];
set(gca, 'FontSize', 12)


color_spec = zeros(16, 3);
for i = 1:16
   if SegPixPerc(i) > (0.8*range(SegPixPerc) + min(SegPixPerc))
       color_spec(i,:) = [0 0 0];
   else
       color_spec(i,:) = [1 1 1];
   end
end

h1 = text(0, 3.5, num2str(round(SegPixPerc(6))), 'HorizontalAlignment', 'center', 'Color', color_spec(6,:), 'FontSize', 16);
h7 = text(0, 2.5, num2str(round(SegPixPerc(12))), 'HorizontalAlignment', 'center', 'Color', color_spec(12,:), 'FontSize', 16);
h13 = text(0, 1.5, num2str(round(SegPixPerc(16))), 'HorizontalAlignment', 'center', 'Color', color_spec(16,:), 'FontSize', 16);
h15 = text(0, -1.5, num2str(round(SegPixPerc(14))), 'HorizontalAlignment', 'center', 'Color', color_spec(14,:), 'FontSize', 16);
h10 = text(0, -2.5, num2str(round(SegPixPerc(9))), 'HorizontalAlignment', 'center', 'Color', color_spec(9,:), 'FontSize', 16);
h4 = text(0, -3.5, num2str(round(SegPixPerc(3))), 'HorizontalAlignment', 'center', 'Color', color_spec(3,:), 'FontSize', 16);
h16 = text(1.5, 0, num2str(round(SegPixPerc(13))), 'HorizontalAlignment', 'center', 'Color', color_spec(13,:), 'FontSize', 16);
h14 = text(-1.5, 0, num2str(round(SegPixPerc(15))), 'HorizontalAlignment', 'center', 'Color', color_spec(15,:), 'FontSize', 16);
h12 = text(2.1, 1.4, num2str(round(SegPixPerc(7))), 'HorizontalAlignment', 'center', 'Color', color_spec(7,:), 'FontSize', 16);
h6 = text(2.9, 2.0, num2str(round(SegPixPerc(1))), 'HorizontalAlignment', 'center', 'Color', color_spec(1,:), 'FontSize', 16);
h8 = text(-2.1, 1.4, num2str(round(SegPixPerc(11))), 'HorizontalAlignment', 'center', 'Color', color_spec(11,:), 'FontSize', 16);
h2 = text(-2.9, 2.0, num2str(round(SegPixPerc(5))), 'HorizontalAlignment', 'center', 'Color', color_spec(5,:), 'FontSize', 16);
h11 = text(2.1, -1.4, num2str(round(SegPixPerc(8))), 'HorizontalAlignment', 'center', 'Color', color_spec(8,:), 'FontSize', 16);
h5 = text(2.9, -2.0, num2str(round(SegPixPerc(2))), 'HorizontalAlignment', 'center', 'Color', color_spec(2,:), 'FontSize', 16);
h9 = text(-2.1, -1.4, num2str(round(SegPixPerc(10))), 'HorizontalAlignment', 'center', 'Color', color_spec(10,:), 'FontSize', 16);
h3 = text(-2.9, -2.0, num2str(round(SegPixPerc(4))), 'HorizontalAlignment', 'center', 'Color', color_spec(4,:), 'FontSize', 16);

end