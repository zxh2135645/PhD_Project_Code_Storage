function [status] = myplotmod(h_gcf, h_gca, h_plot, h_lgnd, noir)

if( noir )
  set(h_gcf,'Color','k');
  set(h_gca,'Color','k','XColor','w','YColor','w');
  set(h_lgnd,'color','none','TextColor','w');
  set(get(h_gca,'Title'),'Color','w');
end

set(h_plot,'LineWidth',1.5);
set(get(h_gca,'Title'),'FontSize',16);
set(get(h_gca,'Xlabel'),'FontSize',16);
set(get(h_gca,'Ylabel'),'FontSize',16);

status = 1;

return;

