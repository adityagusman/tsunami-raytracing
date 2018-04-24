% main raytracing example
% by Aditya Riadi Gusman
% Earthquake Research Institute, the University of Tokyo
% March 21, 2017
clear all
close all

% parameters
h=10;                        % time step sec
mt=7200;                    % maximum time sec
zetaol = 0:1:360;         % initial directions
xo=  141.3871;              % initial location longitude
yo=   37.3737;              % initial location latitude

save rt_parameters.mat

% load bathymetry
load grid_a.mat;
load xya.mat

%% ray tracing program
[gph, gth, gz] = raytracing_sp(xa,ya,grid_a,h,mt,xo,yo,zetaol);

%%
figure
contour(xa,ya,grid_a');
hold on
contour(xa,ya,grid_a',[0 0],'b','linewidth',1);


for iz=1:length(zetaol)
    plot(gph(iz,:),gth(iz,:),'k');
end
plot(xo,yo,'p','markerfacecolor','r','markeredgecolor','k',...
    'markersize',14)
axis equal
hc=colorbar;
ylabel(hc,'Depth, m')
saveas(gcf,'raytracing_sp_c.fig')
print(gcf,'-djpeg','-r300',['raytracing_sp_c.jpg'])