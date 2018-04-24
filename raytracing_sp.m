function [gph, gth, gz] = raytracing_sp(xa,ya,d,h,mt,xo,yo,zetaol)
% xa, ya    : longitude and latitute arrays
% d         : Bathymetry data
% h         : time step
% mt        : maximum simulated time
% xo,yo     : location of point source in lon and lat
% zetaol    : list of initial directions

% by Aditya Riadi Gusman
% Earthquake Research Institute, the University of Tokyo
% March 21, 2017

% parameter
raddeg = acos(-1)/(180); % degree to radian

th=90-ya;           % colatitude
ds=abs(th(2)-th(1));   % grid spacing in degree
dth=ds*raddeg;      % grid spacing rad
dph=dth;            % grid spacing rad
[si, sj]= size(d);  % y size, x size

%ixo=round((xo-min(xa))/(ds));
%iyo=round((max(ya)-yo)/(ds));
cx=abs(xa-xo);
cxi=find(cx==min(cx));
ixo=cxi(1)
cy=abs(ya-yo);
cyi=find(cy==min(cy));
iyo=cyi(1)

d(d<0)=0;
n = 1./sqrt(9.8 * d); % slowness s/m
n(d<=0)=1;

% slowness gradient
dndph0=diff(n)/dph;
dndth0=diff(n,1,2)/dth;
dndph=dndph0(1:si-1,:);
dndth=dndth0(:,1:sj-1);

% time array
t = 0:h:mt;        % Calculates upto mt

% initialized gx and gy
gx = zeros(length(zetaol),length(t));
gy = zeros(length(zetaol),length(t));
gph = zeros(length(zetaol),length(t));
gth = zeros(length(zetaol),length(t));
gz = zeros(length(zetaol),length(t));

%% raytracing calculation

for iz = 1:length(zetaol)
    zetao = zetaol(iz);
    [ph,th,z,ixa,iya] = rt_rungekutta_sp(t,h,dph,dth,n,dndph,dndth,...
        zetao,xo,yo,ixo,iyo,si,sj,d);
    
    gx(iz,1:length(ixa)) = ixa;
    gy(iz,1:length(iya)) = iya;
    
    gph(iz,1:length(ph)) = ph/raddeg;      % longitude
    gth(iz,1:length(th)) = 90-(th/raddeg); % latitude
    gz(iz,1:length(z)) = z/raddeg;        % directions in degree
end

gx(gx==0) = NaN;
gy(gy==0) = NaN;
gph(gph==0) = NaN;
gth(gth==0) = NaN;
gz(gz==0) = NaN;
save tsu_ray_sp.mat gx gy gph gth gz