function [ph,th,z,ixa,iya]=rt_rungekutta_sp(t,h,dph,dth,n,dndph,dndth,...
zetao,xo,yo,ixo,iyo,si,sj,d)
% Ray tracing in Cartesian coordinate system
% t = time array in sec
% h = time interval in sec
% dph = space interval along the longitude in radian
% dth = space interval along the colatitude in radian
% n = slowness
% dndph = slowness gradient in longitude direction
% dndth = slowness gradient in colatitude direction
% zetao = initial ray direction
% xo, yo = source location in longitude and latitude
% ixo, iyo = source location in local coordinate
% si, sj = array size
% d = bathymetry
%
% by Aditya Riadi Gusman
% Earthquake Research Institute, the University of Tokyo
% March 15, 2017
% March 21, 2017
%
raddeg = acos(-1)/(180);    % degree to radian
R = 6371000;                % earth's radius
ph(1) = xo * raddeg;        % initial longitude in radian (at the source)
th(1) = (90-yo) * raddeg;   % initial colatitude in radian (at the source)
z(1) = zetao*raddeg;        % initial wave direction

% souce location in local grid (index)
ixa(1)=ixo;
iya(1)=iyo;
ix=ixo;
iy=iyo;
for i=1:(length(t))                              % calculation loop
    % values of slowness at x,y
    dndphi = dndph(ix,iy);
    dndthi = dndth(ix,iy);
    ni = n(ix,iy);

    % main functions
    F_t = @(t,th,ph,z) 1/ni/R*cos(z);                                  % change the function as you desire
    G_t = @(t,th,ph,z) 1/ni/R/sin(th)*sin(z);
    H_t = @(t,th,ph,z) -sin(z)/ni^2/R*dndthi + ...
        cos(z)/ni^2/R/sin(th)*dndphi - ...
        1/ni/R*sin(z)*cot(th);
    
    % Runge Kuta 4th order
    k_1 = F_t(t(i),th(i),ph(i),z(i));
    l_1 = G_t(t(i),th(i),ph(i),z(i));
    M_1 = H_t(t(i),th(i),ph(i),z(i));
    
    k_2 = F_t(t(i)+0.5*h,th(i)+0.5*h*k_1,ph(i)+0.5*h*l_1,z(i)+0.5*h*M_1);
    l_2 = G_t(t(i)+0.5*h,th(i)+0.5*h*k_1,ph(i)+0.5*h*l_1,z(i)+0.5*h*M_1);
    M_2 = H_t(t(i)+0.5*h,th(i)+0.5*h*k_1,ph(i)+0.5*h*l_1,z(i)+0.5*h*M_1);
    
    k_3 = F_t(t(i)+0.5*h,th(i)+0.5*h*k_2,ph(i)+0.5*h*l_2,z(i)+0.5*h*M_2);
    l_3 = G_t(t(i)+0.5*h,th(i)+0.5*h*k_2,ph(i)+0.5*h*l_2,z(i)+0.5*h*M_2);
    M_3 = H_t(t(i)+0.5*h,th(i)+0.5*h*k_2,ph(i)+0.5*h*l_2,z(i)+0.5*h*M_2);
    
    k_4 = F_t(t(i)+h,th(i)+k_3*h,ph(i)+l_3*h,z(i)+M_3*h); % Corrected  
    l_4 = G_t(t(i)+h,th(i)+k_3*h,ph(i)+l_3*h,z(i)+M_3*h);
    M_4 = H_t(t(i)+h,th(i)+k_3*h,ph(i)+l_3*h,z(i)+M_3*h);

    th(i+1) = th(i) + (1/6)*(k_1+2*k_2+2*k_3+k_4)*h;  % main equation
    ph(i+1) = ph(i) + (1/6)*(l_1+2*l_2+2*l_3+l_4)*h;  % main equation
    z(i+1) = z(i) + (1/6)*(M_1+2*M_2+2*M_3+M_4)*h;  % main equation

    % local coordinate
    ix=ixa(1) + round((ph(i)-ph(1))/dph);        % next grid
    iy=iya(1) + round((th(i)-th(1))/dth);        % next grid
    ixa(i)=ix;
    iya(i)=iy;
    
    % boundaries
    if ix>si-1 || ix<1 || iy>sj-1 || iy<1    
        return  
    end
    
    if d(ix,iy)<10
        return
    end
    if ni >=1
        return  
    end
    
end

