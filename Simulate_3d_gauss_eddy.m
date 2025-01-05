%% Gauss distribution
clear
Save_dir = ['3d_test_pic/'];
% Initialize the parameter
p0 = -1*10000;      % -1 dbar
g = 9.8;
L = 300*10^3;
H = 400;
rho0 = 1.02*1000;
f = 2*2*pi/(23.933*3600)*sind(30);

% Initialize experiment area.
dx = 2; dy = 2; dz =2;  % units: km, km ,m
lon_x = -400:dx:400; lon_x = lon_x.*1000;
lat_y = -400:dy:400; lat_y = lat_y.*1000;
dep_z = -800:dz:0; dep_z = dep_z(end:-1:1);

% Process the grid
[X,Y,Z] = meshgrid(lon_x, lat_y, dep_z);
R = sqrt(X.^2+Y.^2); R_std = R./L; Z_std = Z./H;

% Modelling the Gauss distribution
F = exp(-1.*R_std.^2);
G = exp(-1.*Z_std.^2);

% Caculate the results
P_3d = p0.*F.*G;
rho_3d = 2.*p0./g./H.*F.*Z_std.*G;
N_2_3d = p0./rho0./H.*F.*G.*(4*Z_std.^2-2);
V_3d = -2.*p0./rho0./f./L.*R_std.*F.*G;

% % simple validation
% contour(lon_x,dep_z, squeeze(P_3d(:,201,:)./10000)');colorbar
% contour(lon_x,dep_z, squeeze(V_3d(:,201,:))');colorbar

%% Process the physical field.
u_save = nan(size(V_3d)); v_save = nan(size(V_3d));
eta = P_3d(:,:,1)./rho0./g;

for i = 1:numel(size(V_3d,3))
    [fx, fy] = gradient(P_3d(:,:,i));
    [theta, rho] = cart2pol(fx, fy);
    theta = theta+pi/2;
    u_save(:, :, i) = V_3d(:,:,i).*cos(theta);
    v_save(:, :, i) = V_3d(:,:,i).*sin(theta);
end

%% Save the data
for i = 1:2
    run_time = datestr(datetime(2011,1,1,0,0,0)+hours(24).*(i-1), 'yyyymmdd');
    fnt = [Save_dir, 'gauss_',run_time,'_T.nc']; delete(fnt)
    fnu = [Save_dir, 'gauss_',run_time,'_U.nc']; delete(fnu)
    fnv = [Save_dir, 'gauss_',run_time,'_V.nc']; delete(fnv)

    nccreate(fnt, 'lon', 'Dimensions', {'lon', numel(lon_x)});
    nccreate(fnt, 'lat', 'Dimensions', {'lat', numel(lat_y)});
    nccreate(fnt, 'dep', 'Dimensions', {'dep', numel(dep_z)});
    nccreate(fnt, 'Eta', 'Dimensions', {'lon', 'lat'});
    nccreate(fnt, 'Pre_anomaly', 'Dimensions', {'lon', 'lat', 'dep'});
    ncwrite(fnt, 'lon', lon_x); ncwrite(fnt, 'lat', lat_y); ncwrite(fnt, 'dep', dep_z);
    ncwrite(fnt, 'Pre_anomaly', P_3d);
    ncwrite(fnt, 'Eta', eta);  

    nccreate(fnu, 'lon', 'Dimensions', {'lon', numel(lon_x)});
    nccreate(fnu, 'lat', 'Dimensions', {'lat', numel(lat_y)});
    nccreate(fnu, 'dep', 'Dimensions', {'dep', numel(dep_z)});
    nccreate(fnu, 'U', 'Dimensions', {'lon', 'lat', 'dep'});
    ncwrite(fnu, 'lon', lon_x); ncwrite(fnu, 'lat', lat_y); ncwrite(fnu, 'dep', dep_z); 
    ncwrite(fnu, 'U', u_save);  

    nccreate(fnv, 'lon', 'Dimensions', {'lon', numel(lon_x)});
    nccreate(fnv, 'lat', 'Dimensions', {'lat', numel(lat_y)});
    nccreate(fnv, 'dep', 'Dimensions', {'dep', numel(dep_z)});
    nccreate(fnv, 'V', 'Dimensions', {'lon', 'lat', 'dep'});
    ncwrite(fnv, 'lon', lon_x); ncwrite(fnv, 'lat', lat_y); ncwrite(fnv, 'dep', dep_z); 
    ncwrite(fnv, 'Eta', v_save);  
end; clear fnt fnu fnv i

%% Process the Grid

dx_matrix = repmat(dx*1000, numel(lon_x), numel(lat_y), numel(dep_z));
dy_matrix = repmat(dy*1000, numel(lon_x), numel(lat_y), numel(dep_z));
dz_matrix = repmat(dz, numel(lon_x), numel(lat_y), numel(dep_z));
dep_name = repmat(numel(dep_z), numel(lon_x), numel(lat_y));

%%

fn_topo = [Save_dir, 'topography.nc']; delete(fn_topo);
nccreate(fn_topo, 'lon', 'Dimensions', {'lon', numel(lon_x)});
nccreate(fn_topo, 'lat', 'Dimensions', {'lat', numel(lat_y)});
nccreate(fn_topo, 'dep', 'Dimensions', {'dep', numel(dep_z)});
nccreate(fn_topo, 'U', 'Dimensions', {'lon', 'lat', 'dep'});
ncwrite(fn_topo, 'lon', lon_x); ncwrite(fn_topo, 'lat', lat_y); ncwrite(fn_topo, 'dep', dep_z); 
ncwrite(fn_topo, 'U', u_save);  





