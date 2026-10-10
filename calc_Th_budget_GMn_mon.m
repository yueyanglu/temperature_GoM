clear
homedir = getenv('HOME');
addpath(genpath([homedir '/funcs']));
addpath(genpath([homedir '/mytoolbox']));
addpath(genpath([homedir '/work_GMn']));

%% read model data
grid_dir = '/nexsan/people/xbxu/hycom/GLBb0.08/NWA/grid/';
read_HYCOM_griddep_glob_sub;
g = 9.806;
kk = 1:kdm;
g = 9.806;
plon = plon - 360; % deg E [-180, 180]
scp2 = pscx.*pscy;

%------ times (monthly)
yr_s = 1958;
[mon_s, mon_e, dt] = deal(1, 12*66, 1); % [month]
tunit = ['months since ' num2str(yr_s) '-01-15 12:00:00'];
t_al = mon_s:dt:mon_e;
nt = length(t_al);
dtime_al = cftime(t_al'-1, tunit); %[nt-1]

dt_sec = dt*30*24*3600;
terms_dir = [homedir '/data_GMsub/Th_budgetterms_mon/'];

%% region
fnm = [homedir '/work_GMn/data/ids_regs/index_USregGMn_GLBsub.mat'];

dsreg = load(fnm);
ids_do = dsreg.ids_reg; % ids_GMn

% scp2(isnan(depth)) = NaN;
vol2d = scp2 .* depth; % m^3, vol of each full water col (x,y)
dA = scp2(ids_do); % m^2
vol_GMn = sum(vol2d(ids_do)); % V = int(H*dA)
disp(['vol_GMn=',num2str(vol_GMn)]); % 20812815351158.57 for GMn

if length(dsreg.jj) ~= length(jj) || length(dsreg.ii) ~= length(ii)
    error('jj in this script does not match the shallow region selection !')
end

%%
[dTdt, conv, surf,  diffus,  tend, T_meanV] = deal(zeros(nt,1));

for it = 1:nt
    %------ filename
    dtime = dtime_al(it);
    fnm = [terms_dir 'terms_' num2str(year(dtime),'%04d') '_',num2str(month(dtime),'%02d'),'.mat'];
    fprintf('Reading files from %s...\n',fnm);
    ds = load(fnm);

    % Th budget [c/s] -
    f2d = ds.Th_zsum(ids_do).*dA;
    T_meanV(it) = 1/vol_GMn * sum(f2d,'all');

    f2d = -ds.div_zsum(ids_do).*dA;
    conv(it) = 1/vol_GMn * sum(f2d,'all');
    f2d = ds.surfflx_H(ids_do).*dA;
    surf(it) = 1/vol_GMn * sum(f2d,'all');
    f2d = ds.divD_zsum(ids_do).*dA;
    diffus(it) = 1/vol_GMn * sum(f2d,'all');
end

% time int [deg]
lateral = cumsum(conv*dt_sec);
surf_int = cumsum(surf*dt_sec);
diffus_int = cumsum(diffus*dt_sec);

savenm = [homedir '/work_GMn/data/Th_budget_globsub_reg0_mon.mat'];
save(savenm,'conv','diffus','surf*','T_meanV','lateral','tend',...
      '*_int','fnm','dtime_al','t_al','nt');


%% plot 2d
% 
figure
pcolor(plon,plat,ds.surflx); shading flat
axis square
hold on
colorbar
caxis([-1 1]*2e2)


