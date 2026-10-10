% 
% Calc budget terms for T*h
% 
clear
homedir = getenv('HOME');
workdir = getenv('WORKDIR');
addpath(genpath([homedir '/funcs']));
addpath(genpath([homedir '/mytoolbox']));
addpath(genpath([homedir '/work_GMn']));

%% read grid and depth
grid_dir = '~xbxu/hycom/GLBb0.08/topo/';
jj_NWA = 833:1233;
ii_NWA = 264:766;
jj = jj_NWA + 1140;
ii = ii_NWA + 2348;
read_HYCOM_griddep_glob;
g = 9.806;
kk = 1:kdm;
scp2 = pscx.*pscy;
onesarry = ones(size(pscx));
flg = 1; % 1: T_x = d(Th)/dx* L [T*m]; 0: T_x = d(Th)/dx [T]
% kappa = 20; % m^2/s
%  m^2/s = [m/s] * dx
kappa_x = 0.005 * pscx;
kappa_y = 0.005 * pscy;

%------ times (monthly)
yr_s = 1958;
[mon_s, mon_e, dt] = deal(1, 12*66, 1); % [month]
tunit = ['months since ' num2str(yr_s) '-01-15 12:00:00'];
t_al = mon_s:dt:mon_e;
nt = length(t_al);
dtime_al = cftime(t_al'-1, tunit); %[nt-1]

% data dir and var name
data_dir1 = '/p/cwfs/xbxu/hycom/GLBb0.08/expt_03.9/data/meanstd/';
data_dir2 = '/p/cwfs/xbxu/hycom/GLBb0.08/expt_04.0/data/meanstd/';
fnm_template = 'EXP_archMNA.yyyym.a';

n2d = 24;
n3d = 9;
varid_surflx = 3;
varid_uflx = 3;
varid_vflx = 4;
varid_thck = 6;
varid_temp = 7;
varid_rho = 9;

%  save
save_dir = [workdir '/data_GMsub/Th_budgetterms_mon/'];
if ~exist(save_dir,'dir')
    mkdir(save_dir);
end
confg_fnm = [save_dir 'confg.mat'];
save(confg_fnm,'jj','ii','kk','yr_s','t_al');

% method = 'upwind';
method = '4th_order';

%%
tic
for it = 1:nt
    
    %------ filename for model data
    dtime = dtime_al(it);
    if year(dtime) <= 2019
        EXPstr = '039';
        data_dir = data_dir1;
    elseif year(dtime) >= 2020
        EXPstr = '040';
        data_dir = data_dir2;
    end

    fnm = replace(fnm_template,{'EXP', 'yyyy','m'},{EXPstr, num2str(year(dtime),'%04d'),...
         char('a'+month(dtime)-1)});
    fnm_ful = [data_dir fnm];
    fprintf('...Reading model data from %s\n',fnm_ful);

    %------ save fnm
    savename = [save_dir 'terms_' num2str(year(dtime),'%04d') ...
        '_',num2str(month(dtime),'%02d'),'.mat'];
    if exist(savename,'file')
        fprintf(1,'\nFile already exists, SKIP: %s\n\n',savename);
        continue
    end

    %----- 
    % 2d
    surflx = sub_var2(fnm_ful,[jdm idm],n2d,varid_surflx,jj,ii); % heat flx [W/m^2]: + into oce
    % 3d
    qu = sub_var3(fnm_ful,[jdm idm kdm],n2d,n3d,varid_uflx,jj,ii,kk); % [u*h*L] m^3/s
    qv = sub_var3(fnm_ful,[jdm idm kdm],n2d,n3d,varid_vflx,jj,ii,kk);
    h3d = sub_var3(fnm_ful,[jdm idm kdm],n2d,n3d,varid_thck,jj,ii,kk) / (1000*g);
    rhosurf = sub_var3(fnm_ful,[jdm idm kdm],n2d,n3d,varid_rho,jj,ii,1) + 1034; % kg/m^3
    temp3d = sub_var3(fnm_ful,[jdm idm kdm],n2d,n3d,varid_temp,jj,ii,kk); % [degC]

    %---- vertically summed T*h
    Th_zsum = sum(h3d.*temp3d,3);
    
    %---- surface heat flx  [deg m/s]
    %surfflx_H = surflx./rhosurf./3996; % 3996 is heat capacity [J/(kg*K)]
    surfflx_H = surflx./1029./3996; % 3996 is heat capacity [J/(kg*K)]
    
    %---- lateral flux [uhLT] 
    [uhT,vhT, div_uhT] = deal(zeros(size(h3d)));
    [fuD,fvD, divD] = deal(zeros(size(h3d)));
    for ik = 1:length(kk)
        % [m3/s T]
        [uhT(:,:,ik),vhT(:,:,ik)] = calc_TFluxes_hycom(qu(:,:,ik),qv(:,:,ik),...
            temp3d(:,:,ik),method);
        % grad: d(T*h)/dx*L [T*m]
        [dtdx,dtdy] = calc_GxGy_HYCOM(temp3d(:,:,ik),h3d(:,:,ik),...
            depth,uscx,uscy,vscx,vscy,flg);
        % diffusive flx = ka*grad [m3/s T]
        fuD(:,:,ik) = kappa_x .* dtdx;
        fvD(:,:,ik) = kappa_y .* dtdy;

        % [m/s T]
        div_uhT(:,:,ik) = calc_div_hycom(uhT(:,:,ik),vhT(:,:,ik),onesarry,onesarry,scp2,1);
        divD(:,:,ik) = calc_div_hycom(fuD(:,:,ik),fvD(:,:,ik),onesarry,onesarry,scp2,1);
    end
    % vertical sum
    uhT_zsum = sum(uhT,3, 'omitnan');
    vhT_zsum = sum(vhT,3, 'omitnan');
    div_zsum = sum(div_uhT, 3, 'omitnan');
    %
    fuD_zsum = sum(fuD,3, 'omitnan');
    fvD_zsum = sum(fvD,3, 'omitnan');
    divD_zsum = sum(divD, 3, 'omitnan');

    %------- save snapshot ------
    save(savename,'*_zsum','surfflx_H*','surflx');
    fprintf('Saving terms to %s...\n',savename);

end
toc

