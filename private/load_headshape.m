function  [vol, lead, sens] = load_headshape(cfg, subj)
%LOAD_HEADSHAPE load headshape
%
% Part of EVENTBASED/PRIVATE

if strcmp(cfg.vol.type, 'template')
  load(cfg.vol.template, 'vol', 'lead', 'sens')
  
else
  mod = 'smri';
  cond = 't1';
  mdir = sprintf('%s%04.f/%s/%s/', cfg.data, subj, mod, cond); % mridata dir
  mfile = sprintf('%s_%04.f_%s_%s', cfg.rec, subj, mod, cond); % mridata
  
  load([mdir mfile '_elec.mat'], 'elec')
  sens = elec;
  load([mdir mfile '_vol_' cfg.vol.type '.mat'], 'vol')
  load([mdir mfile '_lead_' cfg.vol.type '.mat'], 'lead')
  
end