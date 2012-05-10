function  [vol, lead, sens] = load_headshape(cfg, subj)
%LOAD_HEADSHAPE load headshape
% Use as:
%   [vol, lead, sens] = load_headshape(cfg, subj)
%
% CFG
%  .vol.type: 'template' or subject-specific ('dipoli' or 'openmeeg')
%    if template
%      .vol.template: file with template containing vol, lead, sens
%    if not template
%      .data: path of /data1/projects/PROJNAME/subjects/
%      .rec: RECNAME in /data1/projects/PROJNAME/recordings/RECNAME/
%
% SUBJ
%   number indicating the subject number
%
% VOL, LEAD, SENS
%   vol, lead and sens to be used for source reconstruction
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