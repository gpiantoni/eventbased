function  [vol, lead, sens] = load_headshape(info, subj)
%LOAD_HEADSHAPE load vol, lead, sens
% Use as:
%   [vol, lead, sens] = load_headshape(cfg, subj)
%
% INFO
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

if strcmp(info.vol.type, 'template')
  load(info.vol.template, 'vol', 'lead', 'sens')
  
else
  mod = 'smri';
  cond = 't1';
  mdir = sprintf('%s%04d/%s/%s/', info.data, subj, mod, cond); % mridata dir
  mfile = sprintf('%s_%04d_%s', info.rec, subj, info.sourcespace); % mridata
  
  load([mdir mfile '_elec.mat'], 'elec')
  sens = elec;
  load([mdir mfile '_vol_' info.vol.type '.mat'], 'vol')
  load([mdir mfile '_lead_' info.vol.type '.mat'], 'lead')
  
end
