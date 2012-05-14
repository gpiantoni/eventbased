function [dataout output] = load_subj(cfg, type, cond)
%LOAD_SUBJ load single-subject data
% Use as:
%   [data] = load_subj(cfg, type, cond)
%
% CFG
%  .data: path of /data1/projects/PROJ/subjects/
%  .rec: REC in /data1/projects/PROJ/recordings/REC/
%  .nick: NICK in /data1/projects/PROJ/subjects/SUBJCODE/MOD/NICK/
%  .mod: modality, MOD in /data1/projects/PROJ/subjects/SUBJCODE/MOD/NICK/
%  .endname: includes preprocessing steps (e.g. '_seldata_gclean_preproc_redef')
%
% TYPE
%   erp, erpsource, pow, powcorr, powsource
%
% COND
%   a string with the name used to read the data. The file name is
%   structured as: TYPE_SUBJ_COND
%   COND should not have leading or trailing underscores
%
% DATA
%   data in the specified condition
%   source-data has three dimensions:
%   1- subject
%   2- baseline source or effect source
%   3- number of peaks
%
% Part of EVENTBASED/PRIVATE

output = '';
dataout = [];

%-----------------%
%-file for each cond
typedir = ['d' type(1:3)]; % derp, dpow or dcon
groupdir = cfg.(typedir);
condname = regexprep(cond, '*', '');

for i = 1:numel(cfg.subjall)
  subj = cfg.subjall(i);
  
  subjfile = sprintf('%s_%04d_%s.mat', type, subj, condname);
  if ~exist([groupdir subjfile], 'file')
    output = [output sprintf('%s does not exist in %s\n', subjfile, groupdir)];
    continue
  end
  
  output = [output sprintf('Loading %s, cond %s, subj %04d: %s\n', type, cond, subj, subjfile)];
  
  load([groupdir subjfile])
  
  switch type
    case {'erp'}
      dataout{i} = timelock;
    case {'pow' ,'powcorr'}
      dataout{i} = freq;
    case {'erpsource'}
      for p = 1:numel(source) % n of peaks
        dataout{i,1,p} = souPre{p};
        dataout{i,2,p} = source{p};
      end
  end
end
%-----------------%