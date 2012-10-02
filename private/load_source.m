function [source] = load_source(info, subj, cond)
%LOAD_SOURCE load source data from virtual electrode
% Use as:
%   [data] = load_source(cfg, subj, cond)
%
% INFO
%   .dsou: directory with SOURCE data
%  
% SUBJ
%   number indicating the subject number
%
% COND
%   condition to read the data in. It's necessary that you run SOURCE_SUBJ
%   before. COND should not have leading or trailing underscores
%
% SOURCE
%   virtual electrode for condition of interest
%
% Part of EVENTBASED/PRIVATE

%---------------------------%
%-file for each cond
type = 'source';
condname = regexprep(cond, '*', '');
subjfile = sprintf('%s_%04d_%s.mat', type, subj, condname);

if exist([info.dsou subjfile], 'file')
  load([info.dsou subjfile], type)
  
else
  source = [];
  
end
%---------------------------%
