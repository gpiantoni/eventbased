function [output data1 data2] = load_subj(info, type, cond)
%LOAD_SUBJ load single-subject data
% Use as:
%   [data] = load_subj(info, type, cond)
%
% INFO
%   .subjall: index of subjects to analyze (vector)
%   .derp: directory with ERP data
%   .dpow: directory with POW data
%   .dcor: directory with correlation data
%
% TYPE
%   erp, erpsource, erpstat, pow, powcorr, powsource, powstat
%
% COND
%   a string with the name used to read the data. The file name is
%   structured as: TYPE_SUBJ_COND
%   COND should not have leading or trailing underscores
%   It can also be a cell with multiple COND
%
% DATA
%   data in the specified condition
%   source-data has three dimensions:
%   1- subject
%   2- (1) source at baseline or (2) source of interest
%   3- number of peaks
%
% Part of EVENTBASED/PRIVATE

if ischar(cond)
  
  %-------------------------------------%
  %-one condition
  %-----------------%
  %-read the data
  [data1 output] = read_data(info, type, cond);
  %-----------------%
  
  %-----------------%
  %-if data is completely empty
  if isempty(data1)
    output = sprintf('%s data in condition %s does not exist!\n', ...
      type, cond);
    return
  end
  %-----------------%
  
  %-----------------%
  %-check if datasets are missing
  nodata = cellfun(@isempty, data1);
  if any(nodata,2)
    output = sprintf('%s!!! WARNING: in condition %s, no data for subjects: %s !!!\n', ...
      output, cond, sprintf(' %d', info.subjall(nodata)));
    data1 = data1(~nodata);
  end
  %-----------------%
  %-------------------------------------%
  
else 
  
  %-------------------------------------%
  %-two conditions
  %-----------------%
  %-read the data
  [data1 output] = read_data(info, type, cond{1});
  [data2 outtmp] = read_data(info, type, cond{2});
  output = [output outtmp];
  %-----------------%
  
  %-----------------%
  %-if data is completely empty
  if isempty(data1) 
    output = sprintf('%s data in condition %s does not exist!\n', ...
      type, cond{1});
    return
  end
  
  if isempty(data2) 
    output = sprintf('%s data in condition %s does not exist!\n', ...
      type, cond{2});
    return
  end
  %-----------------%
  
  %-----------------%
  %-check if datasets are missing
  nodata1 = cellfun(@isempty, data1);
  nodata2 = cellfun(@isempty, data2);
  nodata = nodata1 | nodata2;
  if any(nodata)
    output = sprintf('%s!!! WARNING: in condition %s and %s, missing data for subjects: %s !!!\n', ...
      output, cond{1}, cond{2}, sprintf(' %04d', info.subjall(nodata)));
    data1 = data1(~nodata);
    data2 = data2(~nodata);
  end
  %-----------------%
  
end
%-------------------------------------%

%-------------------------------------%
function [dataout output] = read_data(info, type, cond)
%read_data read the single subject data

%---------------------------%
%-file for each cond
output = '';
clear dataout
typedir = ['d' type(1:3)]; % derp, dpow or dcon
groupdir = info.(typedir);
condname = regexprep(cond, '*', '');

for i = 1:numel(info.subjall)
  subj = info.subjall(i);
  
  subjfile = sprintf('%s_%04d_%s.mat', type, subj, condname);
  if ~exist([groupdir subjfile], 'file')
    output = [output sprintf('%s does not exist in %s\n', subjfile, groupdir)];
    continue
  end
  
  output = [output sprintf('Loading %s, cond %s, subj %04d: %s\n', type, cond, subj, subjfile)];
  
  load([groupdir subjfile])
  
  switch type
    case 'erp'
      dataout{i} = erp_s;
      
    case 'pow'
      dataout{i} = pow_s;

    case 'powcorr'
      dataout{i} = powcorr_s;
      
    case 'conn'
      dataout{i} = conn_s;
      
    case 'erpsource'
      for p = 1:numel(erpsource_s_A) % n of peaks
        dataout(i,1,p,:) = erpsource_s_B; % baseline
        dataout(i,2,p,:) = erpsource_s_A; % of interest
      end
      
    case 'powsource'
      for p = 1:size(powsource_s_A,1) % n of peaks
        dataout(i,1,:,:) = powsource_s_B; % baseline
        dataout(i,2,:,:) = powsource_s_A; % of interest
      end
      
    case 'erpstat'
      for p = 1:numel(erpstat_s_A) % n of peaks
        dataout{i,1,p} = erpstat_s_B{p}; % baseline
        dataout{i,2,p} = erpstat_s_A{p}; % of interest
      end
      
    case 'powstat'
      for p = 1:numel(powstat_s_A) % n of peaks
        dataout{i,1,p} = powstat_s_B{p}; % baseline
        dataout{i,2,p} = powstat_s_A{p}; % of interest
      end
  end
end

if ~exist('dataout', 'var')
  dataout = [];
end
%---------------------------%
%-------------------------------------%