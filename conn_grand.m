function conn_grand(cfg)
%CONN_GRAND connectivity analysis across subjects
%
% CFG
%  .cond: name to be used in projects/PROJNAME/subjects/0001/MOD/CONDNAME/
%  .dcon: directory with connectivity data
%  .log: name of the file and directory with analysis log
%  .conn.test: a cell with the condition defined by redef. 
%              It can, but need not, be identical to cfg.test
%  .conn.method: method used for connectivity
%  .conn.toi: vector with time points to run connectivity on
%  .gconn.freq:
%    if 'all': it takes all the frequency in the data (one value)
%    if 'any': it takes each frequency in the data (can be alot alot)
%    if two scalar: it takes each frequency between the extremes ([8 12], means each frequency between 8 and 12, so 8 9 10 11 12, five values)
%    if a cell with scalar: it takes the average between the two limits ({[8 12]}, means average of all the frequencies between 8 and 12, one value)
%
% OUT
%  [cfg.dcon COND_CONNMETHOD_GRANDCONN]: a matrix with all connectivity measures (chan X chan X time X freq X test X subj)
%
% Part of EVENTBASED group-analysis
% see also ERP_SUBJ, ERP_GRAND, 
% ERPSOURCE_SUBJ, ERPSOURCE_GRAND, ERPSTAT_SUBJ, ERPSTAT_GRAND,
% POW_SUBJ, POW_GRAND, POWCORR_SUBJ, POWCORR_GRAND,
% POWSOURCE_SUBJ, POWSOURCE_GRAND, POWSTAT_SUBJ, POWSTAT_GRAND,
% CONN_SUBJ, CONN_GRAND, CONN_STAT

%---------------------------%
%-start log
output = sprintf('%s started at %s on %s\n', ...
  mfilename,  datestr(now, 'HH:MM:SS'), datestr(now, 'dd-mmm-yy'));
tic_t = tic;
%---------------------------%

%---------------------------%
%-preparation
%-----------------%
%-define name of spectrum
switch cfg.conn.method
  
  %-symmetric
  case 'coh'
    spctrm = 'cohspctrm';
  case 'csd'
    spctrm = 'crsspctrm';
  case 'plv'
    spctrm = 'plvspctrm';
  case 'powcorr'
    spctrm = 'powcorrspctrm';
  case 'amplcorr'
    spctrm = 'amplcorrspctrm';
  case 'ppc'
    spctrm = 'ppcspctrm';
  case 'psi'
    spctrm = 'psispctrm';
    
    %-asymmetric
  case 'granger'
    spctrm = 'grangerspctrm';
  case 'dtf'
    spctrm = 'dtfspctrm';
  case 'pdc'
    spctrm = 'pdcspctrm';
    
    %-cca toolbox
  case 'gc'
    spctrm = 'gc';
end
%-----------------%

%-----------------%
%-create grand file
%-------%
%-load one example 
condname = regexprep(cfg.conn.test{cfg.statconn.ttest2(1)}, '*', ''); % it does not matter
outputfile = sprintf('%s_%s*_%s.mat', cfg.cond, cfg.conn.method, condname);
  
allsub = dir([cfg.dcon outputfile]);
load([cfg.dcon allsub(1).name], 'stat')
%-------%

%-------%
%-freq
if ischar(cfg.gconn.freq) && strcmp(cfg.gconn.freq, 'all')
  connfreq = [stat.freq(1) stat.freq(end)];

elseif ischar(cfg.gconn.freq) && strcmp(cfg.gconn.freq, 'any')
  connfreq = [stat.freq' stat.freq'];
  
elseif isnumeric(cfg.gconn.freq)
  freq1 = nearest(stat.freq, cfg.gconn.freq(1));
  freq2 = nearest(stat.freq, cfg.gconn.freq(2));
  connfreq = [stat.freq(freq1:freq2)' stat.freq(freq1:freq2)'];
  
elseif iscell(cfg.gconn.freq)
  connfreq = cell2mat(cfg.gconn.freq');
  
end
%-------%
%-----------------%

%-----------------%
%-prepare conn structure
% MAT is chan X chan X time X freq X test X subj
gconn.label = stat.label;
gconn.time = cfg.conn.toi;
gconn.freq = mat2cell(connfreq, size(connfreq,1), [2]);
gconn.mat = nan(numel(gconn.label), numel(gconn.label), numel(gconn.time), size(gconn.freq,1), numel(cfg.conn.test), numel(cfg.subjall));
%-----------------%
%---------------------------%

%---------------------------%
%-loop over conditions
for k = 1:numel(cfg.conn.test)
  
  %-----------------%
  %-file for each cond
  condname = regexprep(cfg.conn.test{k}, '*', '');
  subjfile = @(s) sprintf('%s%s_%s_%02.f_%s.mat', cfg.dcon, cfg.cond, cfg.conn.method, s, condname);
  allname = cellfun(subjfile, num2cell(cfg.subjall), 'uni', 0);
  
  allfiles = true(1, numel(allname));
  for i = 1:numel(allname)
    if ~exist(allname{i}, 'file')
      output = [output sprintf('%s does not exist\n', allname{i})];
      allfiles(i) = false;
    end
  end
  allname = allname(allfiles);
  %-----------------%

  %---------------------------%
  %-loop over frequencies
  for f = 1:size(connfreq,1)
    %--------%
    %-freq of interest
    foi = nearest(stat.freq, connfreq(f,1)) : nearest(stat.freq, connfreq(f,2));
    
    if isempty(foi)
      error(['no frequencies are selected'])
    end
    %--------%
    
    %-----------------%
    %-loop over subj
    for i = 1:numel(allname)
      load(allname{i}, 'stat')
      gconn.mat(:,:,:, f, k, i) = mean( stat.(spctrm)(:,:, foi, :), 3);
    end
    %-----------------%
  end
  %---------------------------%
  
end
%---------------------------%

%-----------------%
%-save
save([cfg.dcon cfg.cond '_' cfg.conn.method '_grandconn'], 'gconn')
%-----------------%

%---------------------------%
%-end log
toc_t = toc(tic_t);
outtmp = sprintf('%s ended at %s on %s after %s\n\n', ...
  mfilename, datestr(now, 'HH:MM:SS'), datestr(now, 'dd-mmm-yy'), ...
  datestr( datenum(0, 0, 0, 0, 0, toc_t), 'HH:MM:SS'));
output = [output outtmp];

%-----------------%
fprintf(output)
fid = fopen([cfg.log '.txt'], 'a');
fwrite(fid, output);
fclose(fid);
%-----------------%
%---------------------------%
