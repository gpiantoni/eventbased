function conn_grand(info, opt)
%CONN_GRAND connectivity analysis across subjects
%
% INFO
%  .log: name of the file and directory with analysis log
%  .dcon: directory with CONN data
%
% CFG.OPT
%  .conn.method*: (if .conn.type == 'cca') 'gc'; (if .conn.type == 'ft')
%                'coh', 'csd', 'plv', 'powcorr', 'amplcorr', 'ppc', 'psi'
%                (symmetric) or 'granger', 'dtf', 'pdc' (asymmetric)
%  .cond*: cell with conditions (e.g. {'*cond1' '*cond2'})'
%  .avgfreq: how you want to average frequency bands
%    if 'all': it takes all the frequency in the data (one value)
%    if 'any': it takes each frequency in the data (can be alot alot)
%    if two scalar: it takes each frequency between the extremes ([8 12], means each frequency between 8 and 12, so 8 9 10 11 12, five values)
%    if a cell with scalar: it takes the average between the two limits ({[8 12]}, means average of all the frequencies between 8 and 12, one value)
%    (necessary only for conn.method = 'ft')
%  .conn.toi*: vector with time points to run connectivity on
%
% IN
%  [info.dcon 'conn_SUBJ_COND'] 'conn_subj': power analysis for single-subject
%
% OUT
%  [info.dcon COND_CONNMETHOD_GRANDCONN]: a matrix with all connectivity measures (chan X chan X time X freq X test X subj)
%
% * indicates obligatory parameter
%
% Part of EVENTBASED group-analysis
% see also ERP_SUBJ, ERP_GRAND,
% ERPSOURCE_SUBJ, ERPSOURCE_GRAND, ERPSTAT_SUBJ, ERPSTAT_GRAND,
% POW_SUBJ, POW_GRAND, POW_GRP, POWCORR_SUBJ, POWCORR_GRAND,
% POWSOURCE_SUBJ, POWSOURCE_GRAND, POWSTAT_SUBJ, POWSTAT_GRAND, SOURCE_SUBJ, 
% CONN_SUBJ, CONN_GRAND, CONN_PLOT_TIME, CONN_PLOT_TIMEFREQ, CONN_STAT

%---------------------------%
%-start log
output = sprintf('%s began at %s on %s\n', ...
  mfilename,  datestr(now, 'HH:MM:SS'), datestr(now, 'dd-mmm-yy'));
tic_t = tic;
%---------------------------%

%---------------------------%
%-preparation
%-----------------%
%-define name of spectrum
switch opt.conn.method
  
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
    
  case 'gc'
    %-cca toolbox
    spctrm = 'gc';
    opt.avgfreq = 'all';
end
%-----------------%

%-----------------%
%-freq info
%-------%
%-load one example
condname = regexprep(opt.cond{1}, '*', ''); % it does not matter
outputfile = sprintf('conn_*_%s.mat', condname);

allsub = dir([info.dcon outputfile]);
load([info.dcon allsub(1).name], 'conn_s')
%-------%

if ischar(opt.avgfreq) && strcmp(opt.avgfreq, 'all')
  connfreq = [conn_s.freq(1) conn_s.freq(end)];
  
elseif ischar(opt.avgfreq) && strcmp(opt.avgfreq, 'any')
  connfreq = [conn_s.freq' conn_s.freq'];
  
elseif isnumeric(opt.avgfreq)
  freq1 = nearest(conn_s.freq, opt.avgfreq(1));
  freq2 = nearest(conn_s.freq, opt.avgfreq(2));
  connfreq = [conn_s.freq(freq1:freq2)' conn_s.freq(freq1:freq2)'];
  
elseif iscell(opt.avgfreq)
  connfreq = cell2mat(opt.avgfreq');
  
end
%-----------------%
%---------------------------%

%-----------------------------------------------%
%-read the data
%---------------------------%
%-loop over conditions
for k = 1:numel(opt.cond)
  
  cond     = opt.cond{k};
  condname = regexprep(cond, '*', '');
  
  %-----------------%
  %-file for each cond
  [outtmp, data] = load_subj(info, 'conn', cond);
  output = [output outtmp];
  if isempty(data); continue; end
  %-----------------%
  
  %---------------------------%
  %-prepare conn structure
  % MAT is chan X chan X time X freq X subj
  conn.label = conn_s.label;
  conn.time = opt.conn.toi;
  conn.freq = mat2cell(connfreq, ones(size(connfreq,1),1), [2]);
  conn.mat = nan(numel(conn.label), numel(conn.label), numel(conn.time), size(conn.freq,1), numel(info.subjall));
  %---------------------------%
  
  %---------------------------%
  %-loop over frequencies
  for f = 1:size(connfreq,1)
    
    %--------%
    %-freq of interest
    foi = nearest(data{1}.freq, connfreq(f,1)) : nearest(data{1}.freq, connfreq(f,2));
    
    if isempty(foi)
      error(['no frequencies are selected'])
    end
    %--------%
    
    %-----------------%
    %-loop over subj
    % data{i} has dimord "chan_chan_freq_time"
    for i = 1:numel(data)
      conn.mat(:, :, :, f, i) = permute(mean( data{i}.(spctrm)(:,:, foi, :), 3), [1 2 4 3 5 6]);
    end
    %-----------------%
    
  end
  %---------------------------%
  
  %-----------------%
  %-save
  save([info.dcon 'conn_' condname], 'conn')
  clear conn
  %-----------------%
  
end
%-----------------------------------------------%

%---------------------------%
%-end log
toc_t = toc(tic_t);
outtmp = sprintf('%s ended at %s on %s after %s\n\n', ...
  mfilename, datestr(now, 'HH:MM:SS'), datestr(now, 'dd-mmm-yy'), ...
  datestr( datenum(0, 0, 0, 0, 0, toc_t), 'HH:MM:SS'));
output = [output outtmp];

%-----------------%
fprintf(output)
fid = fopen([info.log '.txt'], 'a');
fwrite(fid, output);
fclose(fid);
%-----------------%
%---------------------------%
