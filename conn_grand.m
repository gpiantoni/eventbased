function grandconn(cfg)
%GRANDCONN grand connectivity analysis
% it only works with frequency-domain Granger data at the moment. It
% creates a matrix (chan_chan_subj_freq) for further analysis
% There are four ways to specify cfg.gconn.freq
%   - char: 'all' : it takes all the frequencies in the data
%   - char: 'any' : it takes every single frequency in the data
%   - numeric: it takes every single frequency between the two limits ([8 12], means each frequency between 8 and 12, so 8 9 10 11 12, five values)
%   - cell: average between the two limits ({[8 12]}, means average of all the frequencies between 8 and 12, one value)

mversion = 8;
%08 12/01/15 includes cca
%07 12/01/13 create huge gconn struct/matrix with all the info we need
%06 12/01/12 can deal with time in the 5th dimension
%05 11/12/01 can use granger, dtf, coh
%04 11/11/20 average across freq bands
%03 11/10/05 simple average for connectivity
%02 11/09/27 cfg.conn.cond -> cfg.test
%01 11/08/17 created

%---------------------------%
%-start log
output = sprintf('%s (v%02.f) started at %s on %s\n', ...
  mfilename,  mversion, datestr(now, 'HH:MM:SS'), datestr(now, 'dd-mmm-yy'));
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
condname = regexprep(cfg.test{1}, '*', ''); % it does not matter
outputfile = sprintf('%s_%s*_%s.mat', cfg.proj, cfg.conn.method, condname);
  
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
gconn.mat = zeros(numel(gconn.label), numel(gconn.label), numel(gconn.time), size(gconn.freq,1), numel(cfg.test), numel(cfg.subjall));
%-----------------%
%---------------------------%

%---------------------------%
%-loop over conditions
for k = 1:numel(cfg.test)
  
  %-----------------%
  %-file for each cond
  condname = regexprep(cfg.test{k}, '*', '');
  outputfile = sprintf('%s_%s*_%s.mat', cfg.proj, cfg.conn.method, condname);
  
  allsub = dir([cfg.dcon outputfile]);
  spcell = @(name) sprintf('%s%s', cfg.dcon, name);
  allname = cellfun(spcell, {allsub.name}, 'uni', 0);
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
save([cfg.dcon cfg.proj '_' cfg.conn.method '_grandconn'], 'gconn')
%-----------------%

%---------------------------%
%-end log
toc_t = toc(tic_t);
outtmp = sprintf('%s (v%02.f) ended at %s on %s after %s\n\n', ...
  mfilename, mversion, datestr(now, 'HH:MM:SS'), datestr(now, 'dd-mmm-yy'), ...
  datestr( datenum(0, 0, 0, 0, 0, toc_t), 'HH:MM:SS'));
output = [output outtmp];

%-----------------%
fprintf(output)
fid = fopen([cfg.log '.txt'], 'a');
fwrite(fid, output);
fclose(fid);
%-----------------%
%---------------------------%