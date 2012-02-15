function seldata(cfg, subj)
%SELECT DATA get data from recordings and put them in subject directory
% it recreates the gosd folder for each subject
% you need to specify your own trialfun

mversion = 15;
%15 12/02/14 subj-specific channels (necessary for neckersd)
%14 12/02/03 read elec from mat as well and can force renaming of the channels
%13 12/02/02 read elec from sfp file directly
%12 12/02/02 renamed, massive changes: use own trialfun
%11 11/11/04 read elec file instead of load elec
%10 11/09/22 write output to cfg.log (don't use cfg.fid)
%09 11/09/13 use consistently upper case
%08 11/09/12 2nd argument for subj (and subj -> subj)
%07 11/08/30 proj is more standard (subj dir are 0001), but rec not yet
%06 11/08/19 more standard
%05 11/07/22 delete
%04 11/07/22 add trial info
%03 11/07/21 all the data are called data (uncomment if you want better names)
%02 11/07/21 renamed seltrl
%01 11/07/20 created

%-----------------%
%-input
if nargin == 1
  subj = cfg.subj;
end
%-----------------%

%---------------------------%
%-start log
output = sprintf('(p%02.f) %s (v%02.f) started at %s on %s\n', ...
  subj, mfilename,  mversion, datestr(now, 'HH:MM:SS'), datestr(now, 'dd-mmm-yy'));
tic_t = tic;
%---------------------------%

%---------------------------%
%-dir and files
rdir = sprintf('%s%04.f/%s/%s/', cfg.recs, subj, cfg.mod, cfg.rawd); % recording raw
ddir = sprintf('%s%04.f/%s/%s/', cfg.data, subj, cfg.mod, cfg.cond); % data
if isdir(ddir); rmdir(ddir, 's'); end
mkdir(ddir)

%-----------------%
%-read sensors (can be sfp or mat)
sens = ft_read_sens(cfg.sens.file);
sens.label = upper(sens.label); % <- EGI labels are uppercase, but the elec file is lowercase
%-----------------%
%---------------------------%

%---------------------------%
%-find raw data
allfile = dir([rdir '*' cfg.rcnd '*']);

for i = 1:numel(allfile)
  
  dataset = [rdir allfile(i).name];
  
  %-----------------%
  %-definetrials
  cfg1 = [];
  cfg1.trialfun = cfg.seldata.trialfun;
  cfg1.dataset = dataset;
  
  cfg2 = ft_definetrial(cfg1);
  
  %--------%
  %-ignore files with no events
  if all(cfg2.trl(1,1:3) == [0 0 0])
    continue
  end
  %--------%
  %-----------------%
    
  %-----------------%
  %-preprocessing
  cfg2.feedback = cfg.seldata.feedback;
  if iscell(cfg.seldata.selchan)
    cfg2.channel = cfg.seldata.selchan{subj};
  else
    cfg2.channel = cfg.seldata.selchan;
  end
  
  data = ft_preprocessing(cfg2);
  %-----------------%
  
  %-----------------%
  %-fix channels
  data.label = cfg.seldata.label;
  data.elec = sens;
  %-----------------%
  
  %-----------------%
  %-save data
  [~, filename] = fileparts(allfile(i).name);
  savename = [cfg.proj '_' filename '_' mfilename];
  save([ddir savename], 'data');
  clear data
  %-----------------%
 
end
%---------------------------%

%---------------------------%
%-end log
toc_t = toc(tic_t);
outtmp = sprintf('(p%02.f) %s (v%02.f) ended at %s on %s after %s\n\n', ...
  subj, mfilename, mversion, datestr(now, 'HH:MM:SS'), datestr(now, 'dd-mmm-yy'), ...
  datestr( datenum(0, 0, 0, 0, 0, toc_t), 'HH:MM:SS'));
output = [output outtmp];

%-----------------%
fprintf(output)
fid = fopen([cfg.log '.txt'], 'a');
fwrite(fid, output);
fclose(fid);
%-----------------%
%---------------------------%