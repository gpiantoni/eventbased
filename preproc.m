function preproc(cfg, subj)
%PREPROC do preprocessing (including rerefering) TODO

mversion = 7;
%07 12/02/02 renamed to preproc (before only reref, now everything)
%06 11/12/01 write output to cfg.log (don't use cfg.fid)
%05 11/11/16 ft_preprocessin does not accept inputfile apparently
%04 11/09/12 2nd argument for subj (and subj -> subj)
%03 11/08/19 ddir -> ddir
%02 11/08/16 compute scalpcurrentdensity
%01 11/07/21 created

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
ddir = sprintf('%s%04.f/%s/%s/', cfg.data, subj, cfg.mod, cfg.cond); % data
allfile = dir([ddir '*' cfg.endname '.mat']); % files matching a preprocessing

sens = ft_read_sens(cfg.sens.file);
sens.label = upper(sens.label);
%---------------------------%

%---------------------------%
%-loop over files
for i = 1:numel(allfile)

  %-----------------%
  %-files
  inputfile = [ddir allfile(i).name];
  [~, filename] = fileparts(allfile(i).name);
  outputfile = [ddir filename '_' mfilename];
  %-----------------%

  %-----------------%
  %-prepr
  load(inputfile); % ft_preprocessing does not really use inputfile
  cfg1 = cfg.preproc;
  cfg1.feedback = 'none';
  cfg1.outputfile = outputfile;
  ft_preprocessing(cfg1, data);
  %-----------------%

  %-----------------%
  %-scalp current density
  if strcmpi(cfg.preproc.csd.do, 'yes')
    cfg1 = [];
    cfg1.method = cfg.preproc.csd.method;
    cfg1.elec = sens;
    cfg1.feedback = 'none';
    cfg1.inputfile = outputfile; % it rewrites the same file
    cfg1.outputfile = outputfile;
    ft_scalpcurrentdensity(cfg1);
  end
  %-----------------%
  
  %-----------------%
  %-clear
  if any(strcmp(mfilename, cfg.step(cfg.clear+1)))
    delete([ddir allfile(i).name])
  end
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