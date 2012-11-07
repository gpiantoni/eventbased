function r_grand(info, opt)
%R_GRAND call R subfunctions (you need own script to write csv)
%
% INFO
%  .log: name of the file and directory to save log
%
% CFG.OPT
%  .rdir*: directory with R functions
%  .rfun(1).name*: name of the R function
%  .rfun(1).args{1}*: arguments for R function
%  .rfun(1).tolog: do you want the output into the log file? (logical)
%  
%  This function calls the R function through Rscript and it writes down
%  the arguments, so formatting is important. To read the arguments in R,
%  at the beginning of the R function, write: 
%    args <- commandArgs(TRUE)
%  then call each argument with
%    args[[1]], args[[2]], args[[3]]
%  The arguments are always formatted as text, so you might need 
%    eval(parse(text=args[[1]]))
%  to actually use it.
%  Passing arguments can be very tricky, because of the transformations 
%  occuring to the strings. For example, if you want to pass:
%    c('E41', 'E32')
%  You need to write
%    cfg.opt.rfun(1).args{1} = '''c("E41", "E41")'''
%  All the marks are necessary and note the difference between ' and "
% 
%  If your function writes an output, use the cfg.opt.rfun(1).tolog option,
%  this will add an extra argument with the name of the files it writes in.
%  You can write, within the R function, with:
%    sink(args[[1]], append=TRUE)
%    # your script
%    sink()
% 
% * indicates obligatory parameter
%
% Part of EVENTBASED group-analysis
% see also ERP_SUBJ, ERP_GRAND, 
% ERPSOURCE_SUBJ, ERPSOURCE_GRAND, ERPSTAT_SUBJ, ERPSTAT_GRAND,
% POW_SUBJ, POW_GRAND, POW_GRP, POWCORR_SUBJ, POWCORR_GRAND,
% POWSOURCE_SUBJ, POWSOURCE_GRAND, POWSTAT_SUBJ, POWSTAT_GRAND,
% SOURCE_SUBJ, CONN_SUBJ, CONN_GRAND, CONN_STAT,
% R_GRAND

%---------------------------%
%-start log
output = sprintf('%s began at %s on %s\n', ...
  mfilename,  datestr(now, 'HH:MM:SS'), datestr(now, 'dd-mmm-yy'));
tic_t = tic;
%---------------------------%

%---------------------------%
%-run R functions
routput = [info.log filesep 'r_output.txt'];

for i = 1:numel(opt.rfun)

  %-----------------%
  %-prepare arguments
  args = sprintf(' %s', opt.rfun(i).args{:});

  if isfield(opt.rfun(i), 'tolog') && ~isempty(opt.rfun(i).tolog) ...
      && opt.rfun(i).tolog
    args = [args ' ' routput];
  end
  %-----------------%
  
  system(['Rscript ' opt.rdir opt.rfun(i).name ' ' args]);
end
%---------------------------%

%---------------------------%
%-read info written by R script
if isfield(opt.rfun, 'tolog') && any([opt.rfun.tolog])
  fid = fopen(routput, 'r');
  outtmp = fread(fid, '*char');
  fclose(fid);
  
  output = [output sprintf('\n%s\n', outtmp)];
end
%---------------------------%

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