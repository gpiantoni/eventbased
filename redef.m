function redef(cfg, subj)
%REDEF redefine trials
% This function should be as flexible bc it uses ft_redefinetrial
% Therefore, it does not do much, but it calls a function "event2trl_xxx"
% and then it calls ft_redefinetrial. It creates a separate dataset for
% each condition. It only keeps trials which are in the good part of the data.
%
% CFG
%  .data: path of /data1/projects/PROJ/subjects/
%  .nick: NICK in /data1/projects/PROJ/subjects/0001/MOD/NICK/
%  .mod: modality, MOD in /data1/projects/PROJ/subjects/0001/MOD/NICK/
%  .endname: includes preprocessing steps (e.g. '_seldata_gclean')
%
%  .log: name of the file and directory to save log
%
%  .step: all the analysis step (for cfg.clear)
%  .clear: cell with the name of preprocessing steps to delete
%
%  .redef.event2trl: function name in NICK_private which creates the correct trl based on events
%
%  .preproc1: struct to pass to ft_preprocessing before cutting trials (if empty, no preprocessing)
%  .preproc2: struct to pass to ft_preprocessing after cutting trials (if empty, no preprocessing)
%  .csd.method: method to do scalp current density ('finite' or 'spline' or 'hjorth')
%
% IN
%  data in /data1/projects/PROJ/subjects/SUBJ/MOD/NICK/
%
% OUT
%  data, after preprocessing, rereferencing and cut in short trials
%  It adds the condition name in the middle of the file, before cfg.endname
%  and it appends '_redef' at the end of the filename
%
% You need to write your own function to create trials. Call the function
% something like "event2trl_XXX" and use as
%   [cond output] = event2trl_gosdtrl(cfg, event)
% where
%   cfg is cfg.redef (it also includes cfg.fsample with the sampling
%   frequency of that specific dataset)
%
%   cond is a struct with
%     .name = 'name of the condition'
%     .trl = a nX3 matrix used by ft_definetrial
%     .trialinfo = extra_trialinfo (optional)
%   output is a text for output
%
% Part of EVENTBASED preprocessing
% see also SELDATA, GCLEAN, REDEF

%---------------------------%
%-start log
output = sprintf('%s (%04d) began at %s on %s\n', ...
  mfilename, subj, datestr(now, 'HH:MM:SS'), datestr(now, 'dd-mmm-yy'));
tic_t = tic;
%---------------------------%

%---------------------------%
%-dir and files
ddir = sprintf('%s%04d/%s/%s/', cfg.data, subj, cfg.mod, cfg.nick); % data dir
allfile = dir([ddir '*' cfg.endname '.mat']); % files matching a preprocessing

%-------%
%-for CSD
if isfield(cfg.sens, 'file') && ~isempty(cfg.sens.file)
  sens = ft_read_sens(cfg.sens.file);
  sens.label = upper(sens.label);
end
%-------%
%---------------------------%

%-------------------------------------%
%-loop over files
for i = 1:numel(allfile)
  
  %-----------------%
  %-
  load([ddir allfile(i).name]) % to get the events
  output = sprintf('%s\n%s\n', output, allfile(i).name);
  
  if isempty(event)
    continue
  end
  %-----------------%
  
  %-----------------%
  %-preprocessing on the full file
  if isfield(cfg, 'preproc1') && ~isempty(cfg.preproc1)
    cfg1 = cfg.preproc1;
    cfg1.feedback = 'none';
    data = ft_preprocessing(cfg1, data);
  end
  %-----------------%
  
  %-----------------%
  %-define the new trl
  cfg.redef.fsample = data.fsample; % pass the sampling frequency as well
  [cond outtmp] = feval(cfg.redef.event2trl, cfg.redef, event);
  output = [output outtmp];
  %-----------------%
  
  %---------------------------%
  %-loop over condition
  dataorig = data;
  for c = 1:numel(cond)
    if ~isempty(cond(c).trl)
      
      %-----------------%
      %-redefine trials
      %---------%
      %-insert condition name into the condition field
      basicname = allfile(i).name(1:strfind(allfile(i).name, cfg.endname)-1); % name without cfg.endname
      outputfile = [basicname '-' cond(c).name cfg.endname '_' mfilename];
      %---------%
      
      %---------%
      %-use only trials which are part of the data
      goodtrl = false(size(cond(c).trl,1), 1);
      for t = 1:size(cond(c).trl,1)
        goodtrl(t) = any(cond(c).trl(t,1) >= dataorig.sampleinfo(:,1) & cond(c).trl(t,2) <= dataorig.sampleinfo(:,2));
      end
      trl = cond(c).trl(goodtrl,:);
      %---------%
      
      %---------%
      %-output
      if numel(find(goodtrl)) == 0
        outtmp = sprintf('   cond ''%s'', no trials left SKIP (total trials:% 4d, discarded: % 4d)\n', ...
          cond(c).name, size(cond(c).trl,1), numel(find(~goodtrl)));
        output = [output outtmp];
        continue
      else
        outtmp = sprintf('   cond ''%s'', final trials:% 4d (total trials:% 4d, discarded: % 4d)\n', ...
          cond(c).name, numel(find(goodtrl)), size(cond(c).trl,1), numel(find(~goodtrl)));
        output = [output outtmp];
      end
      %---------%
      
      cfg2 = [];
      cfg2.trl = trl; % <- after ft_rejectartifact
      data = ft_redefinetrial(cfg2, dataorig);
      
      if isfield(cond(c), 'trialinfo') && ~isempty(cond(c).trialinfo)
        data.trialinfo = cond(c).trialinfo(goodtrl,:);
      end
      
      save([ddir outputfile], 'data')
      %-----------------%
      
      %-----------------%
      %-preprocessing on the full file
      if isfield(cfg, 'preproc2') && ~isempty(cfg.preproc2)
        cfg1 = cfg.preproc2;
        cfg1.feedback = 'none';
        cfg1.inputfile = [ddir outputfile]; % it rewrites the same file
        cfg1.outputfile = [ddir outputfile];
        ft_preprocessing(cfg1);
      end
      %-----------------%
      
      %-----------------%
      %-scalp current density
      if isfield(cfg, 'csd') && isfield(cfg.csd, 'method') && ~isempty(cfg.csd.method)
        cfg1 = [];
        cfg1.method = cfg.csd.method;
        cfg1.elec = sens;
        cfg1.feedback = 'none';
        cfg1.inputfile = [ddir outputfile]; % it rewrites the same file
        cfg1.outputfile = [ddir outputfile];
        ft_scalpcurrentdensity(cfg1);
      end
      %-----------------%
      
    end
    
  end
  %---------------------------%
  
  %-----------------%
  %-clear
  previousstep = cfg.step{find(strcmp(cfg.step, mfilename))-1};
  if any(strcmp(cfg.clear, previousstep))
    delete([ddir allfile(i).name])
  end
  %-----------------%
  
end
%-------------------------------------%

%---------------------------%
%-end log
toc_t = toc(tic_t);
outtmp = sprintf('%s (%04d) ended at %s on %s after %s\n\n', ...
  mfilename, subj, datestr(now, 'HH:MM:SS'), datestr(now, 'dd-mmm-yy'), ...
  datestr( datenum(0, 0, 0, 0, 0, toc_t), 'HH:MM:SS'));
output = [output outtmp];

%-----------------%
fprintf(output)
fid = fopen([cfg.log '.txt'], 'a');
fwrite(fid, output);
fclose(fid);
%-----------------%
%---------------------------%