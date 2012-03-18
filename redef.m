function redef(cfg, subj)
%REDEF redefine trials
% This function should be as flexible bc it uses ft_redefinetrial
% Therefore, it does not do much, but it calls a function "event2trl_xxx"
% and then it calls ft_redefinetrial. It creates a separate dataset for
% each condition. It only keeps trials which are in the good part of the data.
%
% CFG
%  .data: name of projects/PROJNAME/subjects/
%  .mod: name of the modality used in recordings and projects
%  .cond: name to be used in projects/PROJNAME/subjects/0001/MOD/CONDNAME/
%  .endname: includes previous steps '_seldata_gclean_preproc'
%  .log: name of the file and directory with analysis log
%
%  .step: all the analysis step (for cfg.clear)
%  .clear: index of cfg.step to remove from subject directory
%  .redef.event2trl: function name in PROJNAME_private which creates the correct trl based on events
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
% see also SELDATA, GCLEAN, PREPROC, REDEF

%---------------------------%
%-start log
output = sprintf('(p%02.f) %s started at %s on %s\n', ...
  subj, mfilename,  datestr(now, 'HH:MM:SS'), datestr(now, 'dd-mmm-yy'));
tic_t = tic;
%---------------------------%

%---------------------------%
%-dir and files
ddir = sprintf('%s%04.f/%s/%s/', cfg.data, subj, cfg.mod, cfg.cond); % data
allfile = dir([ddir '*' cfg.endname '.mat']); % files matching a preprocessing
%---------------------------%

%-------------------------------------%
%-loop over files
for i = 1:numel(allfile)
  
  %-----------------%
  %-define the new trl
  load([ddir allfile(i).name]) % to get the events
  output = sprintf('%s\n%s\n', output, allfile(i).name);
  
  if isempty(event)
    continue
  end
  
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
        outtmp = sprintf('   cond ''%s'', no trials left SKIP (total trials:% 4.f, discarded: % 4.f)\n', ...
          cond(c).name, size(cond(c).trl,1), numel(find(~goodtrl)));
        output = [output outtmp];
        continue
      else
        outtmp = sprintf('   cond ''%s'', final trials:% 4.f (total trials:% 4.f, discarded: % 4.f))\n', ...
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
      
    end
    
  end
  %---------------------------%
  
  %-----------------%
  %-clear
  clear data
  if any(strcmp(mfilename, cfg.step(cfg.clear+1)))
    delete([ddir allfile(i).name])
  end
  %-----------------%
  
end
%-------------------------------------%

%---------------------------%
%-end log
toc_t = toc(tic_t);
outtmp = sprintf('(p%02.f) %s ended at %s on %s after %s\n\n', ...
  subj, mfilename, datestr(now, 'HH:MM:SS'), datestr(now, 'dd-mmm-yy'), ...
  datestr( datenum(0, 0, 0, 0, 0, toc_t), 'HH:MM:SS'));
output = [output outtmp];

%-----------------%
fprintf(output)
fid = fopen([cfg.log '.txt'], 'a');
fwrite(fid, output);
fclose(fid);
%-----------------%
%---------------------------%