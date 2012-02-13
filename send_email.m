%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SEND EMAIL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function send_email(cfg)

% 11/01/27 with pdf attachments
% 11/01/18 adapted from meg_sleep, now uses cfg instead of loading file

output = sprintf('Analysis ended at %s on %s\n', ...
  datestr(now, 'HH:MM:SS'), datestr(now, 'dd-mmm-yy'));
fprintf(output)

fid = fopen([cfg.log '.txt'], 'a');
fwrite(fid, output);
fclose(fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% read log file
fid = fopen([cfg.log '.txt'], 'r');
mailtext = textscan(fid, '%s', 'whitespace' , '', 'BufSize', 1e6);
fclose(fid);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % find attachments
pngfind  = dir([cfg.log filesep '*.png']);
attachment = [];
cnt = 0;
for k=1:numel(pngfind)
  cnt = cnt + 1;
   attachment{cnt} = [cfg.log filesep pngfind(k).name];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Send the email
if isempty(attachment)
  send_mail_message('gpiantoni.work', cfg.proj, mailtext)
else
  send_mail_message('gpiantoni.work', cfg.proj, mailtext, attachment)
end
