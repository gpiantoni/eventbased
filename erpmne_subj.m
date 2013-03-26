function erpmne_subj(info, opt, subj)
%ERPMNE_SUBJ: prepare FIFF and noise cov for processing in MNE
%
% INFO
%  .log: name of the file and directory to save log
%  .dmne: directory for MNE data
%
% CFG.OPT
%  .cond*: cell with conditions (e.g. {'*cond1' '*cond2'})
%  .erp*: a structure with cfg to pass to ft_timelockanalysis 
%  .rescale: scaling factor, because MNE assumes V instead of uV (use 1e-6)
%  .cov*: two number in s, the covariance window (do
%
% IN:
%  data in /PROJ/subjects/SUBJ/MOD/NICK/
%
% OUT
%  [info.dmne 'mne_SUBJ_COND'] : erp data in MNE format
%  [info.dmne 'mne_SUBJ_COND_cov'] : noise covariance in MNE format
%
% * indicates obligatory parameter
%
% Part of EVENTBASED single-subject
% see also ERP_SUBJ, ERP_GRAND, 
% ERPSOURCE_SUBJ, ERPSOURCE_GRAND, ERPSTAT_SUBJ, ERPSTAT_GRAND,
% POW_SUBJ, POW_GRAND, POW_GRP, POWCORR_SUBJ, POWCORR_GRAND,
% POWSOURCE_SUBJ, POWSOURCE_GRAND, POWSTAT_SUBJ, POWSTAT_GRAND,
% SOURCE_SUBJ, CONN_SUBJ, CONN_GRAND, CONN_STAT

%---------------------------%
%-start log
output = sprintf('%s (%04d) began at %s on %s\n', ...
  mfilename, subj, datestr(now, 'HH:MM:SS'), datestr(now, 'dd-mmm-yy'));
tic_t = tic;
%---------------------------%

%-------------------------------------%
%-loop over conditions
for k = 1:numel(opt.cond)
  cond = opt.cond{k};
  condname = regexprep(cond, '*', '');
  
  %---------------------------%
  %-read data
  [data] = load_data(info, subj, cond);
  
  if isempty(data)
    output = sprintf('%sCould not find any file for condition %s\n', ...
      output, cond);
    continue
  end
  
  outputfile = sprintf('mne_%04d_%s', subj, condname);
  outputfile_cov = [outputfile '_cov.fif'];
  outputfile_elec = [outputfile '_elec.txt'];
  %---------------------------%
  
  %TODO only use non-interpolated channels

  %---------------------------%
  %-rescale
  if isfield(opt, 'rescale') && ~isempty(opt.rescale)
      for i = 1:numel(data.trial)
          data.trial{i} = data.trial{i} * opt.rescale;
      end
  end
  %---------------------------%
  
  %---------------------------%
  %-timelock analysis
  cfg = opt.erp;
  cfg.feedback = 'none';
  cfg.preproc.feedback = 'none';
  
  cfg.covariance = 'yes';
  cfg.covariancewindow = opt.cov;
  
  erp_s = ft_timelockanalysis(cfg, data);
  ft2fiff([info.dmne outputfile], erp_s)
  %---------------------------%
  
  %---------------------------%
  %-covariance
  cov = [];
  cov.diag = 1; % if only eeg
  cov.data = erp_s.cov;
  cov.names = data.label';
  cov.bads = {}; % you can specify bad channels
  cov.kind = 1; % I think this is the type of channel (FIFF.FIFFV_EEG_CH)
  cov.dim = size(cov.data,1);
  cov.nfree = numel(erp_s.cov)/2; % this is probably degrees of freedom, or non-zero elements, or something like that TODO
  cov.eig = [];
  cov.eigvec = [];
  cov.projs = [];
  mne_write_cov_file([info.dmne outputfile_cov], cov)
  %---------------------------%
  
%   %---------------------------%
%   %-elec location (maybe not necessary)
%   electxt = sprintf('%s,', erp_s.elec.label{:});
%   electxt = [electxt(1:end-1) '\n'];
%   electxt = [electxt sprintf('%8.3f,', erp_s.elec.chanpos(:,1))];
%   electxt = [electxt(1:end-1) '\n'];
%   electxt = [electxt sprintf('%8.3f,', erp_s.elec.chanpos(:,2))];
%   electxt = [electxt(1:end-1) '\n'];
%   electxt = [electxt sprintf('%8.3f,', erp_s.elec.chanpos(:,3))];
%   electxt = [electxt(1:end-1) '\n'];
% 
%   fid = fopen([info.dmne outputfile_elec], 'w');
%   fprintf(fid, electxt);
%   fclose(fid);
%   %---------------------------%
  
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
fid = fopen([info.log '.txt'], 'a');
fwrite(fid, output);
fclose(fid);
%-----------------%
%---------------------------%
