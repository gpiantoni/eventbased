function erpmne_subj(info, opt, subj)
%ERPMNE_SUBJ: prepare FIFF and noise cov for processing in MNE
%
% INFO
%  .log: name of the file and directory to save log
%
% CFG.OPT
%  .cond*: cell with conditions (e.g. {'*cond1' '*cond2'})
%  .erp*: a structure with cfg to pass to ft_timelockanalysis 
%  .rescale: scaling factor, because MNE assumes V instead of uV (use 1e-6)
%  .cov*: two number in s, the covariance window
%  .electhr: z-axis threshold in mm to reject neck elctrodes
%
% IN:
%  data in /PROJ/subjects/SUBJ/MOD/NICK/
%  folder with freesurfer and watershed (with fiducials)
%
%-fiducials
% you can only get the fiducials using the transformation in
% .hdr.tkrvox2ras (which is actually a very basic transformation). This is
% because the SURFACES use this transformation. We don't care about
% coordinates NOW, we will do the sphere later
%
% T1 = ft_read_mri('SUBJECTS_DIR/SUBJCODE/mri/T1.mgz');
% T1.transform = T1.hdr.tkrvox2ras;
% 
% cfg = [];
% cfg.interactive = 'yes';
% ft_sourceplot(cfg, T1);
%
% write to file in SUBJECTS_DIR/SUBJCODE/mri/fiducials.txt in this format,
% where the first column refers to the uppercase labels.
% FIDNZ -2 94 -26
% FIDT10 75 -9 -63
% FIDT9 -82 -10 -60
% E257 0 -69 67
% where FIDNZ is nasion, FIDT10 is the right ear, FIDT9 is the left ear,
% E257 is Cz
% You can add more (like E128), but 4 is the minimum
%
% OUT
%  [SUBJECTS_DIR 'eeg/mne_SUBJ_COND'] : erp data in MNE format
%  [SUBJECTS_DIR 'eeg/mne_SUBJ_COND_cov'] : noise covariance in MNE format
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
%-elec location
if ~isfield(opt, 'electhr'); opt.electhr = -Inf; end
elec = elec2surf(info.sens.file, opt.SUBJECTS_DIR, subj, opt.electhr);


eegdir = sprintf('%s%04d/eeg/', opt.SUBJECTS_DIR, subj);
if ~isdir(eegdir); mkdir(eegdir); end
%-------------------------------------%

%-------------------------------------%
%-loop over conditions
for k = 1:numel(opt.cond)
  cond = opt.cond{k};
  condname = regexprep(cond, '*', '');
  
  %---------------------------%
  %-read data
  [data, badchan] = load_data(info, subj, cond);
  
  if isempty(data)
    output = sprintf('%sCould not find any file for condition %s\n', ...
      output, cond);
    continue
  end
  
  outputfile = sprintf('mne_%04d_%s', subj, condname);
  outputfile_cov = [outputfile '_cov.fif'];
  %---------------------------%
  
  %---------------------------%
  %-only use non-interpolated channels
  [channel, ielec] = intersect(elec.label, setdiff(data.label, badchan), 'stable');
  cfg = [];
  cfg.channel = channel;
  data = ft_selectdata(cfg, data);
  
  elec.label = elec.label(ielec);
  elec.chanpos = elec.chanpos(ielec,:);
  elec.elecpos = elec.elecpos(ielec,:);
  %---------------------------%
  
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
  ft2fiff([eegdir outputfile], erp_s)
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
  mne_write_cov_file([eegdir outputfile_cov], cov)
  %---------------------------%
    
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

function elec = elec2surf(sensfile, SUBJECTS_DIR, subj, electhr)

subjdir = sprintf('%s%04d/', SUBJECTS_DIR, subj);

%---------------------------%
%-bem
% use the low resolution
bnd = ft_read_headshape(sprintf('%sbem/watershed/%04d_outer_skin_surface', subjdir, subj));
vol.bnd = bnd;

%-dummy, so we call project_elec in ft_prepare_vol_sens
vol.type = 'dipoli';
vol.mat = 1;
%---------------------------%

%---------------------------%
%-elec
elec = ft_read_sens(sensfile);
elec = ft_convert_units(elec, 'mm');
elec.label = upper(elec.label);
%---------------------------%

%---------------------------%
%-remove electrodes below a certain level (neck eletrodes)
elec.label = elec.label(elec.chanpos(:,3) > electhr);
elec.elecpos = elec.elecpos(elec.chanpos(:,3) > electhr, :);
elec.chanpos = elec.chanpos(elec.chanpos(:,3) > electhr,:);
%---------------------------%

%---------------------------%
%-fiducials
%-----------------%
%-read txt
fid = fopen([subjdir 'mri/fiducials.txt'], 'r');
t = textscan(fid, '%s %f %f %f');
fclose(fid);
fidu.label = t{1};
fidu.pos = [t{2} t{3} t{4}];
%-----------------%

%-----------------%
%-create matrix
[~, i_elec, i_fidu] = intersect(elec.label, fidu.label);
B = fidu.pos(i_fidu,:); % because intersect sorts labels
A = elec.chanpos(i_elec,:);
%-----------------%

%-----------------%
%-linear algebra magic
M = [B'; ones(1, size(B, 1))] / ([A'; ones(1, size(A, 1))]);
fprintf('Transformation matrix\n')
fprintf('%6.3f\t%6.3f\t%6.3f\t%6.3f\n', M') % show matrix
%-----------------%

%-----------------%
%-apply to channels
elec.chanpos = warp_apply(M, elec.chanpos);
elec.elecpos = warp_apply(M, elec.elecpos);
%-----------------%

%-----------------%
%-feedback
h = figure('vis', 'off');
ft_plot_mesh(bnd, 'facecolor', 'none', 'edgecolor', 'white')
hold on; 
ft_plot_sens(elec, 'style', '*r')
% ft_plot_slice(T1.anatomy, 'transform', T1.hdr.tkrvox2ras, 'orientation', [0 1 0])
v = [0 90 180 270];
for i = 1:numel(v)
  view(v(i), 0)
  saveas(h, sprintf([subjdir 'mri/elecbnd%d.png'], i))
end
delete(h)
%-----------------%

%-----------------%
[~, elec] = ft_prepare_vol_sens(vol, elec);
elec.chanpos = elec.elecpos;
%-----------------%
%---------------------------%

