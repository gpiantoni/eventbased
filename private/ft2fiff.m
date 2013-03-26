function ft2fiff(filename, data)
%FT2FIFF, based on FIELDTRIP2FIFF
%
% Use as
%   fieldtrip2fiff(filename, data)

% ensure that the filename has the correct extension
[pathstr,name,ext] = fileparts(filename);
if isempty(ext),
  filename = [filename, '.fif'];
elseif ~strcmp(ext, '.fif')
  error('if the filename is specified with extension, this should read .fif');
end

% ensure the mne-toolbox to be on the path
ft_hastoolbox('mne', 1);
FIFF = fiff_define_constants;

% check the input data
data   = ft_checkdata(data, 'datatype', {'timelock'}, 'feedback', 'yes');

info.meas_id.version = nan;
info.meas_id.machid  = [nan;nan];
info.meas_id.secs    = nan;
info.meas_id.usecs   = nan;
info.meas_date       = [nan;nan];

info.nchan    = numel(data.label);
info.highpass = nan;
info.lowpass  = nan;
info.dev_head_t.from = FIFF.FIFFV_COORD_DEVICE;
info.dev_head_t.to = FIFF.FIFFV_COORD_HEAD;
info.dev_head_t.trans = eye(4); % TODO: test
info.ctf_head_t = [];
info.dig      = [];
info.projs    = [];
info.comps    = [];
info.bads     = []; % TODO: pass electrodes that should not be used for source reconstruction (like interpolated channels)
info.ch_names = data.label(:)';
info.chs      = elec2fiff(data.elec, data.label);
info.sfreq     = 1./mean(diff(data.time));
info.isaverage = 1;
info.isepoched = 0;
info.iscontinuous = 0;

evoked.aspect_kind = 100;
evoked.is_smsh     = 0;
evoked.nave        = max(data.dof(:));
evoked.first       = round(data.time(1)*info.sfreq);
evoked.last        = round(data.time(end)*info.sfreq);
evoked.times       = data.time;
evoked.comment     = sprintf('FieldTrip data averaged (Gio code)');
evoked.epochs      = data.avg;

fiffdata.info   = info;
fiffdata.evoked = evoked;

fiff_write_evoked(filename, fiffdata);

%-------------------
% subfunction
function [chs] = elec2fiff(elec, label)

elec = ft_convert_units(elec, 'cm'); % I think that MNE uses cm

FIFF = fiff_define_constants;

for i = 1:numel(label)
  chs(i).ch_name = label{i};
  chs(i).coord_frame = FIFF.FIFFV_COORD_DEVICE;
  chs(i).kind = FIFF.FIFFV_EEG_CH;
  chs(i).unit = 107; % = FIFF.FIFF_UNIT_V; % volts
  chs(i).unit_mul = -6; % =  FIFF.FIFF_UNITM_MU; % micro
  
  
  chs(i).scanno = i;
  chs(i).logno = i;
  
  chs(i).coil_type = 1; % ?
  chs(i).coil_trans = [];
  chs(i).range = 1; % ?
  
  chs(i).cal = 1; % calibration ?
  
  i_elec = find(strcmp(label{i}, elec.label));
  chs(i).eeg_loc = [elec.chanpos(i_elec,:)' zeros(3,1)] / 100;
  chs(i).loc = [chs(i).eeg_loc(:); 0; 1; 0; 0; 0; 1];
  
end

