function [log] = struct2log(cfg, outtype, lvl)
%STRUCT2LOG get field structure and write log
% use as:
%  [log] = struct2log(cfg)
% where cfg is the normal fieldtrip cfg structure
% There is one optional argument, which can be 
%   'email' log for the emails with newlines (default)
%   'csv' log for csv file, in one line, separated by commas.
% Other differences between the two is that csv does not include paths,
% cuts values longer than 30 characters and fields are sorted alphabetically

% 11/12/21 added csv option
% 11/07/22 also for function_handle
% 11/07/20 created

%-----------------%
%-input check
if nargin < 2
  outtype = 'email';
end

if nargin < 3
  lvl = 0;
end

if strcmp(outtype, 'email')
  sep = sprintf('\n');
elseif strcmp(outtype, 'csv')
  sep = ',';
end
%-----------------%

%-------------------------------------%
%-prepare log
log = '';

%---------------------------%
%-if cfg contains multiple cfg(1), cfg(2), cfg(3)
if numel(cfg) ~= 1
  
  for c = 1:numel(cfg)
    flog = struct2log(cfg(c), outtype, lvl);
    log = sprintf('%s%s%s', log, flog, sep);
  end
  return
end
%---------------------------%

%---------------------------%
%-loop over fieldnames
%-----------------%
%-define fieldnames
fn = fieldnames(cfg);

%-------%
%-sort fields for csv (less meaningful, but more consistent)
if strcmp(outtype, 'csv')
  fn = sort(fn); 
end
%-------%
%-----------------%

for i = 1:numel(fn)
  
  if strcmp(outtype, 'email')
    spaces = repmat(' ', 1, lvl * 3);
  elseif strcmp(outtype, 'csv')
    spaces = '';
  end
  
  if isstruct(cfg.(fn{i}))
    
    flog = struct2log(cfg.(fn{i}), outtype,  lvl + 1);
    log = sprintf('%s%s%s:%s%s', log, spaces, fn{i}, sep, flog);
    
  else
    
    %-----------------%
    %-actual writing
    %-------%
    %-get val
    if iscell(cfg.(fn{i}))
      val = [];
      for k = 1:numel(cfg.(fn{i}))
        val = [val tochar(cfg.(fn{i}){k})];
      end
      
    else
      val = tochar(cfg.(fn{i}));
    end
    %-------%
    
    %-------%
    %-cut if the name is too long
    if strcmp(outtype, 'csv') && numel(val) > 30
      val = sprintf('%s...(%1.f)', val(1:15), numel(val));
    end
    %-------%
    
    %-------%
    %-write if it does not contain a path
    if ~strcmp(outtype, 'csv') || isempty(strfind(val, '/')) % it's a path
      log = sprintf('%s%s%s:%s%s', log, spaces, fn{i}, val, sep);
    end
    %-------%
    %-----------------%
  end
  
end
%-----------------%
%-------------------------------------%

%-------------------------------------%
%-subfunction: to char
function [val] = tochar(fld)

  if ischar(fld) 
    val = sprintf(' %s', fld);

  elseif isnumeric(fld)
    val = sprintf(' %g', fld);
    
  elseif isa(fld, 'function_handle')
    val = sprintf(' %s', func2str(fld));
    
  else
    val = 'unrecognized input';
  end
%-------------------------------------%