function ROI = get_roi(info, roi)
%GET_ROI get the regions of interest
% Use as:
%  roi = get_roi(info, opt.mont.roi)
%
% where opt.mont.roi has
%  -Explicit Regions of Interest
%    .name: 'name of the region of interest'
%    .pos: position of the dipoles, as Nx3 matrix
%
%  -Use Regions from source reconstruction
%    .comp: comparison of interest
%    .type: 'erpsource' or 'powsource'
%    .fun: function to determine which dipoles you should use. You can use
%          all the info in "stat", plus .posneg which tells you if the
%          activity is mostly positive or negative
%      Examples:
%      -PEAK VOX: @(x)abs(x.stat)==max(abs(x.stat))
%


if isfield(roi, 'pos')
  ROI = roi;
  
else
 
  %---------------------------%
  %-read from file
  %-----------------%
  %-comparison of interest
  comp = regexprep(roi.comp{1}, '*', '');
  if numel(roi.comp) == 2;
    comp = [comp '_' regexprep(roi.comp{2}, '*', '')];
  end
  %-----------------%
  
  %-----------------%
  %-read the data
  switch roi.type

    case 'erpsource'
      load([info.derp 'erpsource_' comp], 'erpsource')
      stat = erpsource;

    
    case 'powsource'
      load([info.dpow 'powsource_' comp], 'powsource')
      stat = powsource;
      
  end
  %-----------------%
  %---------------------------%
  
  %---------------------------%
  %-transform into ROI
  for i = 1:numel(stat)
    
    ROI(i).name = stat{i}.name;
    selvox = feval( eval(roi.fun), stat{i});
    ROI(i).pos  = stat{i}.pos(selvox, :);
    
  end
  %---------------------------%
    
end