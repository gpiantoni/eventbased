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
%      -PEAK VOX: @(x)x == max(x)
%      -ABOVE UNCORRECTED THRESHOLD: @(x)x >= abs(norminv(0.05/2))
%    .cluster: select only one cluster
%      - 'peak': cluster with peak
%      - 'sum': cluster with higher summed stat
%      - 'largest': largest cluster
%    .expand: use neighbors of one voxel (it only works with one voxel),
%             possible values are 0, 6, 18, 26
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

    statvalue = stat{i}.stat;
    if ~stat{i}.posneg
      statvalue = -1 * statvalue;
    end
    ROI(i).name = stat{i}.name;
    selvox = feval( eval(roi.fun), statvalue);
    if ~any(selvox)
      [~, selvox] = max(statvalue);
    end
    ROI(i).pos  = stat{i}.pos(selvox, :);

    value = statvalue(selvox,:);
    
    %-----------------%
    %-cluster
    if isfield(roi, 'cluster') && ~isempty(roi.cluster)
      
      x = stat{i}.pos(:,1);
      y = stat{i}.pos(:,2);
      z = stat{i}.pos(:,3);
      
      d_x = unique(diff(sort(x)));
      d_y = unique(diff(sort(y)));
      d_z = unique(diff(sort(z)));
      
      if numel(d_x) ~= 2 || numel(d_y) ~= 2 || numel(d_z) ~= 2
        error('grid is not equally spaced')
      end
      
      x_grid = min(x):d_x(2):max(x);
      y_grid = min(y):d_y(2):max(y);
      z_grid = min(z):d_z(2):max(z);
      
      %-create matrix
      B = false(stat{i}.dim);
      M = zeros(stat{i}.dim);
      
      for p = 1:size(ROI(i).pos,1)
        x_pos = find(x_grid == ROI(i).pos(p, 1));
        y_pos = find(y_grid == ROI(i).pos(p, 2));
        z_pos = find(z_grid == ROI(i).pos(p, 3));
        B(x_pos, y_pos, z_pos) = true;
        M(x_pos, y_pos, z_pos) = value(p);
      end
      
      [L, numC] = bwlabeln(B);
      
      switch roi.cluster
        case 'peak'
          
          [~, i_max] = max(M(:));
          goodC = L(i_max);
          
        case 'sum'
          
          sizeC = zeros(1, numC);
          for c = 1:numC
            sizeC(c) = sum(sum(sum((L == c) .* M)));
          end
          [~, goodC] = max(sizeC);
          
        case 'largest'

          sizeC = zeros(1, numC);
          for c = 1:numC
            sizeC(c) = numel(find(L(:) == c));
          end
          [~, goodC] = max(sizeC);
          
      end
      C = L == goodC;
      
      inroi = find(C(:));
      ROI(i).pos = zeros(numel(inroi), 3);
      for r = 1:numel(inroi)
        [x_roi, y_roi, z_roi] = ind2sub(stat{i}.dim, inroi(r));
        ROI(i).pos(r, 1) = x_grid(x_roi);
        ROI(i).pos(r, 2) = y_grid(y_roi);
        ROI(i).pos(r, 3) = z_grid(z_roi);
      end
      
    end
    %-----------------%
    
    %-----------------%
    %-expand
    if isfield(roi, 'expand') && ~isempty(roi.expand)
      if size(ROI(i).pos,1) ~= 1
        error('You can only use the neighbors of a single voxel by expanding it, not a cluster')
      end
      
      x = unique(stat{i}.pos(:,1));
      y = unique(stat{i}.pos(:,2));
      z = unique(stat{i}.pos(:,3));

      x_up = x(find(x > ROI(i).pos(1,1),1));
      x_down = x(find(x < ROI(i).pos(1,1),1,'last'));
      y_up = y(find(y > ROI(i).pos(1,2),1));
      y_down = y(find(y < ROI(i).pos(1,2),1,'last'));
      z_up = z(find(z > ROI(i).pos(1,3),1));
      z_down = z(find(z < ROI(i).pos(1,3),1,'last'));

      %-if empty, create duplicate values, they'll be removed using
      %"unique" at the end
      if isempty(x_up); x_up = ROI(i).pos(1,1); end
      if isempty(x_down); x_down = ROI(i).pos(1,1); end
      if isempty(y_up); y_up = ROI(i).pos(1,2); end
      if isempty(y_down); y_down = ROI(i).pos(1,2); end
      if isempty(z_up); z_up = ROI(i).pos(1,3); end
      if isempty(z_down); z_down = ROI(i).pos(1,3); end
      
      if roi.expand >= 6
         ROI(i).pos(2:7,:) = repmat(ROI(i).pos(1,:), 6, 1);
         ROI(i).pos(2,1) = x_up;
         ROI(i).pos(3,1) = x_down;
         ROI(i).pos(4,2) = y_up;
         ROI(i).pos(5,2) = y_down;
         ROI(i).pos(6,3) = z_up;
         ROI(i).pos(7,3) = z_down;
      end
      
      if roi.expand >= 18
         ROI(i).pos(8:19,:) = repmat(ROI(i).pos(1,:), 12, 1);
         ROI(i).pos(8, [1 2]) = [x_up y_up];
         ROI(i).pos(9, [1 2]) = [x_up y_down];
         ROI(i).pos(10, [1 2]) = [x_down y_up];
         ROI(i).pos(11, [1 2]) = [x_down y_down];
         ROI(i).pos(12, [1 3]) = [x_up z_up];
         ROI(i).pos(13, [1 3]) = [x_up z_down];
         ROI(i).pos(14, [1 3]) = [x_down z_up];
         ROI(i).pos(15, [1 3]) = [x_down z_down];
         ROI(i).pos(16, [2 3]) = [y_up z_up];
         ROI(i).pos(17, [2 3]) = [y_up z_down];
         ROI(i).pos(18, [2 3]) = [y_down z_up];
         ROI(i).pos(19, [2 3]) = [y_down z_down];
      end
      
      if roi.expand >= 26
         ROI(i).pos(20, :) = [x_up y_up z_up];
         ROI(i).pos(21, :) = [x_up y_up z_down];
         ROI(i).pos(22, :) = [x_up y_down z_up];
         ROI(i).pos(23, :) = [x_up y_down z_down];
         ROI(i).pos(24, :) = [x_down y_up z_up];
         ROI(i).pos(25, :) = [x_down y_up z_down];
         ROI(i).pos(26, :) = [x_down y_down z_up];
         ROI(i).pos(27, :) = [x_down y_down z_down];
      end
      
      % remove duplicated rows
      ROI(i).pos = unique(ROI(i).pos, 'rows');
    end
    %-----------------%
    
  end
  %---------------------------%
  
end