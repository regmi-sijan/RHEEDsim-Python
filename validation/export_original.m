function export_original(source_dir, output_dir, compat_dir)
% Execute the unmodified numerical functions distributed with AEJC_v1_0.
% The isolated compatibility directory replaces UI-only waitbar/close calls
% and supplies the standard normal PDF without requiring a statistics package.
original_path = path;
restore_path = onCleanup(@() path(original_path));
addpath(source_dir);
addpath(compat_dir, '-begin');
if ~exist(output_dir, 'dir'), mkdir(output_dir); end
files = {'GaN', 'coordinates1', 'coordinates2', 'sqrt3BL', 'sqrt3HD'};
global Recip coordinates;
for f = 1:numel(files)
  name = files{f};
  fid = fopen(fullfile(source_dir, [name '.txt']), 'r');
  if fid < 0, error('Cannot open example'); end
  fgetl(fid);
  cell = [sscanf(fgetl(fid), '%f')'; sscanf(fgetl(fid), '%f')'];
  factors = sscanf(fgetl(fid), '%f');
  coordinates = fscanf(fid, '%f', [4 Inf])';
  fclose(fid);
  coordinates(:,1) = factors(coordinates(:,1));
  Recip = local_recip(cell);
  for angle = [0 90]
    prefix = fullfile(output_dir, sprintf('%s_angle%d', name, angle));
    S = local_distfinder([Recip [angle; 0]]);
    if S(3) == 0, error('No reciprocal vector'); end
    IM = local_IM(S);
    w = S(3)/5;
    IMfine = local_gbroaden(IM, w);
    if all(coordinates(:,4) == 0)
      streaks = repmat(IMfine(:,2)', 400, 1);
    else
      streaks = local_IM2D(IM, S, w);
    end
    streaks = flipud(streaks);
    dlmwrite([prefix '_reciprocal.txt'], Recip, 'precision', '%.17g');
    dlmwrite([prefix '_vector.txt'], S, 'precision', '%.17g');
    dlmwrite([prefix '_peaks.txt'], IM, 'precision', '%.17g');
    dlmwrite([prefix '_profile.txt'], IMfine, 'precision', '%.17g');
    dlmwrite([prefix '_streaks.txt'], streaks, 'precision', '%.17g');
    fprintf('Exported %s angle %d\n', name, angle);
  end
  [peaks, density] = local_LEEDgen(3, norm(Recip(1,:))/5, 32);
  dlmwrite(fullfile(output_dir, [name '_leed_peaks.txt']), peaks, 'precision', '%.17g');
  dlmwrite(fullfile(output_dir, [name '_leed.txt']), density, 'precision', '%.17g');
end
end
