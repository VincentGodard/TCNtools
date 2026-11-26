% Generates cronuscalc_lmrc.csv (reference values for test-geomagnetism.R).
% Run from a folder next to a clone of CRONUScalc at commit cdbc8f2:
%   git clone https://bitbucket.org/cronusearth/cronus-calc
%   octave-cli cronuscalc_lmrc.m
% then prepend the comment header of the CSV file.
addpath('../cronus-calc/common/production');
c = load('../cronus-calc/common/production/pmag_consts.mat'); consts = c.pmag_consts;
lat  = [45.2  -33.9  64.1   0    -77.5  32.0  -16.5 89.0];
lon  = [5.7   18.4  -21.9   0    161.8 278.0  -68.1 -45.0];
fid = fopen('cronuscalc_lmrc_raw.csv','w'); fprintf(fid,'latitude,longitude,time,rc\n');
for k = 1:numel(lat)
  s.lat = lat(k); s.long = lon(k); s.pressure = 1013.25; s.elevation = 0; s.scaling = 'lm';
  o = get_tdsf(s, consts);
  idx = [1:80, 81:150:numel(o.tv), numel(o.tv)-[150 50 0]];
  idx = unique(idx);
  for i = idx, fprintf(fid,'%.4f,%.4f,%.1f,%.12g\n', lat(k), lon(k), o.tv(i), o.Rc_Lm(i)); end
end
fclose(fid);
printf('%d tv values, max %g\n', numel(o.tv), max(o.tv));
