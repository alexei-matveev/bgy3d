#!/usr/bin/octave -q
##
## Without -q octave will print the GNU banner on every startup. Is
## there a way to tell octave not to open a graphics window?
##
function x = cubic_root (n)
  ## x = int8 (n ^ (1 / 3)); # this  does not work for 128^3 already
  x = 0;
  while (x^3 < n)
    x = x + 1;
  endwhile
endfunction

function bgy3d_plot (path)

  ## Default interval, there is no way to extract it from the vector:
  interval = [-10.0, 10.0];

  ## Vectors from BGY (arent they in *.m format?):
  vec = load (path);

  ## Assuming all the vectors are in the same size:
  [N3, M] = size (vec);
  N = cubic_root (N3);

  ## Reshape all the vectors:
  vec3d = reshape (vec, [N, N, N]);

  ## Cut the yz plane, squeeze(v(x0, 1:y, 1:z)) actually get vx(1:z,
  ## 1:y) at plane x = x0
  vec2d = squeeze (vec3d(N / 2, :, :));

  ## Axis for border:
  x = linspace (interval(1), interval(2), N);

  ## Contourf plot bgy vectors:
  figure (1);

  contourf (x, x, vec2d, "linestyle", "none"); ## "levelsep", 0.05)
  title (path);
  caxis ([0.0, 2.0]);
  colorbar ();

  ## Append .png to the input path:
  out = [path, ".png"];

  ## This saves the figure in a format derived from the file extension:
  print (out);

  ## Now plot 1d projections along x-, y-, and z-axes throuth the middle
  ## of the grid:
  zray = squeeze (vec3d(N/2, N/2, :));
  yray = squeeze (vec3d(N/2, :, N/2));
  xray = squeeze (vec3d(:, N/2, N/2));

  plot (x, xray, "-+;x;", x, yray, "-x;y;", x, zray, "-o;z;");
  title (path);

  ## Append .png to the input path:
  out = [path, ".rays.png"];

  ## This saves the figure in a format derived from the file extension:
  print (out);
endfunction

## Each argument is a path:
arg_list = argv ();
for i = 1:nargin
  bgy3d_plot (arg_list{1});
endfor
