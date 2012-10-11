#!/usr/bin/octave -q
##
## Without -q octave will print the GNU banner on every startup. Is
## there a way to tell octave not to open a graphics window?
##
function bgy3d_plot (path)

  ## Default interval, there is no way to extract it from the vector:
  interval = [-10.0, 10.0];

  ## Vectors from BGY (arent they in *.m format?):
  vec = load (path);

  ## Assuming all the vectors are in the same size:
  [N3, M] = size (vec);
  N = int8 (N3 ^ (1 / 3));

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
  len = size (path, 2);         # FIXME: is there a len() function?
  out = blanks (len + 4);
  out (1:len) = path;
  out (len+1:len+4) = ".png";

  ## This saves the figure in a format derived from the file extension:
  print (out);
endfunction

## Each argument is a path:
arg_list = argv ();
for i = 1:nargin
  bgy3d_plot (arg_list{1});
endfor
