function pqr = loadPqr(filename)
fh = fopen(filename,'r');
C = textscan(fh, '%s%d%s%s%d%f%f%f%f%f');

if length(C{5})>0
  xyz = [C{6} C{7} C{8}];
  q = C{9};
  r = C{10};
  resnum = C{5};
  resid = C{4};
  atomid = C{3};
  atomnum = C{2};
else
  frewind(fh);
  C = textscan(fh,'%s%d%s%s%s%d%f%f%f%f%f');
  xyz = [C{7} C{8} C{9}];
  q = C{10};
  r = C{11};
  resnum = C{6};
  resid = C{4};
  atomid = C{3};
  atomnum = C{2};
end
fclose(fh);
pqr = struct('xyz', xyz, 'q', q, 'r', r, 'resnum',resnum,'resid',0, ...
	     'atomid', 0, 'atomnum',atomnum);
pqr.resid = resid;
pqr.atomid = atomid;
