%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Shimon Panfil: Industrial Physics and Simulations                   %%
% http://industrialphys.com                                           %%
% THE SOFTWARE IS PROVIDED "AS IS", USE IT AT YOUR OWN RISK           %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function a  = read_arr1_f32(filename)
fid=fopen(filename);
n1=fread(fid,1,'int32');
a=fread(fid,n1,'float32');
fclose(fid);
