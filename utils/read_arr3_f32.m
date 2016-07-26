%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Shimon Panfil: Industrial Physics and Simulations                   %%
% http://industrialphys.com                                           %%
% THE SOFTWARE IS PROVIDED "AS IS", USE IT AT YOUR OWN RISK           %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function a  = read_arr3_f32(filename)
fid=fopen(filename);
n1=fread(fid,1,'int32');
n2=fread(fid,1,'int32');
n3=fread(fid,1,'int32');
a=fread(fid,n1*n2*n3,'float32');
a=reshape(a,[n1,n2,n3]);
fclose(fid);
