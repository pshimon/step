function elt=toa2(mag_r,mag_h,M,R,zmin,zmax,z)
mag_mzT=1.48; %T
mag_mz=mag_mzT*10^4/4.0/pi; %gauss

nz=101; % step 1mm
dt=0.1; %s time step
Rfac=0.1; % dHz/dz*r<Rfac*Hz;

dz=(zmax-zmin)/(nz-1);
zz=linspace(zmin, zmax,nz); %distance from  the magnet

Hz=zeros(size(zz));
Hz=field_calc1(zz,mag_r,mag_h); %not normalized
dHdz=zeros(size(zz));
dHdz(1:nz-1)= -(Hz(2:nz)-Hz(1:nz-1))/dz;
dHdz(nz)=dHdz(nz-1);
j0=floor((z-zmin)/dz)+1;
Rmax=Rfac*Hz(j0)/dHdz(j0);
if (R>Rmax)
    error('particle radius is too large, calculation unreliable\n');
end

acc=R^3*Hz.*dHdz*mag_mz^2/M;
v=-acc(floor((z-zmin)/dz)+1)*dt*0.5;
elt=0.0;
while 1
    elt=elt+dt;
    z=z+v*dt;
    if z<=zmin
        break
    end;
    v=v-acc(floor((z-zmin)/dz)+1)*dt;
end
end



