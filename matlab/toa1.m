%magnet parameters
mag_r=2.0; %cm
mag_h=1.0; %cm
fprintf('Magnet volume(cm^3) is %f\n',pi*mag_r^2*mag_h);
mag_mzT=1.48; %T
mag_mz=mag_mzT*10^4/4.0/pi; %gauss
%particle
M=5.0; %g mass of the paricle 
fprintf('mass of the particle(g): %f\n',M);
R=0.1; %cm radious of the particle
fprintf('radius of the particle(cm): %f\n',R);
%calculation region
zmin=5;
fprintf('near boundary(cm): %f\n',zmin);
zmax=15; %cm zin<=z<=zmax
fprintf('far boundary(cm): %f\n',zmax);
z=10.0; % cm initial point
fprintf('initial position(cm): %f\n',z);
if(z>=zmax)|(z<=zmin)
    error('initial position should be between boundaries');
end
nz=101; % step 1mm
dt=0.1; %s time step
Rfac=0.1; % dHz/dz*r<Rfac*Hz;

dz=(zmax-zmin)/(nz-1);
zz=linspace(zmin, zmax,nz); %distance from  the magnet

Hz=zeros(size(zz));
Hz=field_calc1(zz,mag_r,mag_h); %not normalized
plot(zz,Hz*mag_mz);
xlabel('Distance (cm)');
ylabel('Magnetic field (Gs)');
title('Field of magnetic cylinder');
dHdz=zeros(size(zz));
dHdz(1:nz-1)= -(Hz(2:nz)-Hz(1:nz-1))/dz;
dHdz(nz)=dHdz(nz-1);
%figure
%plot(z,dHdz,'g')
%Rmax=Rfac*min(Hz./dHdz);
j0=floor((z-zmin)/dz)+1;
Rmax=Rfac*Hz(j0)/dHdz(j0);
fprintf('maximal particle radius allowed(cm): %f\n',Rmax);
if (R>Rmax)
    error('particle radius is too large, calculation unreliable\n');
end
%R=Rmax; %radius of the particle

 figure
 acc=R^3*Hz.*dHdz*mag_mz^2/M;
 plot(zz,acc)
 xlabel('Distance (cm)');
 ylabel('accelaration of small sphere(cm/s^2)');
 title('accelaration of small sphere');
%acc=sphere_acc1(10,M,z,f)

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
fprintf('time of arrival(s): %f\n',elt);
fprintf('maximal velocity (cm/s): %f\n',v);
    



