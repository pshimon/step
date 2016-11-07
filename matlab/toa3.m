function toa3

   f = figure('Visible','off','Position',[100,200,600,400]);
   p1 = uipanel('Parent',f,'Title','Magnet','Units','normalized','Position',[0.0,0.72,0.2,0.25]);
   t1=uicontrol('Parent',p1,'style','text', 'string','height','Units','normalized','position',[0 0.0 0.6 0.3]);
   e1=uicontrol('Parent',p1,'style','edit','string','10.0','Units','normalized','position',[0.6 0.0 0.4 0.3]);
   t2=uicontrol('Parent',p1,'style','text', 'string','radius','Units','normalized','position',[0 0.3 0.6 0.3]);
   e2=uicontrol('Parent',p1,'style','edit','string','2.0','Units','normalized','position',[0.6 0.3 0.4 0.3]);   
   t3=uicontrol('Parent',p1,'style','text', 'string','magZ','Units','normalized','position',[0 0.6 0.6 0.3]);
   e3=uicontrol('Parent',p1,'style','edit','string','1.48','Units','normalized','position',[0.6 0.6 0.4 0.3]);
   p2 = uipanel('Parent',f,'Title','Sphera','Units','normalized','Position',[0.0,0.55,0.2,0.15]);
   t4=uicontrol('Parent',p2,'style','text', 'string','mass','Units','normalized','position',[0 0.0 0.6 0.5]);
   e4=uicontrol('Parent',p2,'style','edit','string','10.0','Units','normalized','position',[0.6 0.0 0.4 0.5]);
   t5=uicontrol('Parent',p2,'style','text', 'string','radius','Units','normalized','position',[0 0.5 0.6 0.5]);
   e5=uicontrol('Parent',p2,'style','edit','string','0.1','Units','normalized','position',[0.6 0.5 0.4 0.5]);
   p3 = uipanel('Parent',f,'Title','Region','Units','normalized','Position',[0.0,0.25,0.2,0.25]);  
   t6=uicontrol('Parent',p3,'style','text', 'string','pos','Units','normalized','position',[0 0.0 0.6 0.3]);
   e6=uicontrol('Parent',p3,'style','edit','string','10.0','Units','normalized','position',[0.6 0.0 0.4 0.3]);
   t7=uicontrol('Parent',p3,'style','text', 'string','far','Units','normalized','position',[0 0.3 0.6 0.3]);
   e7=uicontrol('Parent',p3,'style','edit','string','15.0','Units','normalized','position',[0.6 0.3 0.4 0.3]);   
   t8=uicontrol('Parent',p3,'style','text', 'string','near','Units','normalized','position',[0 0.6 0.6 0.3]);
   e8=uicontrol('Parent',p3,'style','edit','string','5.0','Units','normalized','position',[0.6 0.6 0.4 0.3]);   
   ax1 = axes('Parent',f,'Units','normalized','Position',[0.26,0.16,0.7,0.34]);
   ax2 = axes('Parent',f,'Units','normalized','Position',[0.26,0.60,0.7,0.34]);  
   pb= uicontrol('style','push','Units','normalized','position',[0.0 0.0 0.2 0.08],'string','Calculate','call',{@pb_call});
   r=uicontrol('style','edit','Units','normalized','position',[0.26,0.0 0.7 0.08],'string','result');
   set(f,'Visible','on'); 
   function pb_call(source,eventdata) 
    mag_r=str2double(get(e2,'string'));
    mag_h=str2double(get(e1,'string'));
    mag_mz=str2double(get(e3,'string'))*10^4/4.0/pi;
    
    M=str2double(get(e4,'string'));
    R=str2double(get(e5,'string'));

    zmin=str2double(get(e8,'string'));
    zmax=str2double(get(e7,'string'));  
    z=str2double(get(e6,'string'));
    if(z>=zmax)|(z<=zmin)
        str=sprintf('initial position should be between boundaries');
        set(r,'string',str,'background','red');
        return
    end
    nz=101; % step 1mm
    dt=0.1; %s time step
    Rfac=0.1; % dHz/dz*r<Rfac*Hz;

    dz=(zmax-zmin)/(nz-1);
    zz=linspace(zmin, zmax,nz); %distance from  the magnet

    Hz=zeros(size(zz));
    Hz=field_calc1(zz,mag_r,mag_h); %not normalized
    plot(ax1,zz,Hz*mag_mz);
    xlabel(ax1,'Distance (cm)');
 %   ylabel(ax1,'Magnetic field (Gs)');
    title(ax1,'Magnetic field (Gs)');
    dHdz=zeros(size(zz));
    dHdz(1:nz-1)= -(Hz(2:nz)-Hz(1:nz-1))/dz;
    dHdz(nz)=dHdz(nz-1);
    j0=floor((z-zmin)/dz)+1;
    Rmax=Rfac*Hz(j0)/dHdz(j0);
    if (R>Rmax)
        str=sprintf('particle radius is too large, calculation unreliable');
        set(r,'string',str,'background','red');
        return
    end
    acc=R^3*Hz.*dHdz*mag_mz^2/M;
    plot(ax2,zz,acc)
    %xlabel(ax2,'Distance (cm)');
    %ylabel(ax2,'accelaration of small sphere(cm/s^2)');
    title(ax2,'Acceleration (cm/s^2)');
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
    str=sprintf('time of arrival(s): %f',elt);
    set(r,'string',str,'background','green');
   end 
        
   
end   

