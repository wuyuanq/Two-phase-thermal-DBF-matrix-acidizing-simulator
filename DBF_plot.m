% This is the DBF_plot() function which generates a series of images to demonstrate the results.

% Input parameters:
% model: the model

% Author: Yuanqing Wu. Email: wuyuanq@gmail.com
% Last edited on May 23rd, 2025

function DBF_plot( model )

    % simplify the symbols of the model
    nx = model.nx;
    ny = model.ny;
    xs = model.xs;
    ys = model.ys;
    nt = model.nt;
    NUMFRAME = model.NUMFRAME;
    soludoc = model.soludoc;
    
    % load porosities
    poro = zeros(ny, nx);
    temp = load(model.fporotxt);
    c = 0;
    for j = 1 : ny
        for i = 1 : nx
            c = c + 1;
            poro(j,i) = temp(c);
        end
    end
    
    % load saturations
    Sw = zeros(ny, nx);
    temp = load(model.fSwtxt);
    c = 0;
    for j = 1 : ny
        for i = 1 : nx
            c = c + 1;
            Sw(j,i) = temp(c);
        end
    end
    
    % load permeabilities in x-direction
    Kxx = zeros(ny, nx);
    temp = load(model.fKxxtxt);
    c = 0;
    for j = 1 : ny
        for i = 1 : nx
            c = c + 1;
            Kxx(j,i) = temp(c);
        end
    end
    
    % load water-phase x-direction velocities
    velxw = zeros(ny, nx+1);
    temp = load(model.fvxwtxt);
    c = 0;
    for j = 1 : ny
        for i = 1 : nx+1 
            c = c + 1;
            velxw(j,i) = temp(c);
        end
    end
    
    % load oil-phase x-direction velocities
    velxn = zeros(ny, nx+1);
    temp = load(model.fvxntxt);
    c = 0;
    for j = 1 : ny
        for i = 1 : nx+1 
            c = c + 1;
            velxn(j,i) = temp(c);
        end
    end
    
    % load water-phase y-direction velocities
    velyw = zeros(ny+1, nx);
    temp = load(model.fvywtxt);
    c = 0;
    for j = 1 : ny+1
        for i = 1 : nx
            c = c + 1;
            velyw(j,i) = temp(c);
        end
    end
    
    % load oil-phase y-direction velocities
    velyn = zeros(ny+1, nx);
    temp = load(model.fvyntxt);
    c = 0;
    for j = 1 : ny+1
        for i = 1 : nx
            c = c + 1;
            velyn(j,i) = temp(c);
        end
    end
    
    % load pressures
    p = zeros(ny, nx);
    temp = load(model.fptxt);
    c = 0;
    for j = 1 : ny
        for i = 1 : nx
            c = c + 1;
            p(j,i) = temp(c);
        end
    end
    
    % load concentrations
    Cf = zeros(ny, nx);
    temp = load(model.fCftxt);
    c = 0;
    for j = 1 : ny
        for i = 1 : nx
            c = c + 1;
            Cf(j,i) = temp(c);
        end
    end
    
    % load temperatures
    Tem = zeros(ny, nx);
    temp = load(model.fTemtxt);
    c = 0;
    for j = 1 : ny
        for i = 1 : nx
            c = c + 1;
            Tem(j,i) = temp(c);
        end
    end
    
    % move the water-phase velocities on the edges into the center of the cell
    xwvec = zeros(ny, nx);
    ywvec = zeros(ny, nx);
    for j = 1 : ny
        for i = 1 : nx
            xwvec(j,i) = (velxw(j,i)+velxw(j,i+1))/2;
            ywvec(j,i) = (velyw(j,i)+velyw(j+1,i))/2;
        end
    end
    
    % compute the modulus of the velocity in the water phase 
    vwmodulus = zeros(ny, nx);
    for j = 1 : ny
        for i = 1 : nx
            vwmodulus(j,i) = sqrt(xwvec(j,i)^2+ywvec(j,i)^2);
        end
    end
    
    % move the oil-phase velocities on the edges into the center of the cell
    xnvec = zeros(ny, nx);
    ynvec = zeros(ny, nx);
    for j = 1 : ny
        for i = 1 : nx
            xnvec(j,i) = (velxn(j,i)+velxn(j,i+1))/2;
            ynvec(j,i) = (velyn(j,i)+velyn(j+1,i))/2;
        end
    end
    
    % compute the modulus of the velocity in the oil phase 
    vnmodulus = zeros(ny, nx);
    for j = 1 : ny
        for i = 1 : nx
            vnmodulus(j,i) = sqrt(xnvec(j,i)^2+ynvec(j,i)^2);
        end
    end
    
    % construct a grid
    xcenter = zeros(nx, 1);
    ycenter = zeros(ny, 1);
    for i = 1 : nx
        xcenter(i) = (xs(i) + xs(i+1))/2;
    end
    for i = 1 : ny
        ycenter(i) = (ys(i) + ys(i+1))/2;
    end
    
    % begin to draw the images
    fh = figure();
    h = title('Porosity field in porous media');
    set(h, 'fontsize', 14, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
    h = xlabel('X(m)');
    set(h, 'fontsize', 14, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
    h = ylabel('Y(m)');
    set(h, 'fontsize', 14, 'FontWeight', 'bold', 'VerticalAlignment', 'middle');
    axis equal
    hold on;
    [X,Y] = meshgrid(xcenter, ycenter); 
    contourf(X, Y, poro, 200, 'linecolor', 'none');
    t = colorbar;
    set(get(t,'title'), 'string', 'Porosity', 'Fontsize', 12);
    axis([0,xs(nx+1),0,ys(ny+1)]);
    hold off;
    saveas(fh, [soludoc, '/Figure 1_porosity.fig']);
    
   %create porosity movie
    v = VideoWriter([soludoc, '/poroMovie.avi']);
    open(v);
    isBT = 'false';
    for f = 1 : NUMFRAME+1
        if(f == 1)
            temp = load(['soln_poro_raw_',num2str(2),'.txt']);
        else
            fname = ['soln_poro_raw_',num2str(1+(f-1)*nt/NUMFRAME),'.txt'];
            if(exist(fname,'file') == 2) 
                temp = load(fname);
            else
                temp = load(model.fporotxt);
                isBT = 'true';
            end
        end
        c = 0;
        for j = 1 : ny
            for i = 1 : nx
                c = c + 1;
                poro(j,i) = temp(c);
            end
        end
        contourf(X, Y, poro, 200, 'linecolor', 'none');
        axis equal
        axis([0,xs(nx+1),0,ys(ny+1)]);
        writeVideo(v,getframe(gcf));
        if(strcmp(isBT,'true'))
            break
        end
    end 
    close(v);
    
    fh = figure();
    h = title('Saturation field in porous media');
    set(h, 'fontsize', 14, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
    h = xlabel('X(m)');
    set(h, 'fontsize', 14, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
    h = ylabel('Y(m)');
    set(h, 'fontsize', 14, 'FontWeight', 'bold', 'VerticalAlignment', 'middle');
    axis equal
    hold on;
    [X,Y] = meshgrid(xcenter, ycenter); 
    contourf(X, Y, Sw, 200, 'linecolor', 'none');
    t = colorbar;
    set(get(t,'title'), 'string', 'Saturation', 'Fontsize', 12);
    axis([0,xs(nx+1),0,ys(ny+1)]);
    hold off;
    saveas(fh, [soludoc, '/Figure 2_saturation.fig']);
    
    fh = figure();
    h = title('Permeability field in x-direction in porous media');
    set(h, 'fontsize', 14, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
    h = xlabel('X(m)');
    set(h, 'fontsize', 14, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
    h = ylabel('Y(m)');
    set(h, 'fontsize', 14, 'FontWeight', 'bold', 'VerticalAlignment', 'middle');
    axis equal
    hold on;
    [X,Y] = meshgrid(xcenter, ycenter); 
    contourf(X, Y, Kxx, 200, 'linecolor', 'none');
    t = colorbar;
    set(get(t,'title'), 'string', 'Perm', 'Fontsize', 12);
    axis([0,xs(nx+1),0,ys(ny+1)]);
    hold off;
    saveas(fh, [soludoc, '/Figure 3_permeability.fig']);

    fh = figure();
    h = title('Water-phase velocity component in x direction in porous media (unit: m/s)');
    set(h, 'fontsize', 14, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
    h = xlabel('X(m)');
    set(h, 'fontsize', 14, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
    h = ylabel('Y(m)');
    set(h, 'fontsize', 14, 'FontWeight', 'bold', 'VerticalAlignment', 'middle');
    axis equal
    hold on;
    [X,Y] = meshgrid(xs, ycenter); 
    contourf(X, Y, velxw, 200, 'linecolor', 'none');
    t = colorbar;
    set(get(t,'title'), 'string', 'velXw', 'Fontsize', 12);
    axis([0,xs(nx+1),0,ys(ny+1)]);
    hold off;
    saveas(fh, [soludoc, '/Figure 4_vxw.fig']);
    
    fh = figure();
    h = title('Water-phase velocity component in y direction in porous media (unit: m/s)');
    set(h, 'fontsize', 14, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
    h = xlabel('X(m)');
    set(h, 'fontsize', 14, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
    h = ylabel('Y(m)');
    set(h, 'fontsize', 14, 'FontWeight', 'bold', 'VerticalAlignment', 'middle');
    axis equal
    hold on;
    [X,Y] = meshgrid(xcenter, ys); 
    contourf(X, Y, velyw, 200, 'linecolor', 'none');
    t = colorbar;
    set(get(t,'title'), 'string', 'velYw', 'Fontsize', 12);
    axis([0,xs(nx+1),0,ys(ny+1)]);
    hold off;
    saveas(fh, [soludoc, '/Figure 5_vyw.fig']);
 
    fh = figure();
    h = title('Modulus of water-phase velocity in porous media (unit: m/s)');
    set(h, 'fontsize', 14, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
    h = xlabel('X(m)');
    set(h, 'fontsize', 14, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
    h = ylabel('Y(m)');
    set(h, 'fontsize', 14, 'FontWeight', 'bold', 'VerticalAlignment', 'middle');
    axis equal
    hold on;
    [X,Y] = meshgrid(xcenter, ycenter); 
    contourf(X, Y, vwmodulus, 200, 'linecolor', 'none');
    t = colorbar;
    set(get(t,'title'), 'string', '|vw|', 'Fontsize', 12);
    axis([0,xs(nx+1),0,ys(ny+1)]);
    hold off;
    saveas(fh, [soludoc, '/Figure 6_vwmodulus.fig']);

    fh = figure();
    h = title('Water-phase streamlines in porous media');
    set(h, 'fontsize', 14, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
    h = xlabel('X(m)');
    set(h, 'fontsize', 14, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
    h = ylabel('Y(m)');
    set(h, 'fontsize', 14, 'FontWeight', 'bold', 'VerticalAlignment', 'middle');
    axis equal
    hold on;
    [X,Y] = meshgrid(xcenter, ycenter);
    h = streamslice(X, Y, xwvec, ywvec, 2);
    axis([0,xs(nx+1),0,ys(ny+1)]);
    set(h, 'color', 'black');
    hold off;
    saveas(fh, [soludoc, '/Figure 7_water_streamline.fig']);
 
    fh = figure();
    h = title('Oil-phase velocity component in x direction in porous media (unit: m/s)');
    set(h, 'fontsize', 14, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
    h = xlabel('X(m)');
    set(h, 'fontsize', 14, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
    h = ylabel('Y(m)');
    set(h, 'fontsize', 14, 'FontWeight', 'bold', 'VerticalAlignment', 'middle');
    axis equal
    hold on;
    [X,Y] = meshgrid(xs, ycenter); 
    contourf(X, Y, velxn, 200, 'linecolor', 'none');
    t = colorbar;
    set(get(t,'title'), 'string', 'velXn', 'Fontsize', 12);
    axis([0,xs(nx+1),0,ys(ny+1)]);
    hold off;
    saveas(fh, [soludoc, '/Figure 8_vxn.fig']);
    
    fh = figure();
    h = title('Oil-phase velocity component in y direction in porous media (unit: m/s)');
    set(h, 'fontsize', 14, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
    h = xlabel('X(m)');
    set(h, 'fontsize', 14, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
    h = ylabel('Y(m)');
    set(h, 'fontsize', 14, 'FontWeight', 'bold', 'VerticalAlignment', 'middle');
    axis equal
    hold on;
    [X,Y] = meshgrid(xcenter, ys); 
    contourf(X, Y, velyn, 200, 'linecolor', 'none');
    t = colorbar;
    set(get(t,'title'), 'string', 'velYn', 'Fontsize', 12);
    axis([0,xs(nx+1),0,ys(ny+1)]);
    hold off;
    saveas(fh, [soludoc, '/Figure 9_vyn.fig']);
 
    fh = figure();
    h = title('Modulus of oil-phase velocity in porous media (unit: m/s)');
    set(h, 'fontsize', 14, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
    h = xlabel('X(m)');
    set(h, 'fontsize', 14, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
    h = ylabel('Y(m)');
    set(h, 'fontsize', 14, 'FontWeight', 'bold', 'VerticalAlignment', 'middle');
    axis equal
    hold on;
    [X,Y] = meshgrid(xcenter, ycenter); 
    contourf(X, Y, vnmodulus, 200, 'linecolor', 'none');
    t = colorbar;
    set(get(t,'title'), 'string', '|vn|', 'Fontsize', 12);
    axis([0,xs(nx+1),0,ys(ny+1)]);
    hold off;
    saveas(fh, [soludoc, '/Figure 10_vnmodulus.fig']);

    fh = figure();
    h = title('Oil-phase streamlines in porous media');
    set(h, 'fontsize', 14, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
    h = xlabel('X(m)');
    set(h, 'fontsize', 14, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
    h = ylabel('Y(m)');
    set(h, 'fontsize', 14, 'FontWeight', 'bold', 'VerticalAlignment', 'middle');
    axis equal
    hold on;
    [X,Y] = meshgrid(xcenter, ycenter);
    h = streamslice(X, Y, xnvec, ynvec, 2);
    axis([0,xs(nx+1),0,ys(ny+1)]);
    set(h, 'color', 'black');
    hold off;
    saveas(fh, [soludoc, '/Figure 11_oil_streamline.fig']);
 
    fh = figure();
    h = title('Pressure field in porous media (unit: Pa)');
    set(h, 'fontsize', 14, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
    h = xlabel('X(m)');
    set(h, 'fontsize', 14, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
    h = ylabel('Y(m)');
    set(h, 'fontsize', 14, 'FontWeight', 'bold', 'VerticalAlignment', 'middle');
    axis equal
    hold on;
    [X,Y] = meshgrid(xcenter, ycenter); 
    contourf(X, Y, p, 200, 'linecolor', 'none');
    t = colorbar;
    set(get(t,'title'), 'string', 'pres', 'Fontsize', 12);
    axis([0,xs(nx+1),0,ys(ny+1)]);
    hold off;
    saveas(fh, [soludoc, '/Figure 12_pressure.fig']);
    
    fh = figure();
    h = title('Concentration field in porous media (unit: mol/m^3)');
    set(h, 'fontsize', 14, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
    h = xlabel('X(m)');
    set(h, 'fontsize', 14, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
    h = ylabel('Y(m)');
    set(h, 'fontsize', 14, 'FontWeight', 'bold', 'VerticalAlignment', 'middle');
    axis equal
    hold on;
    [X,Y] = meshgrid(xcenter, ycenter); 
    contourf(X, Y, Cf, 200, 'linecolor', 'none');
    t = colorbar;
    set(get(t,'title'), 'string', 'Concen', 'Fontsize', 12);
    axis([0,xs(nx+1),0,ys(ny+1)]);
    hold off;
    saveas(fh, [soludoc, '/Figure 13_concentration.fig']);
    
    fh = figure();
    h = title('Temperature field in porous media (unit: K)');
    set(h, 'fontsize', 14, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
    h = xlabel('X(m)');
    set(h, 'fontsize', 14, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
    h = ylabel('Y(m)');
    set(h, 'fontsize', 14, 'FontWeight', 'bold', 'VerticalAlignment', 'middle');
    axis equal
    hold on;
    [X,Y] = meshgrid(xcenter, ycenter); 
    contourf(X, Y, Tem, 200, 'linecolor', 'none');
    t = colorbar;
    set(get(t,'title'), 'string', 'Tempe', 'Fontsize', 12);
    axis([0,xs(nx+1),0,ys(ny+1)]);
    hold off;
    saveas(fh, [soludoc, '/Figure 14_temperature.fig']);
    
end

