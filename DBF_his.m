% This is the DBF_his() function which generates a series of history images to demonstrate the results.

% Input parameters:
% model: the model

% Author: Yuanqing Wu. Email: wuyuanq@gmail.com
% Last edited on May 23rd, 2025

function DBF_his( model )

    fh = figure();
    strtitle = 'The average porosity progress';
    h = title(strtitle);
    set(h, 'fontsize', 16, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
    h = xlabel('Time(s)');
    set(h, 'fontsize', 16, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
    h = ylabel('Average porosity');
    set(h, 'fontsize', 16, 'FontWeight', 'bold', 'VerticalAlignment', 'middle');
    hold on;
    temp = load(model.fporohistxt);
    nt = size(temp,1);
    ts = (0:nt)*model.timeEnd/model.nt; 
    poroh = zeros(nt+1, 1);
    k = 1;
    for i = 2 : nt+1
        poroh(i) = temp(k);
        k = k + 1;
    end
    poroh(1) = poroh(2);
    plot(ts, poroh);
    saveas(fh, [model.soludoc, '/average_poro_history.fig']);
    
    fh = figure();
    strtitle = 'The average permeability progress';
    h = title(strtitle);
    set(h, 'fontsize', 16, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
    h = xlabel('Time(s)');
    set(h, 'fontsize', 16, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
    h = ylabel('Average permeability (m^2)');
    set(h, 'fontsize', 16, 'FontWeight', 'bold', 'VerticalAlignment', 'middle');
    hold on;
    temp = load(model.fKxxhistxt);
    Kxxh = zeros(nt+1, 1);
    k = 1;
    for i = 2 : nt+1
        Kxxh(i) = temp(k);
        k = k + 1;
    end
    Kxxh(1) = Kxxh(2);
    plot(ts, Kxxh);
    saveas(fh, [model.soludoc, '/average_Kxx_history.fig']);
    
    fh = figure();
    strtitle = 'The average av progress';
    h = title(strtitle);
    set(h, 'fontsize', 16, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
    h = xlabel('Time(s)');
    set(h, 'fontsize', 16, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
    h = ylabel('Average interfacial surface area per unit volume (m^(-1))');
    set(h, 'fontsize', 16, 'FontWeight', 'bold', 'VerticalAlignment', 'middle');
    hold on;
    temp = load(model.favhistxt);
    avh = zeros(nt+1, 1);
    k = 1;
    for i = 2 : nt+1
        avh(i) = temp(k);
        k = k + 1;
    end
    avh(1) = avh(2);
    plot(ts, avh);
    saveas(fh, [model.soludoc, '/average_av_history.fig']);
    
    fh = figure();
    strtitle = 'The average pressure progress';
    h = title(strtitle);
    set(h, 'fontsize', 16, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
    h = xlabel('Time(s)');
    set(h, 'fontsize', 16, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
    h = ylabel('Average pressure (Pa)');
    set(h, 'fontsize', 16, 'FontWeight', 'bold', 'VerticalAlignment', 'middle');
    hold on;
    temp = load(model.fphistxt);
    ph = zeros(nt+1, 1);
    k = 1;
    for i = 2 : nt+1
        ph(i) = temp(k);
        k = k + 1;
    end
    ph(1) = ph(2);
    plot(ts, ph);
    saveas(fh, [model.soludoc, '/average_p_history.fig']);
    
    fh = figure();
    strtitle = 'The average concentration progress';
    h = title(strtitle);
    set(h, 'fontsize', 16, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
    h = xlabel('Time(s)');
    set(h, 'fontsize', 16, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
    h = ylabel('Average concentration (mol/m^3)');
    set(h, 'fontsize', 16, 'FontWeight', 'bold', 'VerticalAlignment', 'middle');
    hold on;
    temp = load(model.fCfhistxt);
    Cfh = zeros(nt+1, 1);
    k = 1;
    for i = 2 : nt+1
        Cfh(i) = temp(k);
        k = k + 1;
    end
    Cfh(1) = Cfh(2);
    plot(ts, Cfh);
    saveas(fh, [model.soludoc, '/average_Cf_history.fig']);
    
    
    fh = figure();
    strtitle = 'The average temperature progress';
    h = title(strtitle);
    set(h, 'fontsize', 16, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
    h = xlabel('Time(s)');
    set(h, 'fontsize', 16, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
    h = ylabel('Average temperature (K)');
    set(h, 'fontsize', 16, 'FontWeight', 'bold', 'VerticalAlignment', 'middle');
    hold on;
    temp = load(model.fTemhistxt);
    Temh = zeros(nt+1, 1);
    k = 1;
    for i = 2 : nt+1
        Temh(i) = temp(k);
        k = k + 1;
    end
    Temh(1) = Temh(2);
    plot(ts, Temh);
    saveas(fh, [model.soludoc, '/average_Tem_history.fig']);
    
    fh = figure();
    strtitle = 'The average flux progress';
    h = title(strtitle);
    set(h, 'fontsize', 16, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
    h = xlabel('Time(s)');
    set(h, 'fontsize', 16, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
    h = ylabel('Average flux (m/s)');
    set(h, 'fontsize', 16, 'FontWeight', 'bold', 'VerticalAlignment', 'middle');
    hold on;
    temp = load(model.fqhistxt);
    qh = zeros(nt+1, 1);
    k = 1;
    for i = 2 : nt+1
        qh(i) = temp(k);
        k = k + 1;
    end
    qh(1) = qh(2);
    plot(ts, qh);
    saveas(fh, [model.soludoc, '/average_q_history.fig']);
    
    fh = figure();
    strtitle = 'The average input pressure progress';
    h = title(strtitle);
    set(h, 'fontsize', 16, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
    h = xlabel('Time(s)');
    set(h, 'fontsize', 16, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
    h = ylabel('Average input pressure (Pa)');
    set(h, 'fontsize', 16, 'FontWeight', 'bold', 'VerticalAlignment', 'middle');
    hold on;
    temp = load(model.flphistxt);
    lph = zeros(nt+1, 1);
    k = 1;
    for i = 2 : nt+1
        lph(i) = temp(k);
        k = k + 1;
    end
    lph(1) = lph(2);
    plot(ts, lph);
    saveas(fh, [model.soludoc, '/average_lp_history.fig']);
      
end
