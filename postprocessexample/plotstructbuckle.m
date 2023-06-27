clc;clear;close all;cd 'C:\Users\leixi\OneDrive\Desktop';
FONT = 'Arial';
FONTSIZE = 10;
pWidth = 4; % inches
pHeight = 3;
colpos = [247 148 30;0 166 81;237 28 36;0 174 239; 0 0 0]/255; % colors

includecl = 1;
Fmax = [2]
Tcycle = [1,1,.5,.3];
Tcycle(100)= 5;
casenum =100%100%100%4;
Ltott= (95.4+2)*2/1000;
%% forced vibration freq2
add = '';
% add='_u015';
% add = '_force2';
%  add = '_u025';
locc=1;
textnum =['uffreq',num2str(casenum),add,'.txt'];
% cll = load(['cluffreq',num2str(2),'.txt']);
% cdd = load(['cduffreq',num2str(casenum),add,'.txt']);
% disp = load(['dispuffreq',num2str(casenum),add,'.txt']);
loccdd = load(['loccduffreq',num2str(casenum),add,'.txt']);
if includecl ==1
    loccll = load(['loccluffreq',num2str(casenum),add,'.txt']);
end
locdisp = load(['locdispuffreq',num2str(casenum),add,'.txt']);
locx = load(['locxffreq',num2str(casenum),add,'.txt']);
locy = load(['locyffreq',num2str(casenum),add,'.txt']);
timeall = locx(:,1);
indexstart = find(abs(timeall-1.01) ==min(abs(timeall-1.01)));
FF=  Fmax*2*(2*mod(floor(timeall /Tcycle(casenum)),2)-1);
% return
selnode =[37:-4:5,3,6:4:38]+1;
selt =1;

Tcycled = 2*Tcycle(casenum);
selt1= indexstart;
seltall = [selt1, round(selt1+Tcycled/4*100), round(selt1+2*Tcycled/4*100), round(selt1+3*Tcycled/4*100)];
sampleselnode = selnode([3,6,9]);

h1 = figure;
plot(locx(selt,selnode),locy(selt,selnode),'o-')
xlabel('x (mm)', 'Fontname',FONT,'FontSize',FONTSIZE);
ylabel('y (mm)', 'Fontname',FONT,'FontSize',FONTSIZE);
AX = legend('Original Pose');
LEG = findobj(AX,'type','text');
set(LEG,'Fontname',FONT,'FontSize',FONTSIZE);
set(gca,'Fontname', FONT,'FontSize',FONTSIZE);
set(gcf, 'PaperUnits','inches', 'PaperPosition',[0 0 pWidth pHeight], ...
    'PaperSize', [pWidth pHeight]);
xlim([0 65]);
ylim([-100 100]);
saveas(h1, ['origpose',textnum,'.pdf']);
% return

h1 = figure;
count = 1;
for selt = seltall
    h = plot(locx(selt,selnode),locy(selt,selnode), 'Color', colpos(count,:));
    hold on
    count = 1+count;
end
xlabel('x (mm)', 'Fontname',FONT,'FontSize',FONTSIZE);
ylabel('y (mm)', 'Fontname',FONT,'FontSize',FONTSIZE);
legend('t/T = 0','t/T = \pi/4','t/T = \pi/2','t/T = 3\pi/4' );
LEG = findobj(AX,'type','axes');
set(LEG,'Fontname',FONT,'FontSize',FONTSIZE);
set(gca,'Fontname', FONT,'FontSize',FONTSIZE);
set(gcf, 'PaperUnits','inches', 'PaperPosition',[0 0 pWidth pHeight], ...
    'PaperSize', [pWidth pHeight]);
xlim([0 65]);
ylim([-100 100]);
saveas(h1, ['1cyclepose',textnum,'.pdf']);

%% the purpose is to compare the displacement and drag
for selt= seltall
    h1= figure;
    subplot(1,2,1)
    h= plot(locx(selt,selnode),locy(1,selnode), 'Color', colpos(1,:));
    hold on;
    quiver(locx(selt,selnode),locy(1,selnode),loccdd(selt,selnode),loccll(1,selnode));
    hold on;
    h= plot(0*locx(selt,selnode),locy(1,selnode), 'k--');
    xlabel('x (mm)', 'Fontname',FONT,'FontSize',FONTSIZE);
    ylabel('y (mm)', 'Fontname',FONT,'FontSize',FONTSIZE);
    xlim([0 65])
    ylim([-100 100]);
    subplot(1,2,2)
    h=plot(loccdd(selt,selnode),locy(1,selnode), 'Color', colpos(1,:));
    hold on
    h= plot(0*locdisp(selt,selnode),locy(1,selnode), 'k--');
    hold on
    plot(0*locdisp(selt,selnode)+FF(selt)*10,locy(1,selnode), 'k-');
    xlabel('Drag (N/m)', 'Fontname',FONT,'FontSize',FONTSIZE);
    ylabel('y (mm)', 'Fontname',FONT,'FontSize',FONTSIZE);
    xlim([-125, 125])
end
legend('t/T = 0','t/T = \pi/4','t/T = \pi/2','t/T = 3\pi/4' );
LEG = findobj(AX,'type','axes');
set(LEG,'Fontname',FONT,'FontSize',FONTSIZE);
set(gca,'Fontname', FONT,'FontSize',FONTSIZE);
set(gcf, 'PaperUnits','inches', 'PaperPosition',[0 0 pWidth pHeight], ...
    'PaperSize', [pWidth pHeight]);

return
figure
for ii = 1:3
    subplot(3,1,ii)
    plot(timeall,locdisp(:,sampleselnode(ii))*10);
    hold on
    plot(timeall,loccdd(:,sampleselnode(ii)));
    hold on;
    h=plot(timeall,FF*10, 'Color', colpos(1,:));
    hold on;
    plot(timeall,locdisp(:,sampleselnode(ii))*0,'k--');
    xlabel('Time (s)')
    legend('Displacement x *10', 'Stress x', 'Actuation force x')
    
end

Power = FF(1:end-1).*vel;
Powerf = cdd(1:end-1,2)*Ltott.*vel;
Powertot = sum(Powerf)./sum(Power+Powerf)

return
for snapt=10:10:200
    plot(locx(snapt,:),locy(snapt,:),'-')
    hold on
    xlim([-40,40])
    ylim([-100,100])
end
return
FF=  2*2*(2*mod(floor(cdd(:,1)/Tcycle(casenum)),2)-1);
h1 = figure(1);

h=plot(cdd(:,1),FF, 'Color', colpos(1,:));
hold on;
h=plot(cdd(:,1),cdd(:,2)*Ltott, 'Color', colpos(2,:));
hold on;
h=plot(disp(:,1),(disp(:,2)-30)/10, 'Color', colpos(3,:) );

box on
xlabel('Time (s)', 'Fontname',FONT,'FontSize',FONTSIZE);
ylabel('Output', 'Fontname',FONT,'FontSize',FONTSIZE);
AX = legend('Actuation force (N)', 'Total horizontal hydrodynamic force (N)', 'Displacement (cm)');
LEG = findobj(AX,'type','text');
set(LEG,'Fontname',FONT,'FontSize',FONTSIZE);
set(gca,'Fontname', FONT,'FontSize',FONTSIZE);
set(gcf, 'PaperUnits','inches', 'PaperPosition',[0 0 pWidth pHeight], ...
    'PaperSize', [pWidth pHeight]);
saveas(h1, ['Forcedisp',textnum,'.pdf']);
vel = diff(disp(:,2))/1000;
figure
plot(disp(1:end-1,2),vel)
% return
h1=figure;
h=plot(disp(1:end-1,1),Power,'-',disp(1:end-1,1),Powerf,'--');
%     hold on;
TotPower = sum(Power)
totnum = length(vel);
dtt = .01;
maxvel = max(vel(round(totnum/3):round(totnum/3*2)))*1000/dtt;
maxdisp = max(disp(round(totnum/3):round(totnum/3*2),2));
minvel = min(vel(round(totnum/3):round(totnum/3*2)))*1000/dtt;
mindisp = min(disp(round(totnum/3):round(totnum/3*2),2));

approximateupt=(maxdisp-mindisp)/maxvel
approximatedownt=(maxdisp-mindisp)/abs(minvel)
% return
box on
xlabel('Time (s)', 'Fontname',FONT,'FontSize',FONTSIZE);
ylabel('Power (W)', 'Fontname',FONT,'FontSize',FONTSIZE);
AX = legend('Actuation force (N)', 'Total horizontal force (N)', 'Displacement (cm)');
LEG = findobj(AX,'type','text');
set(LEG,'Fontname',FONT,'FontSize',FONTSIZE);
set(gca,'Fontname', FONT,'FontSize',FONTSIZE);
set(gcf, 'PaperUnits','inches', 'PaperPosition',[0 0 pWidth pHeight], ...
    'PaperSize', [pWidth pHeight]);
saveas(h1, ['Power',textnum,'.pdf']);

figure
contourf(disp(:,1),1:20,locdisp')
figure
contourf(disp(:,1),1:20,loccdd')

if locc==1
    snapt=100;
    h1 = figure;
    h=plot(1:20,loccdd(snapt,:), 'Color', colpos(1,:),'o');
%     hold on;
    h=plot(locdisp(snapt,:),(locdisp(snapt,:)-30)/10, 'Color', colpos(3,:) );
%     box on
    xlabel('Time (s)', 'Fontname',FONT,'FontSize',FONTSIZE);
    ylabel('Output', 'Fontname',FONT,'FontSize',FONTSIZE);
    AX = legend('Actuation force (N)', 'Total horizontal hydrodynamic force (N)', 'Displacement (cm)');
    LEG = findobj(AX,'type','text');
    set(LEG,'Fontname',FONT,'FontSize',FONTSIZE);
    set(gca,'Fontname', FONT,'FontSize',FONTSIZE);
    set(gcf, 'PaperUnits','inches', 'PaperPosition',[0 0 pWidth pHeight], ...
        'PaperSize', [pWidth pHeight]);
    saveas(h1, ['Locdisp',textnum,'.pdf']);
end
return
for symm = 0;
    if symm ==1
        stringinp = {'05','1','2','5','10'};vers= 'v1';
    else
        stringinp = {'100','150','200','300','500'};vers= 'v2';
    end
    niiall = 1:length(stringinp);
    h1 =figure;
    for nii = niiall
        data = load([vers,stringinp{nii},'n.txt']);
        plot(data(:,1),data(:,2)-30,'-');
        totl = length(data(:,1));
        limtt(nii) = data(find(data(1:totl/2,2)== max(data(1:totl/2,2))),1);
%         limtt2(nii) = data(find(data(1:totl/2,2)== max(data(1:totl/2,2))),1);
        maxpos(nii)= data(find(data(1:totl/2,2)== max(data(1:totl/2,2))),2)/1000;
        maxvel(nii)=  maxpos(nii)/limtt(nii);
        Tcycle(nii) = (limtt(nii)-.5);
        hold on;
        FF=  str2num(stringinp{nii})*(2*mod(floor(data(:,1)/.5),2)-1);
        vel = diff(data(:,2))/1000;
        Power = FF(1:end-1).*vel;
    %     plot(data(1:end-1,1),Power,'o-');
    %     hold on;
        TotPower(nii) = sum(Power);
    %     ratp(nii) = str2num(stringinp{nii})*(limpos(nii)-.5)^2
      if nii==1 && symm ==1
          legendText{nii} = ['F_{max}= ', num2str(1), 'N'];
      else
          legendText{nii} = strcat('F_{max}= ', num2str(str2num(stringinp{nii})*2), 'N');
      end
      AX= legend(legendText,'Location','northwest');
%       LEG = findTobj(AX,'type','text');
%       set(LEG,'Fontname',FONT,'FontSize',FONTSIZE);
    hold all;
    end
maxvel
TotPower
end
hold on;
plot([-.5;data(:,1)],[0;data(:,1)*0+0],'k--');
if symm==0
ylim([-40 50])
else
ylim([-40 50])
end
xlim([-.2 2])

box on
xlabel('Time(s)', 'Fontname',FONT,'FontSize',FONTSIZE);
ylabel('Displacement (mm)', 'Fontname',FONT,'FontSize',FONTSIZE);
set(gca,'Fontname', FONT,'FontSize',FONTSIZE);
set(gcf, 'PaperUnits','inches', 'PaperPosition',[0 0 pWidth pHeight], ...
    'PaperSize', [pWidth pHeight]);
saveas(h1, ['Structbuckle',num2str(symm),'.pdf']);
return
h= 2
h^3*(30-h)^2

h= 20
h^3*(30-h)^2