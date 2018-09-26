

clc
clear all; 
close all;
z=zeros();
H=zeros();
v=zeros();


hl2=line(2,3);

wl=input('What is you use He-Ne(1) or SLD(2)? : ');
if (wl==2)   
   wavelength=848.80e-9;                                       %SLD
disp('You use SLD');
elseif (wl==1)
   wavelength=632.991432e-9;
   % wavelength=632.8e-9;                                         %He-Ne
disp('You use He-Ne Lasor');
end

   tra=wavelength/(4*pi);




ImRgb=imread('HeNe2_g2_2_0409_barcode2_1.jpg');
%ImRgb=imread('20.png','png');   %Input image0

figure(1), imshow(ImRgb), axis on;

alt=input('Do you want to alternate figure? YER(1)/NO(0) : ');

while alt==1
    figure(1), imshow(ImRgb), axis on;
 disp('What do you want to farmat?');
    farm=input('Alternate UP/DOWN(1) RIGHT/LEFT(2) BOTH(3) : ');
    if (farm==1)
% alternate up/down
ImRgb=flipud(ImRgb);
figure(1), imshow(ImRgb), axis on;
    elseif (farm==2)
% alternate right/left
ImRgb=fliplr(ImRgb);
figure(1), imshow(ImRgb), axis on;
    elseif (farm==3)
        ImRgb=flipud(ImRgb);
        ImRgb=fliplr(ImRgb);
figure(1), imshow(ImRgb), axis on;
    end
    alt=input('Do you ok? YER(0) NO(1)');
end

ImRgb = imgaussfilt(ImRgb,2);

J=rgb2gray(ImRgb);  
ImRgb=J;
J=ImRgb;

xy=0;
xxyy=0;

 while xy==0
    
x1=1;
x2=640;
y1=1;
y2=480;

if xxyy==1
    delete(rec);
end
rec=rectangle('Position',[x1 y1 x2-x1 y2-y1],'EdgeColor','r');

J=ImRgb([y1:y2], [x1:x2]);
figure(2), imshow(J),axis on;;

%xy=input('Are you suffice to figuer? : YER(1) NO(0) : ');
xy=1;
xxyy=1;
 end
%++++++++++++++++++++++++++++++++++++++++++
%ImGray=rgb2gray(ImRgb);
ImGray=ImRgb;

% figure(1), imshow(J);

   
ImGray=medfilt2(J);                                          %median 
%filter ???? noise ??????
DispGray=ImGray;                                             %????????????????????


DispGray=ImGray;
J=ImGray;



w=size(ImGray,2);
ww=size(ImGray,1);


graphrange=[10 size(J,2)];
 

beg=1;
las=size(J,2);


% ??????? scan ??????????

fprintf('+++++++++++++++++++++++++++++++++++++++++++++++++++++++\n');
fprintf('Your figure has size (Colume)X: %d ,(Row)Y: %d \n',size(J,2),size(J,1));
     %  disp('How much do you want to scan rows range (Y) of figure?');
     %  ini=input('ROW NUMBLE ini : ');
     %  fin=input('ROW NUMBLE fin : ');
     ini=1;
     fin=size(J,1);
rn=fin-ini+1;
% ??????????????

n=input('How many do you use average row valume? : ');
%fprintf('You are using %d row for create data\n',n);
%n=4;                                                         %??? n ???????????
%?????????????????????? 
if n==1                                                      %???????????????????=1
   n_Row = 0;                                                %??????????=?????????????
   ininew = ini;
   finnew = fin;
elseif mod(n,2)>0                                            %???????????????????????????????? ????????????????
       n_Row = floor(n/2); 
       ininew = ini+n_Row;  
       finnew = fin-n_Row;
       n_Row = -n_Row:n_Row; 
    elseif mod(n,2)==0                                       %???????????????????????????????? ????????????????
           n_Row = floor(n/2); 
           ininew = ini+n_Row-1;
           finnew = fin-n_Row;
           n_Row = -(n_Row-1):n_Row;
        else 
end                                                          %Return ???????? n_Row ?????????


A = ones(size(n_Row));                                       %??? A=matrix 1????=n_Row


ad = las-beg+1;
wm=0;
FW1=0;
FW2=0;

bkk=0;

back=0;
ex=0;
Wp1=0;
Wp2=0;


Ws1=0;
Ws2=0;
Ws=0;
Wp=0;
FW1=0;
FW2=0;
ALG=0;
reg=1;
for u = ininew:finnew
  
      
    close all;
    ul = u-ini+1;
    fprintf('NumBer %d',ul);
    clear B;
    B = (A*u)+(n_Row);
    
    %???????????????
   %figure(1), imagesc(ImRgb);                               %?????????????????? original(figure1)
    % figure(2), imshow(ImGray);                              %???????????????????????? gray????(figure2)
    
    DispGray(B,:)=255;                                       %???????? D ???????????????????????? (255=?????)
    
    %??????????????????????1
   %figure(3), imshow(DispGray),axis on;                             %??????????????????????????????gray(figure3)
%hAx = axes();
   
   % x = round( axes2pix(sz(2), [1 sz(2)], p(1)) );
    %    y = round( axes2pix(sz(1), [1 sz(1)], p(2)) );
   
    
    %************************ Row white line *************
  % yL = get(gca,'YLim');   line([beg beg],yL,'Color','r','LineWidth',0.5),...
%line([las las],yL,'Color','r','LineWidth',0.5); %vertical line
%********************************************************

   
   clear Row;
    Row = ImGray(B,:);                                       
    Row = double(Row);                                       
    SizeRow = size(Row);                                    
    SumRow = zeros(1,SizeRow(2));
   
    if (n>1)
        SumRow=sum(Row);
        Row = SumRow/n;
    end
         
                                           %????? 1 ??????????????????????????
    
    %????
    figure(4), plot(Row);                                    %??????????????????
    ylabel('Intensity');
    xlabel('x(pixel)');

    FT = fft(Row);                                         
    FTS = fftshift(FT);                                       
    FTS = FTS/max(abs(FTS));                               
    absFTS = abs(FTS);                                      

    clear ('Wp2','Wp1','Ws1','Ws2','Ws','Wp','FW1','FW2','cI','cI2');
    AF=abs(FTS);
%figure(1),plot(AF);
[rF,cF]=find(AF==max(AF));
    
    % Find peak interference
for cc=cF+20:w-100
    cAF1=AF(1,cc);
    cAF2=AF(1,cc+1);
    bxfor=1;                       % condition for non break for loop (if bxfor=2 will break)
    while (cAF1<cAF2)
    [rA,cA]=find(AF==cAF1);
    ITF=max(AF(1,[cA(1,2):w]));
    fprintf('cAF1= %f , cAF2= %f, ITF= %f\n',cAF1,cAF2,ITF);
    [rI,cI]=find(AF==ITF);
    % cI(1,2) is peak interference 
    fprintf('Peak inter = %d\n',cI(1,2));
 
    bxfor=2;
    break;
    end      
    if (bxfor==2)
        break;
    end    
end
    

% Generate Filter++++++++++++++++++

cI2=cI(1,2);





 while ALG==0;
    AL=input ('Your filter determine by Standard(1) or Algorithm(0) : ');
    if (AL==0)
     Wp1=posWp1;
     Wp2=posWp2;
     Ws1=Wp1-1;
     Ws2=Wp2+1;
     Rp=.001;
     Rs=.5;
    elseif(AL==1)
       % WideWp=input('Wide of Wp = ');
       % WideWs=input('Wide of Ws = '); 
       % Rp=input('Rp = ');
       % Rs=input('Rs = ');
      
        WideWp=3;
        WideWs=2; 
        Rp=.001;
        Rs=.5;
        
        
        
    end
    ALG=1;
 end
        Wp2=cI2+WideWp;
        Wp1=cI2-WideWp;
        Ws1=Wp1-WideWs;
        Ws2=Wp2+WideWs;
Ws=[Ws1 Ws2]/w;
Wp=[Wp1 Wp2]/w;

FW1=Ws1;
FW2=Ws2;


    [n,Ws] = cheb1ord(Wp,Ws,Rp,Rs);
    [d,c]=cheby1(n,Rp,Wp,'bandpass');
    [x,t] = freqz(d,c,length(Row),2*length(Row));            %filter signal in frequency domain
    Filter = (x);                                            %???????????????????????
    
    
    
    %*************************** filter line ******************************
 %  ax1 = gca;
 %    set(ax1,'XColor','k','YColor','k','xlim',(graphrange))
     
 %    ax2 = axes ('Position',get(ax1,'Position'),...
 %           'xlim',(graphrange),...
 %           'XAxisLocation','top',...
 %           'YAxisLocation','right',...
 %           'Color','none',...
 %           'XColor','k','YColor','k');
 %   hl2 = line(t,abs(Filter),'Color','r','Parent',ax2);

% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    
    Cut11=FTS'.*Filter;                                        
    [ back,Cut13,Cut12] = function_b(FW1,FW2,Cut11,w ); 
    figure(10), plot(abs(Cut13));                            
    Inv=ifft(Cut13);                                         
    figure(11), plot(abs(Inv));    
    Un=unwrap(angle(Inv));                                      
    tc=(beg:1:las);
    tt=tc';
    l = zeros(ad,1);
    l(:,1)=Un(beg:1:las,1);
    
   
     while (reg==1)
         
         linear=input('Linear Regrassion : Straight(1) or Curve(0) : ');
        reg=0; 
     
     end
     
     
     
    %?????????? 15
    %figure(15), plot(tt,l);
    
    
    
    
    if (linear==1)
    P = polyfit(tt,l,1);
    FitUn=P(1)*tt+P(2);
    %hold on;
    %?????????? 16 
    %figure(16), plot(tt,FitUn,'r');
    
    
    d=zeros();
    de=zeros();
    
          for mm=1:length(tt)
              d(mm)=l(mm)-FitUn(mm);
              de(mm)=tra*d(mm);
              z(ul,mm)=d(mm);
              H(ul,mm)=tra*z(ul,mm);               %hight profile base - to +
          end
    v=H-min(min(H));                             %hight profile base 0 to +a
    
    
    elseif (linear==0)
    fitpoly2=fit(tt,l,'poly2');
    %plot(fitpoly2,tt,l);
    %title('Wrapping')
    %hold on;
    %grid on;
    for mm=1:length(tt)
              d(mm)=l(mm)-fitpoly2(mm);
              de(mm)=tra*d(mm);
              z(ul,mm)=d(mm);
              H(ul,mm)=tra*z(ul,mm);               %hight profile base - to +
          end
    
    end
    %figure(21), plot(de);
    %grid on;
    suml(u,:)=l;
    Sum(u,:)=de;
% break;   
end

std2(H)                                          %finding roughness

figure(22), mesh(suml), axis on;
figure(20), mesh(Sum), axis on;
figure(3), imshow(J);
axis on;
save('Surface3Dprofile01.mat','suml','Sum','tra','w','ww');