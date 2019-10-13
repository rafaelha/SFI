classdef cons
    properties (Constant=true)
        Avo=6.02214129e23;
        R=8.314462145468951;
        F=96485.3365;
        el=1.60217657e-19;
        epsilon=8.854187817e-12;
        epsilon2d=7.323564369075211e-18;
        kB=1.3806488e-23;
        mi=4.981733643307871e-26; % mass of NO in kg
        me=9.10938291e-31;
        mp=1.67262178e-27;
        h=6.62606957e-34;
        hbar=cons.h/(2*pi);
        G=6.67384e-11;
        c=299792458; %m/s
        a0=5.2917721092e-5;
        kBau=0.012374764324710;
        Ryd=10973731.6;
        Rydhc=2.179872000000000e-18;%2.181381270723222e-18;
        RydhcAU=1.955173475509261e+03;
        RydkB=1.579968251682269e+05;
        NOrot=1.67195; %cm-1
        NOprot=1.9971945;
        NOIP=30522.45;
        NOIPN2=30522.45+11.9;
    end
    methods (Static)
        function y=NN(r,den)
            y=4*pi*den*r.^2.*exp(-4*pi*den*r.^3/3);
        end
        function [pf, eden, rden]=penningfraction(n,den)
            Rn0=n.^2*cons.a0; 
            % radius of Rydb. by bohr model using semi-classical method
            Rmax=1.8*(Rn0*2); 
            % Robicheaux paper, within this distance, 90% penning ionize
            pf=1-exp(-4*pi*den*Rmax.^3/3); 
            % proportion between 0 and Rmax
            eden=pf/2*den; 
            % the den of electron produce is half the proportion 
            % (1e- per partner)
            rden=(1-pf)*den; 
            % this is remaining density of rydbergs
        end
        function [r,rp]=randonsphere
            x=(rand-.5)*2;
            y=(rand-.5)*2*sqrt(1-x.^2);
            z=sqrt(1-x.^2-y.^2);
            if rand>.5
                z=-z;
            end
            r=[x;y;z];
            r=r(randperm(3));
            
            x=(rand-.5)*2;
            y=(rand-.5)*2*sqrt(1-x.^2);
            z=sqrt(1-x.^2-y.^2);
            if rand>.5
                z=-z;
            end
            rp=[x;y;z];
            rp=rp(randperm(3));
            rp=rp-r*sum(rp.*r);
            rp=rp/sqrt(sum(rp.^2));
        end
        function [r,rp]=randoncircle
            % create a vector pointing randomly
            x=(rand-.5)*2;
            y=sqrt(1-x.^2);
            if rand>.5
                y=-y;
            end
            r=[x;y];
            r=r(randperm(2));
            % do it one more time
            x=(rand-.5)*2;
            y=sqrt(1-x.^2);
            if rand>.5
                y=-y;
            end
            rp=[x;y];
            rp=rp(randperm(2));
            
            % subtract the common part of the two vector to create normal
            % ones that are prependicular
            rp=rp-r*sum(rp.*r);
            % normalize the new prependicular vector
            rp=rp/sqrt(sum(rp.^2));
        end
        function GHz=cmtoGHz(cm)
            GHz=cons.c*cm*100;
        end
        function mm=GHztomm(GHz)
            mm=cons.c/(GHz*1e9)*1000;
        end
        function cm=nmtocm(nm)
            cm=1e7./nm;
        end
        function y=lambdanm(n1,n2)
            y=1e9./(10968800*abs(1./n1.^2-1./n2.^2));
        end
        function y=yukawa(r,l)
            y=cons.el.^2./(4*pi*cons.epsilon*r).*exp(-r/l);
        end
        function y=aws(den)
            y=(3./(4*pi*den)).^(1/3);
        end
        function y=a2d(den)
            y=(1./(4*pi*den)).^(1/2);
        end
        function y=awsden(aws)
            y=1./(4*pi*aws.^3/3);
        end
        function y=debye(Te,ne)
            y=sqrt(cons.epsilon*cons.kB.*Te./(ne*cons.el^2));
        end
        function y=debnum(Te,ne)
            y=4*pi*ne*cons.debye(Te,ne)^3/3;
        end
        function y=yukawaTene(r,Te,ne)
            y=cons.yukawa(r,cons.debye(Te,ne));
        end
        function y=ncritical(T)
            y=round(sqrt(cons.Rydhc/cons.kB./T));
        end
        function y=scaledT(Ti,Te,ne)
            y=cons.kB*4*pi*cons.epsilon*Ti*cons.debye(Te,ne)/cons.el^2;
        end
        function y=scaledn(ni,Te,ne)
            y=4*pi*ni*cons.debye(Te,ne).^3/3;
        end
        function y=g(den,Te)
            y=cons.el^2./(4*pi*cons.epsilon*cons.kB*Te*cons.aws(den));
        end
        function y=we(den)
            y=sqrt(den*cons.el^2/(cons.me*cons.epsilon));
        end
        function y=wpi(den)
            y=sqrt(den*cons.el^2/(cons.mi*cons.epsilon));
        end
        function y=av(r,N)
            y=sqrt(sum(r.^2.*N,2)./sum(N,2)/3);
        end
        function y=scaledtoEn(scEn,den)
            y=scEn*cons.me*cons.aws(den).^2.*cons.we(den).^2/cons.kB;
        end
        function y=scaledtot(sct,den)
            y=sct/cons.we(den);
        end
        function y=encm(n)
            y=cons.Ryd/100/n^2;
        end
        function y=dErot(J1,J2,B)
            y=B*(J1*(J1+1)-J2*(J2+1));
        end
        function y=EF(n,J,a)
            y=((cons.encm(n)-cons.dErot(2,J,2))/a)^2;
        end
        function n=boundton(En,den)
            n=(-cons.scaledtoEn(En,den)*cons.kB/cons.Rydhc).^(-1/2);
        end
        function str=gettimedate()
            str=datestr(now);
            str=strrep(str,':','_');
            str=strrep(str,'-','_');
            str=strrep(str,'.','_');
            str=strrep(str,' ','_');
        end
        function cc=fitgauss(t,data,tlim,tag)
            [a,~]=size(data);
            ind=tlim(1):tlim(2);
            t=t(ind)';t=t(:);
            data=data(:,ind);
            ft=fittype('gauss1');
            cc=zeros(a,3);
            
            parfor i=1:a
                y=data(i,:);y=y(:);
                fo(i)=fitoptions(ft);
                fo(i).Lower=[0 min(t) 0];
                ub=max(y)-min(y)+eps;
                fo(i).Upper=[ub*1.2 max(t) max(t)-min(t)];
                if strcmp(tag,'ramp')
                    [~,b]=findpeaks(y,'minpeakdistance',100,...
                    'minpeakheight',mean(y)/2,'sortstr','descend');
                    if length(b)>=1
                        b=b(1);
                        fo(i).StartPoint=[y(b) t(b) 1];
                    else
                        fo(i).StartPoint=[mean(y)/2 mean(t) 1];
                    end
                elseif strcmp(tag,'fixed')
                    fo(i).StartPoint=[mean(y)/2 mean(t) max(t)-min(t)];
                end
                cf=fit(t,y,ft,fo(i));
                cc(i,:)=coeffvalues(cf);
            end
        end
        function EF=ntoEF(n,tp)
            % n=(5.14E+9./(16*EF)).^(1/4) ->  EF=5.14E+9./(16*n.^4)
            if strcmp(tp,'a')
                EF=5.14E+9./(16*n.^4);
            elseif strcmp(tp,'d');
                EF=5.14E+9./(9*n.^4);
            end
        end
        function y=hgauss(x,cc)
            a=cc(1);
            b=cc(2);
            c=cc(3);
            y=a*exp(-(x-b).^2/c^2);
        end
        function tau = rydperiod(n)
             %\tau^2 = {4\pi^2\mu \over kZe^2}a^3 
             tau = sqrt(16*pi^3*cons.epsilon*cons.me/(cons.el^2)...
                 *(cons.a0*1e-6*n^2)^3);
        end
        function tau = rydperiodscaled(n,den)
             %\tau^2 = {4\pi^2\mu \over kZe^2}a^3 
             tau = cons.rydperiod(n)*cons.we(den);
        end
    end
end