function [DeconvMat,FiltMat]=CaDeconvNew(inmat,tau,thrup,thrdc,pfilt,StrongSm,maxcount,tbin,butterfilt,plotit);
% EXAMPLE THAT WORK FOR TRONDHEIM DATA:[DeconvMat,FiltMat]=CaDeconvNew(c,3,30,0,100,0,200,1/31,0.1,0);

% inmat: input data: time x cells x stimuli 
% EWELINA for your data you will not have stimuli, so you need to enter a 2D matrix, and first dimension is time, SO YOU MAY NEED TO ROATTE YOUR MATRIX
% tau: exponential time constant of kernel for deconv
% thrup: threshold for filtering, appropximately noise level, in % DF/F,
% normally between 0.5 and 5
% thrdc: set to 0
% pfilt: number of times for post-filtering, ie, filtering of whole trace
% after iterative filtering of segments around peaks. Normally 0.
% StrongSm: if 1, a slightly different smoothing algorkithm is used. This
% causes more distortions. Normally, set to 0.
% maxcount: maximum number of iterations in filter procedure. Normally 200.
% tbin: frame time of input data in sec.
% butterfilt: cutoff frequency for butterworth filter. For tbin = 128 ms,
% set to 0.4.  When tbin = 0.256, set to 0.2, and so on.
% plotit: if 1, program outputs some plots for control purposes.
ntraces=size(inmat,2);
nodors=size(inmat,3);

FiltMat=inmat;

for a0=1:nodors,
    disp(['Processing odor ',int2str(a0)]);
    Dmat=inmat(:,:,a0);
    if plotit,
        figure(5);clf;
        plot(Dmat); hold on;
    end
    
    if butterfilt>0,
        [B,A]=butter(4,butterfilt);
        for a1=1:ntraces,
            Dmat(:,a1)=filtfilt(B,A,Dmat(:,a1));
            if plotit,
                plot(Dmat,'g-');
            end
        end
    end
    FiltMat(:,:,a0)=Dmat;
    
    if thrup>0,
        disp('smoothing...')
        for a1=1:ntraces,
            changed=1;
            count=1;
            lastpk=0;
            
            while ( (changed==1) & (count<maxcount) ),
                [pospeaks,negpeaks]=peakdetect(Dmat(:,a1));
                if ( (length(pospeaks)>0) & (length(negpeaks)>0) )
                    if (pospeaks(1)<negpeaks(1))
                        peaks1=negpeaks;
                        peaks2=pospeaks;
                    else
                        peaks1=pospeaks;
                        peaks2=negpeaks;
                    end
                else
                    changed=0;
                end
                posindex=ones(size(pospeaks));
                negindex=ones(size(negpeaks));
                negindex=negindex*(-1);
                allindex=[posindex;negindex];
                allpeaks=[pospeaks;negpeaks];
                if (length(allpeaks)==0)
                    allpeaks=[1;allpeaks];
                    if (Dmat(2,a1)<Dmat(1,a1)),
                        allindex=[1;allindex];
                    else
                        allindex=[-1;allindex];
                    end
                end
                
                if (allpeaks(1)>1),
                    allpeaks=[1;allpeaks];
                    if (Dmat(2,a1)<Dmat(1,a1)),
                        allindex=[1;allindex];
                    else
                        allindex=[-1;allindex];
                    end
                end
                if (allpeaks(length(allpeaks))<size(Dmat,1)),
                    allpeaks=[allpeaks;size(Dmat,1)];
                    if (Dmat(size(Dmat,1),a1)>Dmat(size(Dmat,1)-1,a1)),
                        allindex=[allindex;1];
                    else
                        allindex=[allindex;-1];
                    end
                end
                sort(allpeaks,1);
                %sort(allindex,1);
                
                allamps=zeros(size(allpeaks));
                for b1=2:(length(allamps)-1),
                    allamps(b1)=( abs(Dmat(allpeaks(b1),a1) - Dmat(allpeaks(b1-1),a1) ) + abs(Dmat(allpeaks(b1),a1) - Dmat(allpeaks(b1+1),a1) ) ) /2;
                end
                allamps(1)=abs(Dmat(1,a1)-Dmat(allpeaks(2),a1));
                allamps(length(allpeaks))=abs(Dmat(size(Dmat,1))-Dmat(allpeaks(length(allpeaks)-1)));
                %allamps(1)=min(allamps(2:(length(allamps)-1)));
                %allamps(length(allamps))=allamps(1);
                %allamps(1)=0;
                %allamps(length(allamps))=0;
                
                allamps2=allamps;
                [Y,allamps2ind]=sort(allamps2);
                allamps2=Y;
                
                count2=1;
                help3=0;
                take=0;
                %[allpeaks,allamps]
                
                while ( (count2<length(allamps2)) & (take==0) );
                    help4=allamps2ind(count2);
                    
                    if StrongSm==0,
                        if ( (allamps(help4)<thrup) & (lastpk ~= allpeaks(help4)) & (allamps(help4) ~=0) ),
                            take=1;
                        end
                    else
                        if ( (allamps(help4)<thrup) & (lastpk ~= allpeaks(help4)) & (allindex(help4)==1) ),
                            take=1;
                        elseif (allindex(help4)==(-1))
                            take=1;
                        end
                    end
                    count2=count2+1;
                end %2nd while loop
                if (take==1),
                    if (help4>1),
                        if help4>length(allpeaks),
                            stretch=Dmat(allpeaks(help4-1):allpeaks(help4+1),a1);
                        else
                            stretch=Dmat(allpeaks(help4-1):allpeaks(help4),a1);
                        end
                    else
                        stretch=Dmat(allpeaks(help4):allpeaks(help4+1),a1);
                    end
                    if length(stretch)>0,
                        stretch=itersmooth2(stretch,3);
                    end
                    if (help4>1),
                        if help4>length(allpeaks),
                            Dmat(allpeaks(help4-1):allpeaks(help4+1),a1)=stretch;
                        else
                            Dmat(allpeaks(help4-1):allpeaks(help4),a1)=stretch;
                        end
                    else
                        Dmat(allpeaks(help4):allpeaks(help4+1),a1)=stretch;
                    end
                    
                    
                    
                    lastpk=allpeaks(help4);
                    changed=1;
                else
                    changed=0;
                end
                count=count+1;
            end % 1st while loop
            if (pfilt>0),
                Dmat(:,a1)=itersmooth2(Dmat(:,a1),pfilt);
            end
        end %for a1 loop
    end
    if plotit,
     plot(Dmat,'c-','linewidth',2);
    end
    if tau>0,
        Dmat=deconvtrace2(tau,tbin,Dmat,0);
    end
    inmat(:,:,a0)=Dmat;
end %for a0 loop
DeconvMat=inmat;
if plotit,
    figure(5);
    if tau>0,   
        plot(2*Dmat,'r-','linewidth',2);
    else
        plot(Dmat,'r-','linewidth',2);
    end
end

% %butterfilt=1;
% butterfilt2=0;
% stepfilt=1;
% interfactor=1;
% stepfthrup=stepup;
% stepfthrdwn=10;
% tbin=0.128;
% %stimonfr=16;
% stepiter=5;
% boxlength=3;
% 
% intup=3;
% intdown=round(tconst/tbin);
% if intdown<1,
%     intdown=1;
% end
% 
% 
% [B,A]=butter(4,[0.2]);
% 
% nodors=size(Dmat,2)
% for a0=1:nodors,
%     disp(['Processing odor # ',int2str(a0)]);
%     Dmat1odor=squeeze(Dmat(:,a0,:))';
% 
%     Dmat1odorMax=max(Dmat1odor,[],1);
%     stretch=find( (Dmat1odorMax>thr1) & (Dmat1odorMax<thr2));
%     Dmat1odor2=Dmat1odor(:,stretch);
%     nselected=size(Dmat1odor2,2);
%     if a0==1,
%         DeconvMat=zeros(size(Dmat1odor2,1),size(Dmat1odor2,2),nodors);
%         FiltMat=zeros(size(Dmat1odor2,1),size(Dmat1odor2,2),nodors);
%     end    
%     time=((1:size(Dmat1odor2,1))*tbin)';
%     time=time*ones(1,size(Dmat1odor,2));
%     
%     if plotit,
%         figure(3); clf;
%         subplot(321);
%         plot(time,Dmat1odor)
%         title('All raw','fontsize',8);
%         v=axis;
%         if (min(time(:,1)) < max(time(:,1))),
%             axis([min(time(:,1)) max(time(:,1)) v(3) v(4)]);
%         end
%         line([stimonfr*tbin,stimonfr*tbin],[v(3),v(4)],'color',[0.7 0.7 0.7]);
%     end
%     
%     time=((1:size(Dmat1odor2,1))*tbin)';
%     time=time*ones(1,size(Dmat1odor2,2));
%     
%     if plotit,
%         figure(3);
%         subplot(322);
%         plot(time,Dmat1odor2)
%         title('Selected raw','fontsize',8);
%         v=axis;
%         if (min(time(:,1)) < max(time(:,1)))
%             axis([min(time(:,1)) max(time(:,1)) v(3) v(4)]);
%         end
%         line([stimonfr*tbin,stimonfr*tbin],[v(3),v(4)],'color',[0.7 0.7 0.7]);
%     end
%     
%     Dmat1odor3=zeros(size(Dmat1odor2,1)*interfactor,size(Dmat1odor2,2));
%     %nselected=size(Dmat1odor3,2);
%     for a1=1:nselected,
%         if butterfilt,
%             if (interfactor > 1)
%                 Dmat1odor3(:,a1)=filtfilt(B, A, interp(Dmat1odor2(:,a1),interfactor));
%             else
%                 Dmat1odor3(:,a1)=filtfilt(B, A, Dmat1odor2(:,a1));
%             end
%         else
%             if (interfactor > 1)
%                 Dmat1odor3(:,a1)=interp(Dmat1odor2(:,a1),interfactor);
%             else
%                 Dmat1odor3(:,a1)=Dmat1odor2(:,a1);
%             end
%         end
%     end
%     
%     if plotit,
%         figure(3);
%         subplot(323);
%         plot(time,Dmat1odor3);
%         title('Butter filtered','fontsize',8);
%         v=axis;
%         if (min(time(:,1)) < max(time(:,1)))
%             axis([min(time(:,1)) max(time(:,1)) v(3) v(4)]);
%         end
%         line([stimonfr*tbin,stimonfr*tbin],[v(3),v(4)],'color',[0.7 0.7 0.7]);
%     end
%     
%     Dmat1odor2butter=Dmat1odor3;
%     for a1=1:nselected,
%         if stepfilt,
%             if butterfilt2,
%                 Dmat1odor3(:,a1)=filtfilt(B, A, stepfilter3(Dmat1odor3(:,a1),stepfthrup,stepfthrdwn,stimonfr,intup,intdown,plotit));
%             else
%                 Dmat1odor3(:,a1)=stepfilter3(Dmat1odor3(:,a1),stepfthrup,stepfthrdwn,stimonfr,intup,intdown,plotit);
%             end
%         elseif butterfilt2,
%             Dmat1odor3(:,a1)=filtfilt(B, A, Dmat1odor3(:,a1));
%         end
%         
%         %Dmat1odor3(:,a1)=stepfilter(Dmat1odor3(:,a1),stepfthr,15);
%     end
%     
%     if plotit,
%         figure(3);
%         subplot(324);
%         plot(time,Dmat1odor3);
%         title('Stepfiltered','fontsize',8);
%         v=axis;
%         if (min(time(:,1)) < max(time(:,1)))
%             axis([min(time(:,1)) max(time(:,1)) v(3) v(4)]);
%         end
%         line([stimonfr*tbin,stimonfr*tbin],[v(3),v(4)],'color',[0.7 0.7 0.7]);
%     end
%         
%     Dmat1odor4=deconvtrace2(tconst,tbin,Dmat1odor3,plotit);
%     
%     if plotit,
%         figure(3);
%         subplot(326);
%         plot(time,Dmat1odor4);
%         title('Deconvolved','fontsize',8);
%         v=axis;
%         if (min(time(:,1)) < max(time(:,1)))
%             axis([min(time(:,1)) max(time(:,1)) v(3) v(4)]);
%         end
%         line([stimonfr*tbin,stimonfr*tbin],[v(3),v(4)],'color',[0.7 0.7 0.7]);
%         
%         figure(4);clf;
%         i1=ceil(sqrt(nselected));
%         i2=i1;
%         for a1=1:nselected,
%             subplot(i1,i2,a1),
%             plot(time,Dmat1odor2(:,a1),'g-');hold on;
%             plot(time,Dmat1odor2butter(:,a1),'c-');hold on;
%             plot(time,Dmat1odor4(:,a1)/max(Dmat1odor4(:,a1))*max(Dmat1odor3(:,a1)),'r-');hold on;
%             plot(time,Dmat1odor3(:,a1),'b-');hold on;
%             axis tight;
%             v=axis;
%             if (min(time(:,1)) < max(time(:,1)))
%                 axis([min(time(:,1)) max(time(:,1)) v(3) v(4)]);
%             end
%             axis([0 6 v(3) v(4)]);
%             line([stimonfr*tbin,stimonfr*tbin],[v(3),v(4)],'color',[0.7 0.7 0.7]);
%             
%         end
%     end
%     DeconvMat(:,:,a0)=Dmat1odor4;
%     FiltMat(:,:,a0)=Dmat1odor3;
% end




