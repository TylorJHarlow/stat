function [MI, MeanAmp]=MI(Pha, Am)

    % Estimates Modulation Index (Tort 2010)between phase of low frequency 
    % signal and the amplitude of a higher frequency signal in a manner similar
    % to the Event-Related PAC of Knight formulation.  
    % 
    % dimensions = [time x trials] if event related
    % dimensions = [trials x time] if global
    
    nbin = 36;
    position=zeros(1,nbin);
    winsize = 2*pi/nbin;
     
    % now we compute the mean amplitude in each phase:
    for i = 1:length(Pha(:,1))
        Phase = Pha(i,:);
        Amp = Am(i,:);
        MeanAmp=zeros(1,nbin); 
        for j = 1:(nbin)   
            position(j) = -pi + (j-1)*winsize;
            I = find(Phase <  position(j) + winsize & Phase >=  position(j));
            MeanAmp(i,j) = nanmean(Amp(I)); 
        end
        
        % the center of each bin (for plotting purposes) is position+winsize/2
         
        % quantifying the amount of amp modulation by means of a
        % normalized entropy index (Tort et al PNAS 2008):
        MI(i) = (log(nbin)-(-sum((MeanAmp(i,:)/sum(MeanAmp(i,:))).*log((MeanAmp(i,:)/sum(MeanAmp(i,:)))))))/log(nbin);
        
    end
end
