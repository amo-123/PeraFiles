function [weights] = FcnWLSweightCalc(Frame,Par,Filt)

variance_LRF = zeros(size(Par.LRF,1),size(Par.LRF,2),size(Par.LRF,3));
for ind_ch = 1:36
    eta = permute(Par.LRF(ind_ch,:,:),[2,3,1]);
    variance = zeros(size(eta,1),size(eta,2));
    i_0 = 1;
    i_end = Filt.Num_rec;
    for i = i_0:i_end
    diff = sum(Frame(i,:))*eta-Frame(i,1);
    variance = variance + diff.^2;
    end
    variance = variance./(i_end-i_0);
%     figure
%     imagesc(variance)
%     drawnow
    variance_LRF(ind_ch,:,:) =  variance;
end
weights = 1./variance_LRF;

end