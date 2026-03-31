function [ S ] = Generate_STD( tpos, networkSize, spikeMatrix, tau, pfactor )

S = zeros(length(tpos), networkSize);

for j = 1:networkSize
    %Create vector of index
    spkd = find(spikeMatrix(:,j)==1);
    R0 = pfactor;
    if ~isempty(spkd)
        S(spkd(1),j) = R0;
        for i = 2:(length(spkd)+1)
            if i < (length(spkd) + 1)
                t0 = spkd(i-1)+1;
                t1 = spkd(i);
                t = (t0:t1)';
                if size(t,1)==5745
                    dddd=0;
                end
                S(t,j) = R0.*exp(-(tpos(t)-tpos(t0-1))./(tau)) ;
                S(spkd(i),j) = pfactor*(1-S(spkd(i), j)) + S(spkd(i), j);
                R0 = S(spkd(i),j);
            else
                t0 = spkd(i-1)+1;
                t1 = length(tpos);
                t = (t0:t1)';
                S(t,j) = R0.*exp(-(tpos(t)-tpos(t0-1))./(tau));
            end
        end
    end
end

S(2:end,:)=S(1:end-1,:);
s(1,:)=0;

end

