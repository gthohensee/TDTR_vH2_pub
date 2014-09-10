function [tdelay_data,dr_data,tdelay_deviation] = ROR_interp(tdelay_data,ratio_data,...
                                                             tdelay_data2,ratio_data2)
Lmax = max(length(tdelay_data),length(tdelay_data2));
t1 = zeros(Lmax,1); t2 = t1;
t1(1:length(tdelay_data)) = tdelay_data;
t2(1:length(tdelay_data2)) = tdelay_data2;

if sum(t1 == t2) ~= Lmax % FALSE if t1 == t2.
    fprintf('Paired data sets have different time delay sampling. Interpolating...\n')
    if length(tdelay_data) >= length(tdelay_data2)
        t1_interp = interp(tdelay_data,10);
        r1_interp = interp(ratio_data,10);

        t1_2 = -1*ones(length(tdelay_data2),1);
        r1_2 = t1_2;
        for tj = 1:length(tdelay_data2)
            [~,index] = min(abs(t1_interp - tdelay_data2(tj)));
            t1_2(tj) = t1_interp(index);
            r1_2(tj) = r1_interp(index);
        end

        dr_data = r1_2 ./ ratio_data2;
        tdelay_data = tdelay_data2;
        tdelay_deviation = sum(abs(t1_2 - tdelay_data2));
    else % reverse case
        t2_interp = interp(tdelay_data2,10);
        r2_interp = interp(ratio_data2,10);

        t2_1 = -1*ones(length(tdelay_data),1);
        r2_1 = t2_1;
        for tj = 1:length(tdelay_data)
            [~,index] = min(abs(t2_interp - tdelay_data(tj)));
            t2_1(tj) = t2_interp(index);
            r2_1(tj) = r2_interp(index);
        end

        dr_data = ratio_data ./ r2_1;
        %tdelay_data = tdelay_data;
        tdelay_deviation = sum(abs(t2_1 - tdelay_data));
    end 
else
    dr_data = ratio_data ./ ratio_data2;
    %tdelay_data = tdelay_data;
    tdelay_deviation = 0;
end

end

