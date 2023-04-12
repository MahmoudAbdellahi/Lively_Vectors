function [stat,p] = lv_mdim_pt_stats(x, parametric)
% calculating the sample wise stats, takes midim difference of the conditions and returns the zvalue and
% pvalue of the wilcoxon test or tstat and pvalue of a ttest according to parametric flag.. difference maybe in mdim and the z/p
% values will be in mdim .. first dim is ppnts/samples so it will be
% compressed when the stats are done..
 
if parametric~=1 
    dims = size(x);
    x = reshape(x, dims(1),[]); % to 2d
    % removing zeros and getting the average rank of ties
    x(x==0)=nan;
    [val,ids] = sort(abs(x),1);
    ranks = repmat((1:size(val,1))',1,size(x,2)); % ranks in: (1:length(val))' because the smallest is when abs(x) is closest to zero

    [tier,tiec] = find(diff(val,1)==0);
    tiec_temp = unique(tiec);
    for c=1:length(tiec_temp)% for each column when there is/are tie(s)
        tie = tier(tiec == tiec_temp(c)); % rows of the respective column, because a column could include many rows

        shifts = [1 (find(diff(tie)~=1)+1)' length(tie)+1]; % diff(tie) to see when they don't come after each other
        for i=1:length(shifts)-1
            id = [tie( shifts(i) : shifts(i+1)-1 ) ; tie(shifts(i+1)-1)+1];
            ranks(id,tiec_temp(c)) = mean( ranks(id,tiec_temp(c)) ); % average rank for ties
        end
    end
    n = sum(~isnan(x),1); % values that don't equal nans

    x(isnan(x))=0; ranks(isnan(x))=0; temp=[];
    
%     for i=1:size(x,2), temp(:,i) = sign(x(ids(:,i), i)); end
    % linear indexing instead of looping
    rows = ids(:)'; columns = reshape(repmat(1:size(ids,2),size(ids,1),1), 1,[]);
    inds = sub2ind(size(ids), rows, columns);
    temp = x(inds);
    temp = sign(reshape(temp, size(x)));

    T = sum(temp .* ranks, 1); 
    
    varT = (2*n.^3 + 3*n.^2 +n)*1/6;
    z = T ./ sqrt(varT); 
    if length(dims)>2 % because if 2d then the first dimesion is compressed 
        stat = reshape(z, dims(2:end) ); else, stat = z; end
    p = (1-normcdf(stat))*2; % two-tailed test
elseif parametric==1
    [~,p,~,stat] = ttest(x);
    stat = squeeze(stat.tstat);
    p = squeeze(p);
end

end