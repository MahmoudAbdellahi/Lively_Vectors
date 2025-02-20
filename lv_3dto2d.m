function data = lv_3dto2d(data,sz)
    % from 3d(trial*ch*time) to 2d(ch*time trials) 
    data = permute(data,[2 3 1]); 
    data = reshape(data,sz(2),[]); %2d
end