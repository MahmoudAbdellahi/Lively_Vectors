function data = lv_2dto3d(data,sz)
    % goes from 2d(ch*time trials) to 3d(trial*ch*time) 
    data = reshape(data,sz(2),sz(3),[]); 
    data = permute(data,[3 1 2]);

end