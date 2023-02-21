function dat = lv_load(path, varargin)
% uses h5 format to load the data seperated data.trial other_data from two
% files the h5 and .mat
% TO USE:
% takes optional argument with the name of the matrix that has double
% because that's the field that has the data and is the biggest so when
% you put its name it will be loaded from h5 and the other parts in other data
% if not specified then the normal load of matlab is used.. 


if ~isempty(varargin)
    array_name = varargin{1};
    
    load(path); % other_data
    dat=other_data; clear other_data; 
    eval(['dat.' array_name '= h5read([path ''.h5''] , [''/'' array_name]);'])  % data.trial .. tricky becaus ewhat's inside the eval and then inside h5read are variables that 
    % h5 read will know they are variables so eval just puts '' to execute
    % h5 and then h5read considers anything inside it as variables not
    % strings because it's a fucntion within a function(eval)
else
    dat = load(path);
end
    
    
end

