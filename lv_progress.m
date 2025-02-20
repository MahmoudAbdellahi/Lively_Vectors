function lv_progress(iterator, progress_end, msg) % normalised_progress is between 0 and 1 
% EXAMPLE:  
% progress_end=100;
% for i=1:progress_end
%    lv_progress(i, progress_end, 'working on something: ')
%    pause(0.1);
% end
 
normalised_progress = iterator/progress_end;
no_squares = 20; 

if iterator==1, fprintf(['\n' msg ]); end

line_length = no_squares+8; % becuase the print of percentage: 2 spaces then 3 decimals and then the '%' symbol and the 2 brackets.. and no_squares is the rest of the line

if iterator>1, fprintf( repmat('\b',1,line_length) ); end % +8

current_squares = round(normalised_progress*no_squares);
fprintf(['[' repmat('|',1,current_squares) repmat('.',1,(no_squares-current_squares)) ']']);
 
fprintf('  %3.0f%%',normalised_progress*100); 

if iterator==progress_end, fprintf('\n'); end

end

 
