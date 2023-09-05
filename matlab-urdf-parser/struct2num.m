function out = struct2num(s)
%
% Structure to number 
% 

F = fieldnames(s);
nsub = length(F);
out = s;
for i=1:nsub
    % if the subitem is struct, do a recursive operation
    if(isstruct(s.(F{i})))
        out.(F{i}) = struct2num(out.(F{i})); % recursive
    else
        out.(F{i}) = str2num(out.(F{i}));
        % if string put it back to string
        if isempty( out.(F{i}) )
            out.(F{i}) = s.(F{i});
        end
    end
end
