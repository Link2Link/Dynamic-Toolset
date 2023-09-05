function val = getfield_safe(opt,field,default_val)
%
% Get the field item if exists, and return default otherwise.
%

if isfield(opt,field)
    val = getfield(opt,field);
else
	fprintf(2,'[getfield_safe] [%s] field does not exist. Using [%s] instead. \n',...
        field,num2str(default_val));
    val = default_val;
end
