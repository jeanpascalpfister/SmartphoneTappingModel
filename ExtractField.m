function list = ExtractField(Data,s)
% JPP, 26.8.2015
%
% Ex: alpha = ExtractField(Data,'fit(1)')
%
% makes a vector (or a matrix) of values from specified fields in Data

n = length(Data);

% extract subfields (e.g. if s= 'RT.RTT1' -> s1 = 'RT')
indsep = strfind(s,'.');
indcurly = strfind(s,'{');

if isempty(indsep)    
    s1 = s;
else
    if isempty(indcurly)
        s1 = s(1:indsep-1);
    else
        s1 = s(1:min(indcurly,indsep)-1);
    end
end

for k=1:n    
    if isfield(Data{k},s1)
        eval(['a = Data{k}.' s ';']) 
    else
        a = NaN;
    end
    if isnumeric(a)
        list(k,:) = a;
    elseif isstr(a)
        list{k} = a;
    else
        error('unknown format')
    end
    
   %eval(['list(k,:) = Data{k}.' s ';']) 
end


end

