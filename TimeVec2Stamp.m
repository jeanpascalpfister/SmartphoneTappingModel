function t = TimeVec2Stamp(v)
% JPP 3.11.2016
% converts a date vectore [y m d h m s] into a time stamp.
% it is the opposite of Stamp2TimeVec


t = (datenum(v)-datenum(1970,1,1))*86400*1000;

end

