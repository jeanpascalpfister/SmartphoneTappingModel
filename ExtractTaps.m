function taps = ExtractTaps(D,v1,v2)
% JPP  3.11.2016
% extracts taps of subject k (D = Data{k}) from date vector v1 to v2

t0 = D.t0;
ITI = D.ITI;

taps = [t0 ; t0 + cumsum(ITI)];
t1 = TimeVec2Stamp(v1);
t2 = TimeVec2Stamp(v2);

taps = taps(t1<taps & taps<t2);


end

