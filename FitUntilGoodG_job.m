function [k,s,q,h,G,D] = FitUntilGoodG_job(goodG, modelgenfun, varargin)

G = Inf;

% Repeatedly fit to stagnation until we get a good enough fit
while G > goodG
    % Fit until stagnation
    [G,k,s,q,h] = DetermineGoodG_job(modelgenfun, varargin{:});
end

% Continue the fit
[m,con,obj,opts] = modelgenfun();
m = m.Update(k);
for i = 1:numel(con)
    con(i) = con(i).Update(s(:,i), q{i}, h{i});
end
[m,con,G,D] = FitObjective(m,con,obj,opts);
k = m.k;
s = [con.s];
q = {con.q};
h = {con.h};

end