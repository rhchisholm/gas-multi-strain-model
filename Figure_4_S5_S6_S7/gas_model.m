function SS = gas_model(R0,DI,sigma,omega,x,NumMigrationEvents,burnin_time)

numsims = length(R0);

SumStats = zeros(numsims,12);

p = gcp('nocreate'); %If no pool, do not create new one.

if isempty(p)
    parpool(16)
end

parfor i = 1:numsims
    r = R0(i);
    d = DI(i);
    s = sigma(i);
    SumStats(i,:) = multistrain_model42(r,d,s,omega,x,NumMigrationEvents,burnin_time);
end

SS = SumStats;