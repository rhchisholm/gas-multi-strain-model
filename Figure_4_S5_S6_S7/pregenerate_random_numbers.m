%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Random numbers %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Generate random numbers for contact events
ContactRand = poissrnd(Cpertimestep,1e6,1);
countCR = 1;

% Generate random numbers for sampling contacts
SamplingContactsRand = rand(1e6,1);
countSCR = 1;

% Generate random numbers for migration events
MRRand = poissrnd(MRpertimestep,1e6,1);
countMR = 1;
