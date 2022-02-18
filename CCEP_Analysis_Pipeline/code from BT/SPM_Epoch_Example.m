D = spm_eeg_load();

% Define the parameters
S.D      = D;
S.bc     = 0;
S.trl    = [1000 2500 -500;3000 4500 -500];
S.conditionlabels = {'Angry','Happy'};
 
D = spm_eeg_epochs(S);

D_Raw = spm_eeg_load()
D.trialonset


D = meeg(D)
D.time

sum(D_Raw(50,[1000:2500],1) - D(50,:,1))

sum(D_Raw(50,[3000:4500],1) - D(50,:,2))
