function [bad_chan] = promptBadChanSpec


%% Block
bad_chan = nan;
while sum(isnan(bad_chan)) >= 1 && ~isempty(bad_chan)
    prompt = 'bad_chan: ';
    bad_chan_tmp = input(prompt,'s');
    if ~isempty(bad_chan_tmp)
        bad_chan = str2double(strsplit(bad_chan_tmp));
    else
        bad_chan = [];
    end
    if sum(isnan(bad_chan)) >= 1 && ~isempty(bad_chan)
        disp('You can only type numbers')
    else
    end
end


end