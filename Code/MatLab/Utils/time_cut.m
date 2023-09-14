function x = time_cut(x,options)

%Sanity check
if(options.timecut_interval(1)<0)
    warning('Time interval specified not correct, the time selection is disabled');
    idx_min=1;
    idx_max=length(x);
elseif(options.timecut_interval(2)>length(x)/options.fs)
    warning('Time interval specified not correct, the time selection is disabled');
    idx_min=1;
    idx_max=length(x);
elseif(options.timecut_interval(1) >= options.timecut_interval(2))
    warning('Time interval specified not correct, the time selection is disabled');
    idx_min=1;
    idx_max=length(x);
else
    if(options.timecut_interval(1)==0)
        idx_min=1;
    else
        idx_min=round(options.fs*options.timecut_interval(1));
    end
    idx_max=round(options.fs*options.timecut_interval(2));
end

x=x(idx_min:idx_max);

end%EOF
