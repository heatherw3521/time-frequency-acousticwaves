function sol = solve_HH_to_wave(sol,zzk)

tt = sol.times;
eval_type_time = sol.eval_type_time;
field_type = sol.field_type;
comp_type_freq = sol.comp_type_freq;
freq = sol.freq;
wts = sol.wts;
W = sol.bandlimit;

Uinc_pre = sol.incf;
if isa(Uinc_pre,'cell')
    Uinc = Uinc_pre;
else
    Uinc = {Uinc_pre};
end

uinc = sol.incf;
Uinc_slow = sol.incF_slow;
windows = sol.windows;
Usf = sol.freqvals;

u = zeros(max(size(zzk)),length(tt));
time2evalIFT = 0;

for i1 = 1:length(Uinc)
    switch eval_type_time
        case 'freq only'
            sol.timeevaltime = [];
            sol.timevals = [];
            sol.times = [];
        case 'complexify'
            tim = tic;
            switch comp_type_freq{1}
                case 'window one'
                    
                    for i2 = 1:length(Uinc_slow)
                    if size(zzk,1) == 3
                        Uinc_sl1 = Uinc_slow{i2}(zzk,freq{1});
                        Uinc_sl2 = Uinc_slow{i2}(zzk,freq{2});
                        Uinc_sl3 = Uinc_slow{i2}(zzk,freq{3});
                    else
                        Uinc_sl1 = Uinc_slow{i2}(zzk,freq{1}');
                        Uinc_sl2 = Uinc_slow{i2}(zzk,freq{2}');
                        Uinc_sl3 = Uinc_slow{i2}(zzk,freq{3}');
                    end
                        u1 = eval_sol_time(freq{1},W,zzk, Usf{1}{1}.* Uinc_sl1, tt - windows(i2),...
                            'GLquad', [], field_type, wts, uinc{i2});
                        u3 = eval_sol_time(freq{3},W,zzk, Usf{1}{3}.* Uinc_sl3, tt - windows(i2),...
                            'GLquad', [], field_type, wts, uinc{i2});
                        %u3 = -u3;
                        delt = sol.delt;
                        u2 = eval_sol_time(freq{2},W,zzk, Usf{1}{2}.* Uinc_sl2, tt - windows(i2),...
                            'complexify', delt, field_type, [], uinc{i2});
                        %u2 = -u2;
                    end
                otherwise
                    u1 = eval_sol_time(freq{1},W,zzk, Usf{i1}{1}, tt- windows(i1),...
                        'GLquad', [], field_type, wts, uinc);
                    u3 = eval_sol_time(freq{3},W,zzk, Usf{i1}{3}, tt- windows(i1),...
                        'GLquad', [], field_type, wts, uinc);
                    %u3 = -u3;
                    delt = sol.delt;
                    u2 = eval_sol_time(freq{2},W,zzk, Usf{i1}{2}, tt- windows(i1),...
                        'complexify', delt, field_type, [], uinc);
                    %u2 = -u2;
            end
            u = u + (u1*1i - u2 - u3*1i);

            timout = toc(tim);
            %u = eval_sol_time_faster(sol,zzk, Usf, tt,...
            %    'complexify', del, field_type, wts, uinc);
            time2evalIFT = timout;
        otherwise
            tim = tic;
            switch comp_type_freq{1}
                case 'window one'
                    for i2 = 1:length(Uinc_slow)
                        Uinc_sl = zeros(size(Usf{1}));
                   if size(zzk,1) == 3
                        Uinc_sl = Uinc_slow{i2}(zzk,freq);
                    else
                        Uinc_sl = Uinc_slow{i2}(zzk,freq');
                   end
                        u = u + eval_sol_time(freq,W,zzk, Usf{1} .* Uinc_sl, tt - windows(i2),...
                            eval_type_time, [], field_type, wts, uinc{i2});
                    end
                otherwise
                    u = u + eval_sol_time(freq,W,zzk, Usf{i1}, tt - windows(i1),...
                        eval_type_time, [], field_type, wts, uinc{i1});
            end

            timout = toc(tim);
            time2evalIFT = time2evalIFT + timout;
    end
end
if strcmpi(field_type, 'total')
    uinc = uinc{:};
    if size(zzk,1) == 3
        for i1 = 1:length(tt)
            uinc_temp = uinc(zzk,tt(i1)); % Needs to be changed later, lazy fix.
            if size(uinc_temp,1) == length(u(:,i1))
                u(:,i1) = u(:,i1) + uinc(zzk,tt(i1)); 
            else
                u(:,i1) = u(:,i1) + transpose(uinc(zzk,tt(i1))); 
            end
        end
    else
        u = u + uinc(zzk,tt.'); 
    end
end

sol.timevals = u;
sol.timeevaltime = time2evalIFT;


