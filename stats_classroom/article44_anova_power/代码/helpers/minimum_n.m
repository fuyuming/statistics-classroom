function answer = minimum_n(target, means, sd, level)
% 连续根用于定位，最终回代整数n，最少每组2个。
effect = sqrt(mean((means-mean(means)).^2))/sd;
if effect == 0 && target > level
    error('总体均数全相同，无法达到高于alpha的功效。');
end
at2 = exact_power(2,means,sd,level);
if at2.power >= target
    continuous_n = NaN;
    n = 2;
else
    upper = 4;
    result = exact_power(upper,means,sd,level);
    while result.power < target
        upper = upper*2;
        result = exact_power(upper,means,sd,level);
    end
    continuous_n = fzero(@(x) power_difference(x,means,sd,level,target),[2 upper]);
    n = ceil(continuous_n);
    result = exact_power(n,means,sd,level);
    while result.power < target
        n = n+1;
        result = exact_power(n,means,sd,level);
    end
    while n > 2
        previous = exact_power(n-1,means,sd,level);
        if previous.power < target
            break
        end
        n = n-1;
    end
end
result = exact_power(n,means,sd,level);
previous_power = NaN;
if n > 2
    previous = exact_power(n-1,means,sd,level);
    previous_power = previous.power;
end
answer = table(target,effect,continuous_n,n,n*numel(means),result.power,previous_power, ...
    'VariableNames',{'target','f','continuous_n','n_per_group','total_n','power','previous_power'});
end

function value = power_difference(n,means,sd,level,target)
result = exact_power(n,means,sd,level);
value = result.power-target;
end
