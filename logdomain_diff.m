function z=logdomain_diff(x,y)
%     x > y
% x MUST be greater than y
%
z = x + log(1 - exp(y-x));
end