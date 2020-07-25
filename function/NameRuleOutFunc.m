  function RuleOutLabel = NameRuleOutFunc(Names, RuleOutCell)
RuleOutLabel = zeros(length(Names), 1);

for i = 1:length(Names)
    switch Names{i}
        case RuleOutCell
            RuleOutLabel(i) = 1;
        otherwise
            RuleOutLabel(i) = 0;
    end
end
end