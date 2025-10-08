function extractVarFromStruct(inputStruct)

vars = fieldnames(inputStruct);

for n = 1:length(vars)
    assignin('caller',vars{n},inputStruct.(vars{n}));
end
