function CI = ConfidenceInterval(A)

N = length(A);
STDmean = std(A)/sqrt(N);
dof = N - 1; %Depends on the problem but this is standard for a CI around a mean.
studentst = tinv([.025 0.975],dof); %tinv is the student's t lookup table for the two-tailed 95% CI ...
CI = studentst*STDmean;

end