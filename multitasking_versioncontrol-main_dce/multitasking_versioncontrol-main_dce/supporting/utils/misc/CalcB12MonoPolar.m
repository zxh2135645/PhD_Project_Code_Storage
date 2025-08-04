function [dB12,dMoment1,dMoment2] = CalcB12MonoPolar(dRampUp,dDuration,dRampDown,dGamp1,dGamp2,dDelta)

% % b = 1400 s/mm^2
% dRampUp   = 900;
% dDuration = 9000;
% dRampDown = dRampUp;
% dGamp1 = 74;
% dGamp2 = 74;
% dDelta = 37220;

% Gyromagnetic ratio
m_dGamma     = 42.5756e6*2*pi;
dGamma2      = (m_dGamma * 1.e-15) * (m_dGamma * 1.e-15);

% Gradient slope
dSlope1 = 0.;
dSlope2 = 0.;

% Gradient amplitude mT/m
dGrad1  = 0.;
dGrad2  = 0.;

% Gradient 0th moment
dM1     = 0.;
dM2     = 0.;

% b-value
dB12    = 0.;

% Time steps
dTimeStep  = 10;    % us
dTimeStep2 = dTimeStep  * dTimeStep;
dTimeStep3 = dTimeStep2 * dTimeStep;
dTimeStep4 = dTimeStep2 * dTimeStep2;
dTimeStep5 = dTimeStep2 * dTimeStep3;

NRampUpSteps   = dRampUp/dTimeStep;
NPlateauSteps  = (dDuration-dRampUp)/dTimeStep;
NRampDownSteps = dRampDown/dTimeStep;

% Start calculation
% RampUp
Slope1 = dGamp1/dRampUp;
Slope2 = dGamp2/dRampUp;
dSlope1 = Slope1;
dSlope2 = Slope2;
for step = 1:NRampUpSteps
    dB12 = dB12       +            dM1     * dM2                                             * dTimeStep...
                      + 1./2.   * (dM1     * dGrad2  + dM2 * dGrad1)                         * dTimeStep2...
                      + 1./6.   * (dM1     * dSlope2 + dM2 * dSlope1 + 2. * dGrad1 * dGrad2) * dTimeStep3...
                      + 1./8.   * (dGrad1  * dSlope2 + dGrad2 * dSlope1)                     * dTimeStep4...
                      + 1./20.  * (dSlope1 * dSlope2)                                        * dTimeStep5;
    
    dM1 = dM1 + dGrad1 * dTimeStep + 1./2. * dSlope1 * dTimeStep2;
    dM2 = dM2 + dGrad2 * dTimeStep + 1./2. * dSlope2 * dTimeStep2;
    dGrad1 = dGrad1 + dSlope1 * dTimeStep;
    dGrad2 = dGrad2 + dSlope2 * dTimeStep;
end

% Gradient plateau
dSlope1 = 0;
dSlope2 = 0;
for step = 1:NPlateauSteps
    dB12 = dB12       +            dM1     * dM2                                             * dTimeStep ...
                      + 1./2.   * (dM1     * dGrad2  + dM2 * dGrad1)                         * dTimeStep2...
                      + 1./6.   * (dM1     * dSlope2 + dM2 * dSlope1 + 2. * dGrad1 * dGrad2) * dTimeStep3...
                      + 1./8.   * (dGrad1  * dSlope2 + dGrad2 * dSlope1)                     * dTimeStep4...
                      + 1./20.  * (dSlope1 * dSlope2)                                        * dTimeStep5;
    
    dM1 = dM1 + dGrad1 * dTimeStep + 1./2. * dSlope1 * dTimeStep2;
    dM2 = dM2 + dGrad2 * dTimeStep + 1./2. * dSlope2 * dTimeStep2;
    dGrad1 = dGrad1 + dSlope1 * dTimeStep;
    dGrad2 = dGrad2 + dSlope2 * dTimeStep;
end

% RampDown
Slope1 = -dGamp1/dRampDown;
Slope2 = -dGamp2/dRampDown;
dSlope1 = Slope1;
dSlope2 = Slope2;
for step = 1:NRampDownSteps
    dB12 = dB12       +            dM1     * dM2                                             * dTimeStep ...
                      + 1./2.   * (dM1     * dGrad2  + dM2 * dGrad1)                         * dTimeStep2...
                      + 1./6.   * (dM1     * dSlope2 + dM2 * dSlope1 + 2. * dGrad1 * dGrad2) * dTimeStep3...
                      + 1./8.   * (dGrad1  * dSlope2 + dGrad2 * dSlope1)                     * dTimeStep4...
                      + 1./20.  * (dSlope1 * dSlope2)                                        * dTimeStep5;
    
    dM1 = dM1 + dGrad1 * dTimeStep + 1./2. * dSlope1 * dTimeStep2;
    dM2 = dM2 + dGrad2 * dTimeStep + 1./2. * dSlope2 * dTimeStep2;
    dGrad1 = dGrad1 + dSlope1 * dTimeStep;
    dGrad2 = dGrad2 + dSlope2 * dTimeStep;
end

% Time between gradient pulses
dB12 = dB12 + dM1 * dM2 * dDelta;

temp = dM1 * dM2;

% RampUp
Slope1 = -dGamp1/dRampUp;
Slope2 = -dGamp2/dRampUp;
dSlope1 = Slope1;
dSlope2 = Slope2;
for step = 1:NRampUpSteps
    dB12 = dB12       +            dM1     * dM2                                             * dTimeStep...
                      + 1./2.   * (dM1     * dGrad2  + dM2 * dGrad1)                         * dTimeStep2...
                      + 1./6.   * (dM1     * dSlope2 + dM2 * dSlope1 + 2. * dGrad1 * dGrad2) * dTimeStep3...
                      + 1./8.   * (dGrad1  * dSlope2 + dGrad2 * dSlope1)                     * dTimeStep4...
                      + 1./20.  * (dSlope1 * dSlope2)                                        * dTimeStep5;
    
    dM1 = dM1 + dGrad1 * dTimeStep + 1./2. * dSlope1 * dTimeStep2;
    dM2 = dM2 + dGrad2 * dTimeStep + 1./2. * dSlope2 * dTimeStep2;
    dGrad1 = dGrad1 + dSlope1 * dTimeStep;
    dGrad2 = dGrad2 + dSlope2 * dTimeStep;
end

% Gradient plateau
dSlope1 = 0;
dSlope2 = 0;
for step = 1:NPlateauSteps
    dB12 = dB12       +            dM1     * dM2                                             * dTimeStep ...
                      + 1./2.   * (dM1     * dGrad2  + dM2 * dGrad1)                         * dTimeStep2...
                      + 1./6.   * (dM1     * dSlope2 + dM2 * dSlope1 + 2. * dGrad1 * dGrad2) * dTimeStep3...
                      + 1./8.   * (dGrad1  * dSlope2 + dGrad2 * dSlope1)                     * dTimeStep4...
                      + 1./20.  * (dSlope1 * dSlope2)                                        * dTimeStep5;
    
    dM1 = dM1 + dGrad1 * dTimeStep + 1./2. * dSlope1 * dTimeStep2;
    dM2 = dM2 + dGrad2 * dTimeStep + 1./2. * dSlope2 * dTimeStep2;
    dGrad1 = dGrad1 + dSlope1 * dTimeStep;
    dGrad2 = dGrad2 + dSlope2 * dTimeStep;
end

% RampDown
Slope1 = dGamp1/dRampDown;
Slope2 = dGamp2/dRampDown;
dSlope1 = Slope1;
dSlope2 = Slope2;
for step = 1:NRampDownSteps
    dB12 = dB12       +            dM1     * dM2                                             * dTimeStep ...
                      + 1./2.   * (dM1     * dGrad2  + dM2 * dGrad1)                         * dTimeStep2...
                      + 1./6.   * (dM1     * dSlope2 + dM2 * dSlope1 + 2. * dGrad1 * dGrad2) * dTimeStep3...
                      + 1./8.   * (dGrad1  * dSlope2 + dGrad2 * dSlope1)                     * dTimeStep4...
                      + 1./20.  * (dSlope1 * dSlope2)                                        * dTimeStep5;
    
    dM1 = dM1 + dGrad1 * dTimeStep + 1./2. * dSlope1 * dTimeStep2;
    dM2 = dM2 + dGrad2 * dTimeStep + 1./2. * dSlope2 * dTimeStep2;
    dGrad1 = dGrad1 + dSlope1 * dTimeStep;
    dGrad2 = dGrad2 + dSlope2 * dTimeStep;
%     dSlope1 = dSlope1 + dSign * Slope1;
%     dSlope2 = dSlope2 + dSign * Slope2;
end

% Convert gradient moments to mT/m ms
dMoment1 = dM1 / 1.e3;
dMoment2 = dM2 / 1.e3;

% Convert b-value to s/mm^2
dB12 = dB12 * dGamma2; 

