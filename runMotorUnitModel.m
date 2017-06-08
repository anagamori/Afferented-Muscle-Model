%---------------------------
% Author: Akira Nagamori
% Last update: 6/8/17
% REFERENCE:
% Fuglevand et al. 1993
% Fuglevand et al. 2006
%---------------------------

clc
close all
clear all

% load the data from closed-loop simulation of afferented muscle
load output
            
% use simulated ND as an input motor unit model
Input= decimate(output.ND,10); % decimate it to sampling rate of 1000 Hz
Input(Input<0) = 0;

% used-defined parameters for motor unit model
pool_parameter.N = 120; % number of motor units 
pool_parameter.gain = 1.5; % gain 
pool_parameter.ISI_CoV = 0.05; % coefficient of variation of interspike intervals

% run motor unit model
[time_MN,spike_train] = MotorUnitModel(pool_parameter,Input);
