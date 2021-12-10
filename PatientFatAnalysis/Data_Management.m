clear all;
close all;

count = 0;
fname = 'D:\PatientFatAnalysis\Yinyin_Patient_data';
subjects_f = glob(cat(2, fname, '/*'));
qual_f = {};

for i = 1:length(subjects_f)
mod_f = glob(cat(2, subjects_f{i}, '*'));
if length(mod_f) > 1
count = count + 1;qual_f{count} = subjects_f{i};
end
end