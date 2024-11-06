%% script to convert all LUT.txt files 

LUT_cnt = dir(fullfile('E:\roehri\Matlab\test\src', '*FreeSurferColorLUT.txt'));

for k = 1:length(LUT_cnt)
    input_fname = fullfile(LUT_cnt(k).folder, LUT_cnt(k).name);
    output_fname = strrep(input_fname, '.txt', '.json');
    LUT2NiiVue_cmap(input_fname, output_fname)
end
