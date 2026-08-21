function result = runValidatedFullStackISMWavefront(varargin)
%RUNVALIDATEDFULLSTACKISMWAVEFRONT Run the slower complete validation suite.
%
%   result = runValidatedFullStackISMWavefront()
%   result = runValidatedFullStackISMWavefront('flatField', flatFieldPtu)

    alignmentCsv = ['D:\Luminosa\Data\PSF_batch_outputs\ISM_Aberation2_73\' ...
        'xz_yz_plots\x0_4uW_0_19collar_80mmlens_20260515_155744_frame_alignment.csv'];
    args = varargin;
    if ~hasNameValue(args,'runProfiles')
        args=[args {'runProfiles',true}];
    end
    if ~hasNameValue(args,'runSignProfiles')
        args=[args {'runSignProfiles',true}];
    end
    if ~hasNameValue(args,'nBootstrap')
        args=[args {'nBootstrap',20}];
    end
    if ~hasNameValue(args,'maxIter')
        args=[args {'maxIter',8}];
    end
    if ~hasNameValue(args,'requireValidationPrerequisites')
        args=[args {'requireValidationPrerequisites',true}];
    end
    if ~hasNameValue(args,'airBeadAxialSamples')
        args=[args {'airBeadAxialSamples',5}];
    end
    if ~hasNameValue(args,'outputDir')
        root=fileparts(fileparts(mfilename('fullpath')));
        args=[args {'outputDir',fullfile(root,'output_matlab', ...
            'full_stack_ism_wavefront_validated')}];
    end
    result=estimateFullStackISMWavefront(alignmentCsv,args{:});
end

function tf = hasNameValue(args,name)
    tf=false;
    for k=1:2:numel(args)
        if strcmpi(args{k},name)
            tf=true;
            return;
        end
    end
end
