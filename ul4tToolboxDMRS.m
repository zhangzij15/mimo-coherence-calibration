function [pilot,groupCount] = ul4tToolboxDMRS(cfg,slotInFrame,startSymbol)
%UL4TTOOLBOXDMRS Type 1 double-symbol pilots across a contiguous research grid.
%   Larger grids use several valid Toolbox windows. Absolute CRB offsets
%   continue the DM-RS sequence across the windows. Channel, precoding and
%   link mapping still operate over the full allocation as one abstraction.
if nargin<3, startSymbol = 0; end
segments = cfg.ofdm.num_rb;
if isfield(cfg.ul,'dmrs_segments_rb') && ~isempty(cfg.ul.dmrs_segments_rb)
    segments = cfg.ul.dmrs_segments_rb;
end
rank = cfg.ul.rank;
groupCount = 2;
ports = {0,[0 1],[0 1 4],[0 2 4 6]};
carrier = nrCarrierConfig;
carrier.SubcarrierSpacing = cfg.ofdm.scs_Hz/1000;
carrier.NSlot = slotInFrame;
pusch = nrPUSCHConfig;
pusch.NumLayers = rank;
pusch.TransmissionScheme = 'nonCodebook';
pusch.TransformPrecoding = 0;
pusch.MappingType = 'B';
pusch.SymbolAllocation = [startSymbol 14-startSymbol];
pusch.DMRS.DMRSConfigurationType = 1;
pusch.DMRS.DMRSLength = 2;
pusch.DMRS.DMRSAdditionalPosition = 0;
pusch.DMRS.NumCDMGroupsWithoutData = groupCount;
pusch.DMRS.DMRSPortSet = ports{rank};
pusch.DMRS.NIDNSCID = cfg.ul.nid;
pusch.DMRS.NSCID = cfg.ul.nscid;
% Set the pilot EPRE to the same value as the base-MATLAB backend.
beta = 10^(3*(groupCount-1)/20);
pilotSymbols = startSymbol+(1:2);
outsideSymbols = setdiff(1:14,pilotSymbols);
pilot = complex(zeros(cfg.ofdm.num_tones,2,rank));
rbOffset = 0;
for k = 1:numel(segments)
    segmentRBs = segments(k);
    firstCRB = cfg.ul.prb_start_crb+rbOffset;
    carrier.NSizeGrid = segmentRBs;
    carrier.NStartGrid = firstCRB;
    pusch.NSizeBWP = segmentRBs;
    pusch.NStartBWP = firstCRB;
    pusch.PRBSet = 0:segmentRBs-1;
    symbols = nrPUSCHDMRS(carrier,pusch);
    indices = nrPUSCHDMRSIndices(carrier,pusch);
    toneCount = 12*segmentRBs;
    assert(numel(indices)==toneCount*rank && numel(symbols)==numel(indices) && ...
        all(abs(symbols(:))>0),'Unexpected Toolbox Type 1 DMRS symbols or indices.');
    grid = complex(zeros(toneCount,14,rank));
    grid(indices(:)) = beta*symbols(:)./abs(symbols(:));
    assert(~any(grid(:,outsideSymbols,:),'all'), ...
        'Expected DMRS in the first two PUSCH symbols only.');
    pilot(12*rbOffset+(1:toneCount),:,:) = grid(:,pilotSymbols,:);
    rbOffset = rbOffset+segmentRBs;
end
end
