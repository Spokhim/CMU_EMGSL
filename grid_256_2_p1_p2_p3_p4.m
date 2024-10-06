function iRemap = grid_256_2_p1_p2_p3_p4(subj)
%GRID_256_2_P1_P2_P3_P4 For re-indexing `uni` in _Synchronized.mat files sessions, by subject.
arguments
    subj (1,1) string {mustBeTextScalar, mustBeMember(subj, ["MCP01", "MCP02", "MCP03", "MCP04", "MCP05", "MCP06", "MCP07"])}
end
switch subj
    case {"MCP01", "MCP02", "MCP03", "MCP04"}
        iRemap = [...
            1:32, ... % P1 EXT
            161:192, ... % P1 FLX
            33:64, ... % P2 EXT
            129:160, ... % P2 FLX
            65:96, ... % P3 EXT
            225:256, ... % P3 FLX
            97:128, ... % P4 EXT
            193:224 ... % P4 FLX
            ];
    case {"MCP05", "MCP06"}
        iRemap = [1:128, 193:256, 129:192]; % Accidental switching  of SAGAs in P3/P4 case.
    case "MCP07"
        iRemap = 1:128; % No remap required
end


end