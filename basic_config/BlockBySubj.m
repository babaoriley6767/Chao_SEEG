
function BN = BlockBySubj(sbj,project)

switch project
    case 'SEEG_test'
switch sbj
    case '2017_FT_SP_018_xiangxiang'
        BN = {'SZ1','SZ4','test'};     
    case '2019_TT_SP_001_houwei'
        BN = {'SZ1'};
    case 'group1_SJY'
        BN = {'SJY_SZ1','SJY_SZ2'};
    case 'group1_WXL'
        BN = {'WXL_SZ9','WXL_SZ11'};
    case 'group1_WYH'
        BN = {'WYH_SZ2','WYH_SZ3'};
    case 'group1_XX'
        BN = {'XX_SZ1','XX_SZ4'};
    case 'group2_CBB'
        BN = {'CBB_SZ3','CBB_SZ4'};
    case 'group2_YZZ'
        BN = {'YZZ_SZ1','YZZ_SZ2'};
    case 'group2_ZS'
        BN = {'ZS_SZ1','ZS_SZ10'};
    case 'test_bowen'
        BN = {'block1'};
        
end
end


if ~exist('BN')
    BN = {'none'};
elseif isempty(BN)
    BN = {'none'};
end

end