function [DOCID,GID] = getGoogleSheetInfo(table,sheet)

if isempty(sheet)
    sheet = 'generic';
else
end

switch table
    case 'math_network'
        DOCID = '1rLFTsiDXyDvckau5epvuUfvfD9LGibWn6p8QJBms2NY';
        
        switch sheet
            case 'MMR'
                GID = '0';
                
            case 'UCLA'
                GID = '0';
                
            case 'MFA'
                GID = '371325198';
                
            case 'Calculia'
                GID = '1950213300';
                
            case 'Memoria'
                GID = '1957515618';
                
            case 'Number_comparison'
                GID = '2047525626';
                
            case 'Calculia_production'
                GID = '629658513';
                
            case 'Calculia_China'
                GID = '1555855901';
                
            case 'Scramble'
                GID = '1045421744';
                
            case 'Logo'
                GID = '1045421744';
                
            case '7_heaven'
                GID = '1838563006';
                
            case 'VTC'
                GID = '2116837816';
                
            case 'AllCateg'
                GID = '2116837816';
                
            case 'Flanker'
                GID = '481390136';
                
            case 'Posner'
                GID = '108675259';
                
            case 'Numblet'
                GID = '371325198';
                
            case 'math_network_cohort'
                GID = '1325879124';
                
            case 'cohort'
                GID = '282063015';
                
        end
        
    case 'subjects_by_task'
        DOCID = '1rH0cpsrMXo9dh68eIVme0wjuKCMca5pHRwjowLJVk_8';
        GID = '0';
        
    case 'block_by_sbj'
        DOCID = '1mOI89bk9FR47h0neVU2voF2EBWrUGFon54AQQiBEbUw';
        GID = '0';
        
    case 'bad_channel_by_block'
        DOCID = '1mOI89bk9FR47h0neVU2voF2EBWrUGFon54AQQiBEbUw';
        GID = '1534541388';
        
    case 'demographics'
        DOCID = '1H_oM_GuIMX5RTDeRZrOSBb3O1zms73xmUd_05de0cqA';
        GID = '683763705';
        
    case 'chan_names_ppt'
        DOCID = '1RSgUAGnngW-GLqtSkse05lZ3JGzJ6iImAQeeUrZzqwE';
        switch sheet
            case 'chan_names_fs_figures'
                GID = '1158910949';
            case 'chan_names_ppt_log'
                GID = '0';
                
        end
        
        
        
        
    case 'simulatenous_coverage_math'
        DOCID = '1UKqc0m6guvx5-tuXtfapWkZun4Hx7lJBbl0y7gfyGOs';
         switch sheet
            case 'math_selective_only'
                GID = '486842448';
         end
         
         
         
        
end

end

