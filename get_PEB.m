function PEB = get_PEB(MT_location,BS_location,toa_crbs,aoa_crbs,scenario,transmission_direction)
% ========================================================================================
% Function that calculates the PEB given the MT and BS locations, and the CRBs for the
% AoA angle and ToA. 
% ---------------------------------------------------------------------------------------
%
% INPUTS:
%
% MT_location               :       MT location (matrix with a dimension [2 x 1])
% BS_location               :       BS location (matrix with a dimension [2 x nBS])
% toa_crbs                  :       CRBs for all the ToA (matrix [1 x nBs])
% aoa_crbs                  :       CRBs for all the AoA (matrix [1 x nBs])
% scenario                  :       One from the set {'ToA','AoA','AoA+ToA'}
% transmission_direction    :       Either 'uplink' or 'downlink'
%
% OUTPUTS:
% 
% PEB                       :       Position error bound for the given MT location
%
% ---------------------------------------------------------------------------------------


% Speed of light
c = physconst('LightSpeed');

% Jacobian matrices in a form of function handles (uses MT and BS locations
% as inputs)
toa_jacobian = @(mt_loc,bs_loc) (mt_loc-bs_loc)./sqrt(sum((mt_loc-bs_loc).^2,1))/c;
aoa_jacobian_ul = @(mt_loc,bs_loc) [-mt_loc(2)+bs_loc(2,:);mt_loc(1)-bs_loc(1,:)]./(sum((mt_loc-bs_loc).^2,1));
aoa_jacobian_dl = @(mt_loc,bs_loc) [-bs_loc(2,:)+mt_loc(2);bs_loc(1,:)-mt_loc(1)]./(sum((mt_loc-bs_loc).^2,1));

% Different Jacobian matrix in the case AoA based on the transmission
% direction (only the sign in the AoA-based Jacobian is different)
switch transmission_direction
    case 'uplink'
        aoa_jacobian = aoa_jacobian_ul;
    case 'downlink'
        aoa_jacobian = aoa_jacobian_dl;
end

% Position error bounds for the selected scenario
switch scenario
    case 'ToA'
        J_toa = toa_jacobian(MT_location,BS_location)';
        toa_FIM_in_MT_pos = J_toa'/diag(toa_crbs)*J_toa;
        PEB = sqrt(trace(eye(2)/toa_FIM_in_MT_pos));
    case 'AoA'
        J_aoa = aoa_jacobian(MT_location,BS_location)';
        aoa_FIM_in_MT_pos = J_aoa'/diag(aoa_crbs)*J_aoa;
        PEB = sqrt(trace(eye(2)/aoa_FIM_in_MT_pos));
    case 'AoA+ToA'
        J_aoa_plus_toa = [aoa_jacobian(MT_location,BS_location)';toa_jacobian(MT_location,BS_location)'];
        aoa_plus_toa_FIM = J_aoa_plus_toa'/diag([aoa_crbs,toa_crbs])*J_aoa_plus_toa;
        PEB = sqrt(trace(eye(2)/aoa_plus_toa_FIM));
end



end
