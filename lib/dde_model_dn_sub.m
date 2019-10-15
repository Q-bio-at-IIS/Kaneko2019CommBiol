function dn = dde_model_dn_sub(t,n,Z, params)
% 1~4 : DN1~DN4
% 5 : cTEC
% 6 : DP

cTEC_tau = Z(5,1);

dn1 = params.phidn + (params.deltadn1 + params.mudn1 * n(5) ) * n(1);
dn2 = (params.mudn1 * params.r1 * n(5)) * n(1) + (params.deltadn2+ params.mudn2* n(5)) * n(2);
dn3 = (params.mudn2* params.r2 * n(5)) * n(2) +(params.deltadn3 + params.mudn3* n(5))* n(3);
dn4 = (params.mudn3 * params.r3*  n(5)) * n(3) +(params.deltadn4 + params.mudn4* n(5))* n(4);
ctec = params.phictec+(params.deltac + ...
    params.muc1*n(1) + params.muc2*n(2) + params.muc3*n(3) + params.muc4*n(4) ... 
    )* n(5);
dp = params.mudn4*params.r4*n(5)*n(4) + params.r_2 *n(6)*(1-n(6)/params.k_2) - params.mu_2 * cTEC_tau * n(6);


dn = [dn1; dn2; dn3; dn4; ctec; dp];
end