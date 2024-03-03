##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2023, Sören Jakob Gerke

## calculate normalized concentration from two point calibration
## c = c0 + (csat - c0) * cn
##
## Author: Sören J. Gerke
##

function [cn ext] = calc_cn (phi_dat, cref, method, sig, testplots)

  ext = [];
  conc = cn = [];
  c_sat_ref = cref(2);
  c_des_ref = cref(1);
  cscale = c_sat_ref - c_des_ref;

  ## quenching measurement Ic: desorbed liquid in contact with air
  phi = (phi_dat{1}); #
  phi(isnan (phi)) = 0;
##  phi(phi <= 0) = 1e-6;

  ## minimal quenching reference Ic0: desorbed liquid in equilibrium with nitrogen
  phi_des = (phi_dat{2}); #
  phi_des(isnan (phi_des)) = 0;
##  phi_des(phi_des <= 0) = 1e-6;

  ## maximum quenching reference Ic1: saturated liquid in equilibrium with air
  phi_sat = (phi_dat{3}); #
  phi_sat(isnan (phi_sat)) = 0;
##  phi_sat(phi_sat <=0 ) = 1e-6;

  if (! isempty (sig))
    ## smoothing of calibration reference points
    phi_des = imsmooth (phi_des, sig);
    phi_sat = imsmooth (phi_sat, sig);
  endif

  if (testplots)
    plot_map (phi);
    title ("phi");
    plot_map (phi_des);
    title ("phi_des");
    plot_map (phi_sat);
    title ("phi_sat");
  endif

  ##  linear and SV2 give the same result, SV1 only if phi_des == phi_zero and phi_sat really is at saturation
  switch (method)

    case "linear" ## linear quenching

      K1 = (c_sat_ref - c_des_ref) .* (phi_des .* phi_sat ./ (phi_des - phi_sat));
      K0 = (c_sat_ref - c_des_ref) .* (         - phi_sat ./ (phi_des - phi_sat)) + c_des_ref;

      if ((c_sat_ref==1) && (c_des_ref==0))
        K1 = phi_des .* phi_sat ./ (phi_des - phi_sat);
        K0 = - K1 ./ phi_des;
      endif

      cn = K1 .* (1 ./ phi) + K0;

      ext.K1 = K1;
      ext.K0 = K0;

    case "SV1" ## weak excitation Stern-Volmer formulation

      KSV = 1 ./ cscale .* (phi_des ./ phi_sat - 1); ## pixel wise scale factor from ref data
      conc = 1 ./ KSV .* (phi_des ./ phi - 1) + cref(1);

    case "SV_const" ## weak excitation Stern-Volmer formulation with given KSV

      KSV = cref(3);
      phi_zero = 1.02 * phi_des; # TODO: find zero quenching level from KSV and valid reference
      phi_zero = phi_sat * (KSV * c_sat_ref + 1);
      conc = 1. / KSV .* (phi_zero ./ phi - 1);
      ext.phi_zero = phi_zero;

    case "SV2" ## weak excitation Stern-Volmer from linear relation

      ## may be used for KSV analysis when physical concentrations are known !remove black response!
      ## 1 / phi_cal = ( m * c_cal + 1 ) / phi_zero ## inverse photon flux is sum of quenching part and base inersve photon flux at maximum quenching
      m = 1 ./ cscale .* (1./phi_sat - 1./phi_des); # KSV / phi_zero
      ##y = m * x + n # 1/phi_cal
      ##n = y1 - m * x1 # 1/phi_zero
      ##n = 1/phi_cal0 - m * cref(1)
      ##n = 1/phi_cal1 - m * cref(2)
  ##    phi_zero = 1 ./ (1./phi_des - m*cref(1));
      phi_zero = 1 ./ (1./phi_sat - m*cref(2));
      KSV = m .* phi_zero;

      conc = (phi_zero ./ phi - 1) ./ KSV;
      ext.KSV = KSV;
      ext.phi_zero = phi_zero;

      if testplots
        ## check for one pixel
        idx_x = 2429; idx_y = 423;
        figure ()
        plot ([cref(1) cref(2)], [1/phi_des(idx_y,idx_x) 1/phi_sat(idx_y,idx_x)], "bx-"); hold on;
        plot ([cref(1) cref(2)], 1 ./ [phi(idx_y,idx_x) phi(idx_y,idx_x)] ,"r-");
        xlabel ("c in mg/l")
        ylabel ("1 / phi")
        hold off
        figure ()
        plot ([cref(1) cref(2)], [phi_des(idx_y,idx_x) phi_sat(idx_y,idx_x)], "bx-"); hold on;
        plot ([cref(1) cref(2)], [phi(idx_y,idx_x) phi(idx_y,idx_x)] ,"r-");
        xlabel ("c in mg/l")
        ylabel ("phi")
        ##
        plot_map (m)
        caxis ([0 10])
        colorbar
        plot_map (phi_zero)
        caxis ([0 0.05])
        colorbar
        plot_map (KSV)
        caxis ([0 0.4])
        colorbar
      endif

    case "nonlinear"

      ##  estimate phi0 from linear c-phi-cal
      m1 = 1 ./ cscale .* (phi_sat - phi_des); # KSV / phi_zero
      phi0 = (phi_sat - m1 .* cref(2));
      ##  combined quenching
      if (cref(1) == 0)
        psi1 = zeros (size (phi0));
      else
        psi1 = (phi0 ./ phi_des-1) ./ cref(1);
      endif
      psi2 = (phi0 ./ phi_sat - 1) ./ cref(2);
      m2 = (psi2 - psi1) / cscale; ## m = KS * KD
      psi0 = (psi2 - m2 * cref(2)); ## n = KS + KS

      conc = 1 ./ (2 .* m2) .* (- psi0 + sqrt (psi0 .^ 2 - 4 .* m2 .* (- (1 * phi0 ./ phi - 1)))); # pos. solution
      conc = real (conc);

      if testplots
        ## check for one pixel
        idx_x = 2429; idx_y = 423;
        mp = 1 ./ cscale .* (phi_sat(idx_y,idx_x) - phi_des(idx_y,idx_x)); # KSV / phi_zero
        F0 = phi_sat(idx_y,idx_x) - mp .* cref(2); # estimate from linear
        F = phi(idx_y,idx_x);
        pcal = [[cref(1) (F0./phi_des(idx_y,idx_x)-1)./cref(1)]; [cref(2) (F0./phi_sat(idx_y,idx_x)-1)./cref(2)]];
        pline = createLine (pcal(1,:), pcal(2,:));
        figure ();
        plot (pcal(:,1), pcal(:,2),"bx");
        hold on;
        drawLine (pline, "k");
        xlabel ("c in mg/l");
        ylabel ("(F0 / F -1) / Q");
        mpp = (pcal(2,2) - pcal(1,2)) / cscale;
        n = pcal(2,2) - mpp .* cref(2);
        coeff = [mpp, n, - (F0 / F - 1)];
        c = roots (coeff);
        c = c(c > 0);
        plot (c, (F0 / F - 1) / c, "ro")
        ##
        plot_map (real (conc));
        caxis ([0 2]);
        colorbar;
        plot_map (psi2);
        caxis ([0 0.3]);
        colorbar;
      endif
    otherwise
      error (["calc_cn: method unknown"])

  endswitch

  ##
  if (! isempty (conc))
##    cn = (conc - cref(1)) ./ cscale;
    cn = norm_conc (conc, c_sat_ref, c_des_ref);
  endif

  ## exclude unphysical
  cn(cn < -0.025) = -0.025; # but allow bulk concentration to fluctuate around zero
  cn(cn > 1) = 1;

endfunction
