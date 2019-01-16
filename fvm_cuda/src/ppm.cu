#ifdef PPM
__device__ __inline__ void slope_limiter_ppm(const real qL, const real qC, const real qR, real *wL, real *wR) {
	/* The PPM slope limiter */
	real slopeR = *wR - qC;
	real slopeL = qC - *wL;

	if (slopeR*slopeL <= 0 ) {
		*wL = qC;
		*wR = qC;
	}
	else if ( 3*(slopeL - slopeR) > fabs(slopeL + slopeR)) {
		*wL = qC - 2*slopeR;
	}
	else if ( 3*(slopeL - slopeR) < -fabs(slopeL + slopeR)) {
		*wR = qC + 2*slopeL;
	}

	*wL = MIN2( MAX2( qC, qL ), MAX2( *wL, MIN2( qC, qL )) );
	*wR = MIN2( MAX2( qC, qR ), MAX2( *wR, MIN2( qC, qR )) );

	return;

}
#endif
#ifdef PPM
__global__ void reconstruct(real *cons, real *UL, real *UR, real *dx,
        int dir1,int nx1, int nx2, int nx3, int size_x1, int size_x12,
        int nf,int ntot, int offset, real g1, real dt) {
	/*
	 * Piecewise parabolic reconstruction
	 *
	 * primitive slopes: d_w = (d_rho, d_vx, d_vy, d_vz, d_p)
	 * characteristic slopes:
	 * lam = vx - cs  -->  d_xi = .5*( rho/cs * d_vx - d_p / cs^2)
	 * lam = vx       -->  d_xi = d_rho - d_p/cs^2
	 * lam = vx       -->  d_xi = d_vy
	 * lam = vx       -->  d_xi = d_vz
	 * lam = vx + cx  -->  d_xi = .5*(rho/cs * d_vx + d_p/cs^2)
	 */
    int i,j,k,n,indx,indxm,indxp,indxm2,indxp2;
    int il, iu, jl, ju, kl, ku;
    real dL, uL, uL2,uL3,pL,eL,sL;
    real dLL, uLL, uLL2,uLL3,pLL,eLL,sLL;
    real dC,uC,uC2,uC3,pC,eC,sC;
    real dR,uR,uR2,uR3,pR,eR,sR;
    real dRR,uRR,uRR2,uRR3,pRR,eRR,sRR;
    real cs,cs2;
    real prims_0_L, prims_1_L, prims_2_L, prims_3_L, prims_4_L;
    real prims_0_R, prims_1_R, prims_2_R, prims_3_R, prims_4_R;
    real prims_0_C, prims_1_C, prims_2_C, prims_3_C, prims_4_C;
    real char_0_L, char_1_L, char_2_L, char_3_L, char_4_L;
    real char_0_R, char_1_R, char_2_R, char_3_R, char_4_R;
    real char_0_C, char_1_C, char_2_C, char_3_C, char_4_C;

    real w_0_L, w_1_L, w_2_L, w_3_L, w_4_L;
    real w_0_R, w_1_R, w_2_R, w_3_R, w_4_R;
    real d_6, u_6, u2_6, u3_6, p_6;
    real dw_0_L, dw_1_L, dw_2_L, dw_3_L, dw_4_L;
    real dw_0_R, dw_1_R, dw_2_R, dw_3_R, dw_4_R;
    real dw_0_C, dw_1_C, dw_2_C, dw_3_C, dw_4_C;
    real dxm, dxc, dxp,dxm2,dxp2;
    real slopeL,slopeR,ke,r;
    real dtdx;


    int dir2, dir3;
    /* 1->2->3
     * 2->3->1
     * 3->1->2
     */
    dir2 = (dir1)%3 + 1;
    dir3 = (dir2)%3 + 1;

    if (dir1 == 1) {
        il = -2; iu = nx1+2;
        jl = -NGHX2; ju = nx2+NGHX2;
        kl = -NGHX3; ku = nx3 + NGHX3;
    }
    else if (dir1 == 2) {
        il = -NGHX1; iu = nx1+NGHX1;
        jl = -2; ju = nx2+2;
        kl = -NGHX3; ku = nx3 + NGHX3;
    }
    else {
    	il = -NGHX1; iu = nx1+NGHX1;
		jl = -NGHX2; ju = nx2+NGHX2;
		kl = -2; ku = nx3 + 2;
    }

    for(indx = blockIdx.x*blockDim.x + threadIdx.x; indx<ntot; indx+=blockDim.x*gridDim.x) {
    	unpack_indices(indx,&i,&j,&k,size_x1,size_x12);
        if ((i>=il)&&(i<iu)&&(j>=jl)&&(j<ju)&&(k>=kl)&&(k<ku)) {
            if (dir1 == 1) {
                indxm = indx - 1 ; // (i-1,j,k)
                indxp = indx + 1 ; // (i+1,j,k)
                indxm2 = indx - 2 ; // (i-2,j,k)
                indxp2 = indx + 2 ; // (i+2,j,k)
                dxm = dx[i-1];
                dxc = dx[i];
                dxp = dx[i+1];
                dxm2 = dx[i-2];
                dxp2 = dx[i+2];
            }
            else if (dir1 == 2) {
                indxm = indx - size_x1;
                indxp = indx + size_x1;
                indxm2 = indx - 2*size_x1;
                indxp2 = indx + 2*size_x1;
                dxm = dx[j-1];
                dxc = dx[j];
                dxp = dx[j+1];
                dxm2 = dx[j-2];
                dxp2 = dx[j+2];
            }
            else {
            	indxm = indx - size_x12;
            	indxp = indx + size_x12;
            	indxm2 = indx - 2*size_x12;
            	indxp2 = indx + 2*size_x12;
            	dxm = dx[k-1];
            	dxc = dx[k];
            	dxp = dx[k+1];
            	dxm2 = dx[k-2];
            	dxp2 = dx[k+2];
            }
            dtdx = dt/dxc;

            /* Left neighbor */
            dL  = cons[indxm + 0*ntot];
            uL  = cons[indxm + dir1*ntot]/dL;
            uL2 = cons[indxm + dir2*ntot]/dL;
            uL3 = cons[indxm + dir3*ntot]/dL;
            eL  = cons[indxm + 4*ntot];
            pL = (eL-.5*dL*(uL*uL+uL2*uL2+uL3*uL3))*g1;

            if (dL < PRESSUREFLOOR) dL = PRESSUREFLOOR;
            if (pL < PRESSUREFLOOR) pL = PRESSUREFLOOR;


            /* Left neighbor of left neighbor*/
            dLL  = cons[indxm2 + 0*ntot];
            uLL  = cons[indxm2 + dir1*ntot]/dLL;
            uLL2 = cons[indxm2 + dir2*ntot]/dLL;
            uLL3 = cons[indxm2 + dir3*ntot]/dLL;
            eLL  = cons[indxm2 + 4*ntot];
            pLL = (eLL-.5*dLL*(uLL*uLL+uLL2*uLL2+uLL3*uLL3))*g1;

            if (dLL < PRESSUREFLOOR) dLL = PRESSUREFLOOR;
            if (pLL < PRESSUREFLOOR) pLL = PRESSUREFLOOR;

            /* Center cell */
            dC  = cons[indx + 0*ntot];
            uC  = cons[indx + dir1*ntot]/dC;
            uC2 = cons[indx + dir2*ntot]/dC;
            uC3 = cons[indx + dir3*ntot]/dC;
            eC  = cons[indx + 4*ntot];
            pC = (eC-.5*dC*(uC*uC+uC2*uC2+uC3*uC3))*g1;

            if (dC < PRESSUREFLOOR) dC = PRESSUREFLOOR;
            if (pC < PRESSUREFLOOR) pC = PRESSUREFLOOR;


            /* Right neighbor */
            dR  = cons[indxp + 0*ntot];
            uR  = cons[indxp + dir1*ntot]/dR;
            uR2 = cons[indxp + dir2*ntot]/dR;
            uR3 = cons[indxp + dir3*ntot]/dR;
            eR  = cons[indxp + 4*ntot];
            pR = (eR-.5*dR*(uR*uR+uR2*uR2+uR3*uR3))*g1;
            if (dR < PRESSUREFLOOR) dR = PRESSUREFLOOR;
            if (pR < PRESSUREFLOOR) pR = PRESSUREFLOOR;

            /* Right neighbor of right neighbor*/
            dRR  = cons[indxp2 + 0*ntot];
            uRR  = cons[indxp2 + dir1*ntot]/dRR;
            uRR2 = cons[indxp2 + dir2*ntot]/dRR;
            uRR3 = cons[indxp2 + dir3*ntot]/dRR;
            eRR  = cons[indxp2 + 4*ntot];
            pRR = (eRR-.5*dRR*(uRR*uRR+uRR2*uRR2+uRR3*uRR3))*g1;
            if (dRR < PRESSUREFLOOR) dRR = PRESSUREFLOOR;
            if (pRR < PRESSUREFLOOR) pRR = PRESSUREFLOOR;



            /* Left cell */
            /* C->L, L-> LL, R-> C */
            cs = sqrt((g1+1)*pL/dL);
            cs2 = cs*cs;
            /* Density */
            slopeR = (dC-dL) * 2./(dxm + dxc);
            slopeL = (dL-dLL) * 2./(dxm + dxm2);
            r = slope_limiter(slopeL,slopeR);
            prims_0_C = dxm * r * slopeR;

            /* ux1 */
            slopeR = (uC-uL) * 2./(dxm + dxc);
            slopeL = (uL-uLL) * 2./(dxm + dxm2);
            r = slope_limiter(slopeL,slopeR);
            prims_1_C = dxm * r * slopeR;

            /* ux2 */
            slopeR = (uC2-uL2) * 2./(dxm + dxc);
            slopeL = (uL2-uLL2) * 2./(dxm + dxm2);
			r = slope_limiter(slopeL,slopeR);
			dw_2_L = dxm * r * slopeR;

           /* ux3 */
            slopeR = (uC2-uL2) * 2./(dxm + dxc);
            slopeL = (uL2-uLL2) * 2./(dxm + dxm2);
			r = slope_limiter(slopeL,slopeR);
			dw_3_L = dxm * r * slopeR;

			/* Pressure */
            slopeR = (pC-pL) * 2./(dxm + dxc);
            slopeL = (pL-pLL) * 2./(dxm + dxm2);
			r = slope_limiter(slopeL,slopeR);
			prims_4_C = dxm * r * slopeR;

			/* Characteristics */

            char_0_L = -.5*(dL*prims_1_L/cs - prims_4_L/cs2);
            char_1_L = prims_0_L - prims_4_L/cs2;
            char_4_L = .5*(dL*prims_1_L/cs + prims_4_L/cs2);

            char_0_R = -.5*(dL*prims_1_R/cs - prims_4_R/cs2);
            char_1_R = prims_0_R - prims_4_R/cs2;
            char_4_R = .5*(dL*prims_1_R/cs + prims_4_R/cs2);

            char_0_C = -.5*(dL*prims_1_C/cs - prims_4_C/cs2);
            char_1_C = prims_0_C - prims_4_C/cs2;
            char_4_C = .5*(dL*prims_1_C/cs + prims_4_C/cs2);

            /* Limit */
            // Check dxc and rest
            r = slope_limiter(char_0_L,char_0_R);
            char_0_R = dxm * r * char_0_R;
            r = slope_limiter(char_0_R, char_0_C);
            char_0_C = dxm * r * char_0_C;


            /* Project back */

            dw_0_L = char_0_C + char_1_C + char_4_C;
            dw_1_L = cs/dL*(char_4_C - char_0_C);
            dw_4_L = cs2*(char_0_C + char_4_C);



            /* Right cell */
            /* C->R, L->C, R->RR */

            cs = sqrt((g1+1)*pR/dR);
            cs2 = cs*cs;


            /* Density */
            slopeR = (dRR-dR) * 2./(dxp + dxp2);
            slopeL = (dR-dC) * 2./(dxp + dxc);
            r = slope_limiter(slopeL,slopeR);
            prims_0_C = dxp * r * slopeR;

            /* ux1 */
            slopeR = (uRR-uR) * 2./(dxp + dxp2);
            slopeL = (uR-uC) * 2./(dxp + dxc);
            r = slope_limiter(slopeL,slopeR);
            prims_1_C = dxp * r * slopeR;

            /* ux2 */
            slopeR = (uRR2-uR2) * 2./(dxp + dxp2);
            slopeL = (uR2-uC2) * 2./(dxp + dxc);
			r = slope_limiter(slopeL,slopeR);
			dw_2_R = dxp * r * slopeR;

           /* ux3 */
            slopeR = (uRR3-uR3) * 2./(dxp + dxp2);
            slopeL = (uR3-uC3) * 2./(dxp + dxc);
			r = slope_limiter(slopeL,slopeR);
			dw_3_R = dxp * r * slopeR;

			/* Pressure */
            slopeR = (pRR-pR) * 2./(dxp + dxp2);
            slopeL = (pR-pC) * 2./(dxp + dxc);
			r = slope_limiter(slopeL,slopeR);
			prims_4_C = dxp * r * slopeR;

			/* Characteristics */

            char_0_L = -.5*(dR*prims_1_L/cs - prims_4_L/cs2);
            char_1_L = prims_0_L - prims_4_L/cs2;
            char_4_L = .5*(dR*prims_1_L/cs + prims_4_L/cs2);

            char_0_R = -.5*(dR*prims_1_R/cs - prims_4_R/cs2);
            char_1_R = prims_0_R - prims_4_R/cs2;
            char_4_R = .5*(dR*prims_1_R/cs + prims_4_R/cs2);

            char_0_C = -.5*(dR*prims_1_C/cs - prims_4_C/cs2);
            char_1_C = prims_0_C - prims_4_C/cs2;
            char_4_C = .5*(dR*prims_1_C/cs + prims_4_C/cs2);

            /* Limit */
            // Check dxc and rest
            r = slope_limiter(char_0_L,char_0_R);
            char_0_R = dxp * r * char_0_R;
            r = slope_limiter(char_0_R, char_0_C);
            char_0_C = dxp * r * char_0_C;


            /* Project back */

            dw_0_R = char_0_C + char_1_C + char_4_C;
            dw_1_R = cs/dR*(char_4_C - char_0_C);
            dw_4_R = cs2*(char_0_C + char_4_C);


            /* Center cell */

            cs = sqrt((g1+1)*pC/dC);
            cs2 = cs*cs;


            /* Density */
            slopeR = (dR-dC) * 2./(dxc + dxp);
            slopeL = (dC-dL) * 2./(dxc + dxm);
            r = slope_limiter(slopeL,slopeR);
            prims_0_C = dxc * r * slopeR;

            /* ux1 */
            slopeR = (uR-uC) * 2./(dxc + dxp);
            slopeL = (uC-uL) * 2./(dxc + dxm);
            r = slope_limiter(slopeL,slopeR);
            prims_1_C = dxc * r * slopeR;

            /* ux2 */
			slopeR = (uR2-uC2) * 2./(dxc + dxp);
			slopeL = (uC2-uL2) * 2./(dxc + dxm);
			r = slope_limiter(slopeL,slopeR);
			dw_2_C = dxc * r * slopeR;

           /* ux3 */
			slopeR = (uR3-uC3) * 2./(dxc + dxp);
			slopeL = (uC3-uL3) * 2./(dxc + dxm);
			r = slope_limiter(slopeL,slopeR);
			dw_3_C = dxc * r * slopeR;

			/* Pressure */
			slopeR = (pR-pC) * 2./(dxc + dxp);
			slopeL = (pC-pL) * 2./(dxc + dxm);
			r = slope_limiter(slopeL,slopeR);
			prims_4_C = dxc * r * slopeR;

			/* Characteristics */

            char_0_L = -.5*(dC*prims_1_L/cs - prims_4_L/cs2);
            char_1_L = prims_0_L - prims_4_L/cs2;
            char_4_L = .5*(dC*prims_1_L/cs + prims_4_L/cs2);

            char_0_R = -.5*(dC*prims_1_R/cs - prims_4_R/cs2);
            char_1_R = prims_0_R - prims_4_R/cs2;
            char_4_R = .5*(dC*prims_1_R/cs + prims_4_R/cs2);

            char_0_C = -.5*(dC*prims_1_C/cs - prims_4_C/cs2);
            char_1_C = prims_0_C - prims_4_C/cs2;
            char_4_C = .5*(dC*prims_1_C/cs + prims_4_C/cs2);

            /* Limit */
            // Check dxc and rest
            r = slope_limiter(char_0_L,char_0_R);
            char_0_R = dxc * r * char_0_R;
            r = slope_limiter(char_0_R, char_0_C);
            char_0_C = dxc * r * char_0_C;


            /* Project back */

            dw_0_C = char_0_C + char_1_C + char_4_C;
            dw_1_C = cs/dC*(char_4_C - char_0_C);
            dw_4_C = cs2*(char_0_C + char_4_C);




            real a_m, b_m, a_0, b_0, a_p, b_p;




            /* Interpolate to edges and limit once more */

            w_0_L = .5*(dC + dL) - (dw_0_C  - dw_0_L)/6.;
            w_0_R = .5*(dC + dR) - (dw_0_R  - dw_0_C)/6.;

            w_1_L = .5*(uC + uL) - (dw_1_C  - dw_1_L)/6.;
            w_1_R = .5*(uC + uR) - (dw_1_R  - dw_1_C)/6.;

            w_2_L = .5*(uC2 + uL2) - (dw_2_C  - dw_2_L)/6.;
            w_2_R = .5*(uC2 + uR2) - (dw_2_R  - dw_2_C)/6.;

            w_3_L = .5*(uC3 + uL3) - (dw_3_C  - dw_3_L)/6.;
            w_3_R = .5*(uC3 + uR3) - (dw_3_R  - dw_3_C)/6.;

            w_4_L = .5*(pC + pL) - (dw_4_C  - dw_4_L)/6.;
            w_4_R = .5*(pC + pR) - (dw_4_R  - dw_4_C)/6.;

            slope_limiter_ppm(dL,dC,dR, &w_0_L,&w_0_R);
            slope_limiter_ppm(uL,uC,uR, &w_1_L,&w_1_R);
            slope_limiter_ppm(uL2,uC2,uR2, &w_2_L,&w_2_R);
            slope_limiter_ppm(uL3,uC3,uR3, &w_3_L,&w_3_R);
            slope_limiter_ppm(pL,pC,pR, &w_4_L,&w_4_R);

            /* Evolve */

            d_6  = 6.*(dC  - .5*(w_0_L + w_0_R));
            u_6  = 6.*(uC  - .5*(w_1_L + w_1_R));
            u2_6 = 6.*(uC2 - .5*(w_2_L + w_2_R));
            u3_6 = 6.*(uC3 - .5*(w_3_L + w_3_R));
            p_6  = 6.*(pC  - .5*(w_4_L + w_4_R));

            /* Characteristic wave speeds */




            a_m = (uC-cs)*dtdx;
            a_0 = uC*dtdx;
            a_p = (uC+cs)*dtdx;

            a_m = MIN2(a_m, 0.);
            a_0 = MIN2(a_0, 0.);
            a_p = MIN2(a_p, 0.);

            b_m = MAX2(a_m, 0.);
            b_0 = MAX2(a_0, 0.);
            b_p = MAX2(a_p, 0.);

            dw_0_L = w_0_R - w_0_L;
            dw_1_L = w_1_R - w_1_L;
            dw_2_L = w_2_R - w_2_L;
            dw_3_L = w_3_R - w_3_L;
            dw_4_L = w_4_R - w_4_L;

            /* Evolve over domain of dependence */

            w_0_L = w_0_L - .5 * a_m  * ( dw_0_L + (1. + 2./3*a_m )*d_6);
            w_0_R = w_0_R - .5 * b_p  * ( dw_0_L - (1. - 2./3*b_p )*d_6);

            w_1_L = w_1_L - .5 * a_m  * ( dw_1_L + (1. + 2./3*a_m )*u_6);
            w_1_R = w_1_R - .5 * b_p  * ( dw_1_L - (1. - 2./3*b_p )*u_6);

            w_2_L = w_2_L - .5 * a_m  * ( dw_2_L + (1. + 2./3*a_m )*u2_6);
            w_2_R = w_2_R - .5 * b_p  * ( dw_2_L - (1. - 2./3*b_p )*u2_6);

            w_3_L = w_3_L - .5 * a_m  * ( dw_3_L + (1. + 2./3*a_m )*u3_6);
            w_3_R = w_3_R - .5 * b_p  * ( dw_3_L - (1. - 2./3*b_p )*u3_6);

            w_4_L = w_4_L - .5 * a_m  * ( dw_4_L - w_4_L + (1. + 2./3*a_m )*p_6);
            w_4_R = w_4_R - .5 * b_p  * ( dw_4_L - w_4_L - (1. - 2./3*b_p )*p_6);

            /* Correct for waves which do not reach interface */
            /* Left edge */


            if (uC < 0) {
            	w_0_R += .5*a_m*( dw_0_L - dw_4_L/cs2 + (d_6 - p_6/cs2)*(1 + 2./3 * a_m));
            	w_0_R -= .5*a_0*( dw_0_L - dw_4_L/cs2 + (d_6 - p_6/cs2)*(1 + 2./3 * a_0));

            	w_2_R += .5*a_m*( dw_2_L + u2_6*(1 + 2./3 * a_m));
            	w_2_R -= .5*a_0*( dw_2_L + u2_6*(1 + 2./3 * a_0));
            	w_3_R += .5*a_m*( dw_3_L + u3_6*(1 + 2./3 * a_m));
            	w_3_R -= .5*a_0*( dw_3_L + u3_6*(1 + 2./3 * a_0));

            }
            if (uC + cs < 0) {
            	w_0_R += .5*a_m*( .5*(dC/cs*dw_1_L + dw_4_L/cs2) +  .5*(dC/cs * u1_6 + p_6/cs2)*(1 + 2./3 * a_p));
            	w_0_R -= .5*a_p*( .5*(dC/cs*dw_1_L + dw_4_L/cs2) +  .5*(dC/cs * u1_6 + p_6/cs2)*(1 + 2./3 * a_p));

            	w_1_R += .5*a_m*( .5*(dw_1_L + dw_4_L/(cs*dC)) +  .5*( u1_6 + p_6/(cs*dC))*(1 + 2./3 * a_p));
            	w_1_R -= .5*a_p*( .5*(dw_1_L + dw_4_L/(cs*dC)) +  .5*( u1_6 + p_6/(cs*dC))*(1 + 2./3 * a_p));

            	w_4_R += .5*a_m*( .5*(dC*cs*dw_1_L + dw_4_L) +  .5*( dC*cs*u1_6 + p_6)*(1 + 2./3 * a_p));
            	w_4_R -= .5*a_p*( .5*(dC*cs*dw_1_L + dw_4_L) +  .5*( dC*cs*u1_6 + p_6)*(1 + 2./3 * a_p));

            }
            if (uC > 0) {
            	w_0_L += .5*b_p*( dw_0_L - dw_4_L/cs2 - (d_6 - p_6/cs2)*(1 - 2./3 * b_p));
            	w_0_L -= .5*b_0*( dw_0_L - dw_4_L/cs2 + (d_6 - p_6/cs2)*(1 + 2./3 * b_0));

            	w_2_L += .5*b_p*( dw_2_L - u2_6*(1 - 2./3 * b_p));
            	w_2_L -= .5*b_0*( dw_2_L + u2_6*(1 + 2./3 * b_0));
            	w_3_L += .5*b_p*( dw_3_L - u3_6*(1 - 2./3 * b_p));
            	w_3_L -= .5*b_0*( dw_3_L + u3_6*(1 + 2./3 * b_0));

            }
            if (uC - cs > 0) {
            	w_0_L += .5*b_p*( .5*(-dC/cs*dw_1_L + dw_4_L/cs2) -  .5*(-dC/cs * u1_6 + p_6/cs2)*(1 - 2./3 * b_p));
            	w_0_L -= .5*b_m*( .5*(-dC/cs*dw_1_L + dw_4_L/cs2) -  .5*(-dC/cs * u1_6 + p_6/cs2)*(1 - 2./3 * b_m));

            	w_1_L += .5*b_p*( .5*(dw_1_L - dw_4_L/(cs*dC)) -  .5*( u1_6 - p_6/(cs*dC))*(1 - 2./3 * b_p));
            	w_1_L -= .5*b_m*( .5*(dw_1_L - dw_4_L/(cs*dC)) -  .5*( u1_6 - p_6/(cs*dC))*(1 - 2./3 * b_m));

            	w_4_L += .5*b_p*( .5*(-dC*cs*dw_1_L + dw_4_L) -  .5*( -dC*cs*u1_6 + p_6)*(1 - 2./3 * b_p));
            	w_4_L -= .5*b_m*( .5*(-dC*cs*dw_1_L + dw_4_L) -  .5*( -dC*cs*u1_6 + p_6)*(1 - 2./3 * b_m));


            }

            /* Enforce floors */
            if (w_0_L < PRESSUREFLOOR) w_0_L = PRESSUREFLOOR;
            if (w_0_R < PRESSUREFLOOR) w_0_R = PRESSUREFLOOR;
            if (w_4_L < PRESSUREFLOOR) w_4_L = PRESSUREFLOOR;
            if (w_4_R < PRESSUREFLOOR) w_4_R = PRESSUREFLOOR;



            /* Convert to conservatives and place in final arrays */

            ke = .5*w_0_L*(w_1_L*w_1_L + w_2_L*w_2_L + w_3_L*w_3_L);
            w_1_L *= w_0_L;
            w_2_L *= w_0_L;
            w_3_L *= w_0_L;
            eL = (w_4_L/g1 + ke);


            ke = .5*w_0_R*(w_1_R*w_1_R + w_2_R*w_2_R + w_3_R*w_3_R);
            w_1_R *= w_0_R;
            w_2_R *= w_0_R;
            w_3_R *= w_0_R;
            eR = (w_4_R/g1 + ke);

            UL[indx + 0*ntot]  = w_0_R;
            UR[indxm + 0*ntot] = w_0_L;

            UL[indx + dir1*ntot]  = w_1_R;
            UR[indxm + dir1*ntot] = w_1_L;

            UL[indx + dir2*ntot]  = w_2_R;
            UR[indxm + dir2*ntot] = w_2_L;

            UL[indx + dir3*ntot]  = w_3_R;
            UR[indxm + dir3*ntot] = w_3_L;

            UL[indx + 4*ntot]  = eR;
            UR[indxm + 4*ntot] = eL;




            /* Scalars */
            for(n=5;n<nf;n++) {
                sR = cons[indxp + n*ntot]/dRi;
                sL = cons[indxm + n*ntot]/dLi;
                sC = cons[indx + n*ntot]/dCi;
                slopeR = (sR-sC) * 2./(dxc + dxp);
                slopeL = (sC-sL) * 2./(dxc + dxm);
                r = slope_limiter(slopeL,slopeR);
                sL = sC - .5*dxc * r * slopeR;
                sR = sC + .5*dxc * r * slopeR;

                sL *= dLf;
                sR *= dRf;

                fl = sL*uLf;
                fr = sR*uRf;
                sL += dtdx*(fl-fr);
                sR += dtdx*(fl-fr);

                UL[indx + n*ntot] = sR;
                UR[indxm + n*ntot] = sL;

            }
        }
    }

    return;
}
#endif
