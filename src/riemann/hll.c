

void intercell_flux(double *FL, double *FR, double *SL, double *SR, double *UL, double *UR, int n, double *Fint) {

    int i;

    for(i=0;i<n;i++) {
        if (SL >= 0) {
/* Supersonic left moving wave 
 */
            Fint[i] = FL[i];

        }
        else {
            if (SR <= 0 ) {
/* Supersonic right moving wave 
 */
                Fint[i] = FR[i];

            }
            else {
/* Subsonic case
 */
                Fint[i] = SR[i]*FL[i] - SL[i] * FR[i] + SL[i]*SR[i]*(UR[i]-UL[i]);
                Fint[i] /= (SR[i] - SL[i]);
            }

    }
    return;
}

