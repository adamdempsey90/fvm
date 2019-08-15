#include "Domain.h"


template<size_t N>
Domain<N>::Domain() {
    /* Initialize Domain */

    xmin = 0.;
    xmax = 1.;
    ymin = 0.;
    ymax = 0.;

    /* Create initial root patch */
    morton<N> key;
    key.reset();
    Patch *root = new Patch;

    /* Insert root into l=0 mesh */
    patches[0].insert({key,root});


}

template<size_t N>
void Domain<N>::refine_patch(morton<N> key, Patch *child) {
    /* Refine child patch */

    /* Deactivate patch on current level */

    /* Create siblings */


    /* Copy variables to siblings */


    /* Add patch and new siblings to mesh on finer level */




}

template<size_t N>
void Domain<N>::derefine_patch(morton<N> key, Patch *child) {
    /* Derefine child patch */

    /* Deactivate patch and siblings on current level */

    /* Reactivate patch on coarser level */

    /* Copy variables back to coarse patch */

}

template<size_t N>
void Domain<N>::refine(void) {
    /* Search for patches which need to be refined 
     * and adds them to the queue
     */

    int refine;
    morton<N> key;
    Mesh<N,Patch> *curr;
    Patch *patch;
    Patch *neighbors;
    /* Loop through levels in mesh starting at finest */
    for(int lvl=max_level-1; lvl > 0; lvl--) {
        curr = patches[lvl];
        
        for (Mesh<N,Patch>::iterator it=curr.begin(); it!=curr.end(); ++it) {
            key = it->first;
            patch = it->second;
            get_neighbors(key,patch,neighbors);
            refine = criterion(patch, neighbors);
            patch.refine = refine;
            if (refine == 1) {
                refinelist.push(patch);
            } 

        }
    }


    /* Get neighbor list */


    /* Check if a finer neighbor is already flagged */

    /* If it is then flag for refinement and move on */

    /* If not then check refinement criteria */

}
template<size_t N>
int Domain<N>::criterion(Patch *child, Patch *neighbors) {
    /* Check whether this patch needs to be refined or derefined */
    int ans;

    /* Evaluate derivatives from neighbors */

    return ans;
}


template<size_t N>
void Domain<N>::get_neighbors(morton<N> key, morton<N> *list, uint ext) {
    /* Retrieve neighbor keys */

    /* Start with siblings */


    /* Now cousins */

    /* Check which level the bit is flipped */

    /* Form cousin key to coarser level*/

    /* Check if this is active or not */

    /* If active then done */

    /* If not then form cousin key on our level */


}

template<size_t N>
void Domain<N>::fill_patch(morton<N> key, Patch *child) {
    /* Fill child patch arrays */

    /* Conservative splitting */

}

template<size_t N>
void Domain<N>::exchange_fluxes(morton<N> key, Patch *child) {
    /* At coarse/fine boundary have coarse cell use the sum of fine face fluxes*/


}

template<size_t N>
void Domain<N>::set_dts() {
    /* Set timesteps for mesh */
    
}

template<size_t N>
void Domain<N>::output_domain(void) {
    /* Output mesh */
}




