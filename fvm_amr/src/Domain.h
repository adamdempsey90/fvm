#include <bitset>
#include <map>
#include <vector>
#include "Patch.h"


using uint = unsigned int;
template <size_t N>
using morton = std::bitset<N>;

template <size_t N>
class compare {
    public:
        bool operator()(const morton<N> &a, const morton<N> &b) const 
        {
            return a.to_ullong() < b.to_ullong();
        }

};

template <size_t N, typename mtype>
using Mesh = std::map<morton<N>, mtype *, compare<N>>;

template <size_t N>
uint morton_level(morton<N> key);



template <size_t N>
class Domain {
    /* The Domain class contains the hierarchy of patches and 
     * facilitates communication between the patches.
     * It handles global boundary conditions and fills each patch's
     * ghost zones or active zones when the patch is created.
     * AMR is also controlled by the Domain class. 
     *
     * When creating a Domain class, the maximum level of refinement 
     * must be set via N with N = 2*(max level-1)
     */

    public:

        real xmin, xmax, ymin, ymax; // Domain extent
        uint max_level = 1 + N/2;

        /* Patches are stored with std::map via their Z-order index.
         * The map container is a RB BST. 
         */
        Mesh<N,Patch> patches[1+N/2]; 
    
        std::queue<Patch> refinelist;

    public:
        Domain();
        ~Domain();
        
        /* Functions to handle patch business */
        void get_patch_level();
        //void get_patch_extent();
        void get_neighbors();

        /* Functions to fill patches */
        void fill_patch();
        void fill_ghosts();

        /* AMR functions */
        void refinement();
        void criterion();
        void refine_patch();
        void derefine_patch();
        void enforce_balance();


        /* Sync functions */
        void exchange_fluxes();
        void set_dts();

        /* Output functions */
        void output_domain();
        void output_patch();
    


        

};
