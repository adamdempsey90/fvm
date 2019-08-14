
#define NGHX 2

template<typename real, size_t N>
class Patch {
    /* Patches contain all information for 
     * evolving the fluid equations forward in 
     * time.
     * They are a fixed size determined by the problem/hardware
     * e.g. in 2D they can be 32^2 to fill 1 threadblock (1024 threads)
     * and in 3D they can be 8^3 (512 threads) or 10^3 (1000 threads).
     * These include the ghost cells.
     * 
     * Communication with other patches is done through the Domain class.
     * The Domain class also sets the boundary/initial conditions of
     * a patch.
     */

    private:
        std::array<real,N+1> xm1;
        std::array<real,N+1> xm2;
        std::array<real,N> xc1;
        std::array<real,N> xc2;
        std::array<real,N> dx1;
        std::array<real,N> dx2;

        std::array<real,N*N*5> cons;
        std::array<real,N*N*5> prim;
        std::array<real,N*N> intenergy;
        std::array<bool,N*N> mask;

        real *scalars;

        real *d_cons, *d_prim, *d_xm1, *d_xm2, *d_xc1;
        real *d_xc2, *d_dx1, *d_dx2, *d_intenergy;
        int *d_mask;


        real xmin, xmax, ymin, ymax;
    public:
        Patch2D();
        ~Patch2D();
        void fill_ghosts(); 
        void set_domain(real xmin, real xmax, real ymin, real ymax);
        void fill_arrays();
        void cons_to_prim();
        void prim_to_cons();

/* Copy data to/from device */
        void host_to_dev();
        void dev_to_host();

        real evolve(); // returns new dt 

        void integrate();



};
