#include "defs.h"
#define H5_USE_16_API
#include <hdf5.h>


#define H5CheckError(CMD) {                                          \
 herr_t status= CMD;                                 \
 if(status < 0) {                                              \
   printf("HDF5 failure %s:%d\n",__FILE__,__LINE__);           \
   exit(0); \
 }                                                                 \
}


#ifdef ISFLOAT
#define H5TYPE H5T_NATIVE_FLOAT
#else
#define H5TYPE H5T_NATIVE_DOUBLE
#endif




void write_hdf5_real(real *data, hsize_t *dims, int ndims, hid_t group_path, const char *name);
void read_hdf5_real(real *data, hid_t group_path, const char *name);

void output(int step, GridCons *grid, FluxCons *fluxes, Parameters *params) {
    int i;


    char fname[512];
    char scalarname[512];

    sprintf(fname, "%s_%d.h5",params->outputname,step);


    real *xm1 = &grid->xm1[-NGHX1];
    real *xm2 = &grid->xm2[-NGHX2];


    real *rho       = &grid->cons[0*grid->ntot-grid->offset];
    real *mx1       = &grid->cons[1*grid->ntot-grid->offset];
    real *mx2       = &grid->cons[2*grid->ntot-grid->offset];
    real *mx3       = &grid->cons[3*grid->ntot-grid->offset];
    real *energy    = &grid->cons[4*grid->ntot-grid->offset];
    real *intenergy = &grid->intenergy[-grid->offset];



    hid_t file_id, data_id; 

    hsize_t dims3[2];
    hsize_t dims_single[1] = {1}; 
    hsize_t dims_x1[1], dims_x2[1];
    real time[1] = {grid->time};
    real gamma[1] = {params->gamma};

    dims_x1[0] = grid->size_x[0] + 1;
    dims_x2[0] = grid->size_x[1] + 1;
    dims3[0] = grid->size_x[0];
    dims3[1] = grid->size_x[1];
    //for(i=0;i<3;i++) dims3[i] = grid->nx[i];

    printf("Outputting Results to %s...\n",fname);
    file_id = H5Fcreate(fname, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    data_id = H5Gcreate(file_id,"Data",0);
    //params_id = H5Gcreate(file_id"/Parameters",0);







    write_hdf5_real(time,dims_single,1,data_id,"time");
    write_hdf5_real(gamma,dims_single,1,data_id,"Gamma");
    write_hdf5_real(xm1,dims_x1,1,data_id,"xm1");
    write_hdf5_real(xm2,dims_x2,1,data_id,"xm2");
    write_hdf5_real(rho, dims3, 2, data_id, "Density");
    write_hdf5_real(mx1, dims3, 2, data_id, "Mx1");
    write_hdf5_real(mx2, dims3, 2, data_id, "Mx2");
    write_hdf5_real(mx3, dims3, 2, data_id, "Mx3");
    write_hdf5_real(energy, dims3, 2, data_id, "Energy");
    write_hdf5_real(intenergy, dims3, 2, data_id, "InternalEnergy");

    for(i=5;i<grid->nf;i++) {
    	sprintf(scalarname, "%s%d","Scalar",i-4);
    	write_hdf5_real(&grid->cons[i*grid->ntot - grid->offset],dims3,2,data_id,scalarname);
    }

    //status = H5Gclose(params_id);
    //if (status < 0) printf("HDF5 error\n");
    H5CheckError(   H5Gclose(data_id)    )

    H5CheckError(   H5Fclose(file_id)    )

    return;

}

void read_restart(const char *fname, GridCons *grid, FluxCons *fluxes, Parameters *params) {
	int i;
    char scalarname[512];


    real *xm1 = &grid->xm1[-NGHX1];
    real *xm2 = &grid->xm2[-NGHX2];


    real *rho       = &grid->cons[0*grid->ntot-grid->offset];
    real *mx1       = &grid->cons[1*grid->ntot-grid->offset];
    real *mx2       = &grid->cons[2*grid->ntot-grid->offset];
    real *mx3       = &grid->cons[3*grid->ntot-grid->offset];
    real *energy    = &grid->cons[4*grid->ntot-grid->offset];
    real *intenergy = &grid->intenergy[-grid->offset];



    hid_t file_id, data_id;


    printf("Reading restart file %s...\n",fname);
    file_id = H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);
    data_id = H5Gopen(file_id,"Data");
    //params_id = H5Gcreate(file_id"/Parameters",0);







    read_hdf5_real(&grid->time,data_id,"time");
    read_hdf5_real(&params->gamma,data_id,"Gamma");
    read_hdf5_real(xm1,data_id,"xm1");
    read_hdf5_real(xm2,data_id,"xm2");
    read_hdf5_real(rho,data_id, "Density");
    read_hdf5_real(mx1,data_id, "Mx1");
    read_hdf5_real(mx2, data_id, "Mx2");
    read_hdf5_real(mx3, data_id, "Mx3");
    read_hdf5_real(energy, data_id, "Energy");
    read_hdf5_real(intenergy, data_id, "InternalEnergy");

    for(i=5;i<grid->nf;i++) {
    	sprintf(scalarname, "%s%d","Scalar",i-4);
    	read_hdf5_real(&grid->cons[i*grid->ntot - grid->offset],
    			data_id,scalarname);
    }

    //status = H5Gclose(params_id);
    //if (status < 0) printf("HDF5 error\n");
    H5CheckError(   H5Gclose(data_id)    )
    H5CheckError(   H5Fclose(file_id)    )
    return;

}


void write_hdf5_real(real *data, hsize_t *dims, int ndims, hid_t group_path, const char *name) {
  hid_t dspc_id, dset_id;

  dspc_id = H5Screate_simple(ndims,dims,NULL);
  dset_id = H5Dcreate(group_path,name,H5TYPE,dspc_id,H5P_DEFAULT);

  H5CheckError(   H5Dwrite(dset_id,H5TYPE,H5S_ALL,H5S_ALL,H5P_DEFAULT,data)   )
  H5CheckError(   H5Sclose(dspc_id)   )
  H5CheckError(   H5Dclose(dset_id)   )

  return;
}

void read_hdf5_real(real *data, hid_t group_path, const char *name) {
  hid_t dset_id;
  dset_id = H5Dopen(group_path,name);

  H5CheckError(   H5Dread(dset_id,H5TYPE,H5S_ALL,H5S_ALL,H5P_DEFAULT,data)   )
  H5CheckError(   H5Dclose(dset_id)   )
  return;
}
