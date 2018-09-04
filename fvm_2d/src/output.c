#include "defs.h"
#define H5_USE_16_API
#include <hdf5.h>

#ifdef ISFLOAT
#define H5TYPE H5T_NATIVE_FLOAT
#else
#define H5TYPE H5T_NATIVE_DOUBLE
#endif




void write_hdf5_real(real *data, hsize_t *dims, int ndims, hid_t group_path, char *name);
void output(int step, GridCons *grid, FluxCons *fluxes, Parameters *params) {
    int i,j,k,indx;
    int nx1,nx2,nx3;
    int size_x1,size_x12; 
    nx1 = grid->nx[0];
    size_x1 = grid->size_x[0];

    char fname[512];

    sprintf(fname, "%s_%d.h5",params->outputname,step);


    real *xm1 = grid->xm1;
    real *xm2 = grid->xm2;


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

    dims_x1[0] = grid->nx[0]+1;
    dims_x2[0] = grid->nx[1]+1;
    dims3[0] = grid->size_x[0];
    dims3[1] = grid->size_x[1];
    //for(i=0;i<3;i++) dims3[i] = grid->nx[i];

    herr_t status;
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

    //status = H5Gclose(params_id);
    //if (status < 0) printf("HDF5 error\n");
    status = H5Gclose(data_id);
    if (status < 0) printf("HDF5 error\n");
    status = H5Fclose(file_id);
    if (status < 0) printf("HDF5 error\n");
    return;

}

void write_hdf5_real(real *data, hsize_t *dims, int ndims, hid_t group_path, char *name) {
  hid_t dspc_id, dset_id;
  herr_t status;

  dspc_id = H5Screate_simple(ndims,dims,NULL);
  dset_id = H5Dcreate(group_path,name,H5TYPE,dspc_id,H5P_DEFAULT);

  status = H5Dwrite(dset_id,H5TYPE,H5S_ALL,H5S_ALL,H5P_DEFAULT,data);
	if (status < 0) printf("HDF5 error\n");
  status = H5Sclose(dspc_id);
	if (status < 0) printf("HDF5 error\n");
  status = H5Dclose(dset_id);
	if (status < 0) printf("HDF5 error\n");

  return;
}
