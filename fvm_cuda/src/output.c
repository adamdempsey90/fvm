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


void snapshot_1d(char *name, GridCons *grid, Parameters *params) {

    int size_x1 = grid->size_x1;
    int size_x2 = grid->size_x2;
    int size_x3 = grid->size_x3;
    int ntot = grid->ntot;
    int nf = grid->nf;
    real *xm1 = &grid->xm1[-NGHX1];
    real *xm2 = &grid->xm2[-NGHX2];
    real *xm3 = &grid->xm3[-NGHX3];
    real *prim = &grid->prim[-grid->offset];

    char fname[512];
    sprintf(fname, "%s/%s.h5", params->outputdir,name);
#ifndef SILENT
	printf("Outputting 1D snapshot at t=%.2e to %s...\n",grid->time,fname);
#endif
	int i;
	char scalarname[512];
	hid_t file_id, data_id;


	hsize_t dims3[1];
	hsize_t dims_single[1] = {1};
	hsize_t dims_x1[1],dims_x2[1],dims_x3[1];
	real time[1] = {grid->time};
	real gamma[1] = {params->gamma};

	dims_x1[0] = size_x1 + 1;
	dims_x2[0] = size_x2 + 1;
	dims_x3[0] = size_x3 + 1;

	dims3[0] = size_x1;

	file_id = H5Fcreate(fname, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
	data_id = H5Gcreate(file_id,"Data",0);


	write_hdf5_real(time,dims_single,1,data_id,"time");
	write_hdf5_real(gamma,dims_single,1,data_id,"Gamma");
	write_hdf5_real(xm1,dims_x1,1,data_id,"xm1");
	write_hdf5_real(xm2,dims_x2,1,data_id,"xm2");
	write_hdf5_real(xm3,dims_x3,1,data_id,"xm3");

	write_hdf5_real(&prim[0*ntot], dims3, 1, data_id, "Density");
	write_hdf5_real(&prim[1*ntot], dims3, 1, data_id, "Vx1");
	write_hdf5_real(&prim[2*ntot], dims3, 1, data_id, "Vx2");
	write_hdf5_real(&prim[3*ntot], dims3, 1, data_id, "Vx3");
	write_hdf5_real(&prim[4*ntot], dims3, 1, data_id, "Pressure");

	for(i=5;i<nf;i++) {
		sprintf(scalarname, "Scalar%d",i-4);
		write_hdf5_real(&prim[i*ntot],dims3,1,data_id,scalarname);
	}


	H5CheckError(   H5Gclose(data_id)    )
	H5CheckError(   H5Fclose(file_id)    )

	return;
}

void snapshot_2d(char *name, GridCons *grid, Parameters *params) {

    int size_x1 = grid->size_x1;
    int size_x2 = grid->size_x2;
    int size_x3 = grid->size_x3;
    int ntot = grid->ntot;
    int nf = grid->nf;
    real *xm1 = &grid->xm1[-NGHX1];
    real *xm2 = &grid->xm2[-NGHX2];
    real *xm3 = &grid->xm3[-NGHX3];
    real *prim = &grid->prim[-grid->offset];
    char fname[512];
    sprintf(fname, "%s/%s.h5", params->outputdir,name);
#ifndef SILENT
	printf("Outputting 2D snapshot at t=%.2e to %s...\n",grid->time,fname);
#endif

	int i;
	char scalarname[512];
	hid_t file_id, data_id;

	hsize_t dims3[2];
	hsize_t dims_single[1] = {1};
	hsize_t dims_x1[1],dims_x2[1],dims_x3[1];
	real time[1] = {grid->time};
	real gamma[1] = {params->gamma};

	dims_x1[0] = size_x1 + 1;
	dims_x2[0] = size_x2 + 1;
	dims_x3[0] = size_x3 + 1;

	dims3[0] = size_x1;
	dims3[1] = size_x2;


	file_id = H5Fcreate(fname, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
	data_id = H5Gcreate(file_id,"Data",0);
	write_hdf5_real(time,dims_single,1,data_id,"time");
	write_hdf5_real(gamma,dims_single,1,data_id,"Gamma");
	write_hdf5_real(xm1,dims_x1,1,data_id,"xm1");
	write_hdf5_real(xm2,dims_x2,1,data_id,"xm2");
	write_hdf5_real(xm3,dims_x3,1,data_id,"xm3");

	write_hdf5_real(&prim[0*ntot], dims3, 2, data_id, "Density");
	write_hdf5_real(&prim[1*ntot], dims3, 2, data_id, "Vx1");
	write_hdf5_real(&prim[2*ntot], dims3, 2, data_id, "Vx2");
	write_hdf5_real(&prim[3*ntot], dims3, 2, data_id, "Vx3");
	write_hdf5_real(&prim[4*ntot], dims3, 2, data_id, "Pressure");

	for(i=5;i<nf;i++) {
		sprintf(scalarname, "Scalar%d",i-4);
		write_hdf5_real(&prim[i*ntot],dims3,2,data_id,scalarname);
	}



	H5CheckError(   H5Gclose(data_id)    )
	H5CheckError(   H5Fclose(file_id)    )
	return;
}

void snapshot_3d(char *name, GridCons *grid, Parameters *params) {

    int size_x1 = grid->size_x1;
    int size_x2 = grid->size_x2;
    int size_x3 = grid->size_x3;
    int ntot = grid->ntot;
    int nf = grid->nf;
    real *xm1 = &grid->xm1[-NGHX1];
    real *xm2 = &grid->xm2[-NGHX2];
    real *xm3 = &grid->xm3[-NGHX3];
    real *prim = &grid->prim[-grid->offset];


    char fname[512];
    sprintf(fname, "%s/%s.h5", params->outputdir,name);
#ifndef SILENT
	printf("Outputting 3D snapshot at t=%.2e to %s...\n",grid->time,fname);
#endif

	int i;
	char scalarname[512];
	hid_t file_id, data_id;

	hsize_t dims3[3];
	hsize_t dims_single[1] = {1};
	hsize_t dims_x1[1], dims_x2[1],dims_x3[1];
	real time[1] = {grid->time};
	real gamma[1] = {params->gamma};

	dims_x1[0] = size_x1 + 1;
	dims_x2[0] = size_x2 + 1;
	dims_x3[0] = size_x3 + 1;

	dims3[0] = size_x1;
	dims3[1] = size_x2;
	dims3[2] = size_x3;

	file_id = H5Fcreate(fname, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
	data_id = H5Gcreate(file_id,"Data",0);


	write_hdf5_real(time,dims_single,1,data_id,"time");
	write_hdf5_real(gamma,dims_single,1,data_id,"Gamma");
	write_hdf5_real(xm1,dims_x1,1,data_id,"xm1");
	write_hdf5_real(xm2,dims_x2,1,data_id,"xm2");
	write_hdf5_real(xm3,dims_x3,1,data_id,"xm3");
	write_hdf5_real(&prim[0*ntot], dims3, 3, data_id, "Density");
	write_hdf5_real(&prim[1*ntot], dims3, 3, data_id, "Vx1");
	write_hdf5_real(&prim[2*ntot], dims3, 3, data_id, "Vx2");
	write_hdf5_real(&prim[3*ntot], dims3, 3, data_id, "Vx3");
	write_hdf5_real(&prim[4*ntot], dims3, 3, data_id, "Pressure");

	for(i=5;i<nf;i++) {
		sprintf(scalarname, "Scalar%d",i-4);
		write_hdf5_real(&prim[i*ntot],dims3,3,data_id,scalarname);
	}


	H5CheckError(   H5Gclose(data_id)    )
	H5CheckError(   H5Fclose(file_id)    )

	return;
}



//void output(int step, GridCons *grid, Parameters *params) {
//
//
//    char fname[512];
//
//    sprintf(fname, "%s_%d.h5",params->outputname,step);
//
//#ifdef DIMS3
//    		/* 3d output */
//	output_3d(fname,&grid->xm1[-NGHX1],&grid->xm2[-NGHX2],&grid->xm3[-NGHX3],&grid->prim[-grid->offset],
//			grid->time,params->gamma,grid->ntot,grid->nf,grid->size_x1,grid->size_x2,grid->size_x3);
//#else
//#ifdef DIMS2
//	/* 2d output */
//	output_2d(fname,&grid->xm1[-NGHX1],&grid->xm2[-NGHX2],&grid->prim[-grid->offset],
//						grid->time,params->gamma,grid->ntot,grid->nf,grid->size_x1,grid->size_x2);
//#else
//	/* 1d output */
//	output_1d(fname,&grid->xm1[-NGHX1],&grid->prim[-grid->offset],
//									grid->time,params->gamma,grid->ntot,grid->nf,grid->size_x1);
//#endif
//#endif
//
//
//
//    return;
//
//}


void read_restart(const char *fname, GridCons *grid, Parameters *params) {
	int i,n,ntot,offset;
    char scalarname[512];
    ntot = grid->ntot;
    offset = grid->offset;

    real *xm1 = &grid->xm1[-NGHX1];
    real *xm2 = &grid->xm2[-NGHX2];
    real *xm3 = &grid->xm3[-NGHX3];

    real *vx1       = &grid->prim[1*ntot-offset];
    real *vx2       = &grid->prim[2*ntot-offset];
    real *vx3       = &grid->prim[3*ntot-offset];


    real *rho       = &grid->prim[0*ntot-offset];
    real *pres    = &grid->prim[4*ntot-offset];
    real *intenergy = &grid->intenergy[-offset];



    hid_t file_id, data_id;

#ifndef SILENT
    printf("Reading restart file %s...\n",fname);
#endif
    file_id = H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);
    data_id = H5Gopen(file_id,"Data");
    //params_id = H5Gcreate(file_id"/Parameters",0);







    read_hdf5_real(&grid->time,data_id,"time");
    read_hdf5_real(&params->gamma,data_id,"Gamma");
    read_hdf5_real(xm1,data_id,"xm1");
    read_hdf5_real(xm2,data_id,"xm2");
    read_hdf5_real(xm3,data_id,"xm3");

    read_hdf5_real(vx1,data_id,"Vx1");
    read_hdf5_real(vx2,data_id,"Vx2");
    read_hdf5_real(vx3,data_id,"Vx3");

    read_hdf5_real(rho,data_id, "Density");
    read_hdf5_real(pres,data_id, "Pressure");


    for(i=5;i<grid->nf;i++) {
    	sprintf(scalarname, "Scalar%d",i-4);
    	read_hdf5_real(&grid->prim[i*ntot-offset],
    			data_id,scalarname);
    }

    //status = H5Gclose(params_id);
    //if (status < 0) printf("HDF5 error\n");
    H5CheckError(   H5Gclose(data_id)    )
    H5CheckError(   H5Fclose(file_id)    )


    return;

}

void snapshot(char *name, GridCons *grid, Parameters *params) {
#ifdef DIMS3
    snapshot_3d(name,grid,params);
#else
#ifdef DIMS2
    snapshot_2d(name,grid,params);
#else
    snapshot_1d(name,grid,params);
#endif
#endif
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



