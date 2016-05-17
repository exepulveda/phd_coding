#include <H5Cpp.h>
#include <assert.h>
#include "hdf5_utils.h"
#include <cstdio>
#include <iostream>

using namespace std;
using namespace H5;

int dimensions_hdf5(const std::string name, const std::string dsname) {
    bool load_okay = false;

    H5File file(name,H5F_ACC_RDONLY);
    DataSet dataset = file.openDataSet(dsname);
    DataSpace dataspace = dataset.getSpace();
    int rank = dataspace.getSimpleExtentNdims();
    
    
    dataspace.close();
    file.close();

    return rank;
}

irowvec hdf5_shape(const std::string name, const std::string dsname) {
    bool load_okay = false;

    H5File file(name,H5F_ACC_RDONLY);
    DataSet dataset = file.openDataSet(dsname);
    DataSpace dataspace = dataset.getSpace();
    int rank = dataspace.getSimpleExtentNdims();
    hsize_t *shape = new hsize_t[rank];
    dataspace.getSimpleExtentDims( shape, NULL);
    
    irowvec ret(rank);
    for (int i=0;i<rank;i++) {
        ret[i] = shape[i];
        printf("shape[%d]=%d\n",i,shape[i]);
    }
    
    dataspace.close();
    file.close();

    delete [] shape;
    
    return ret;
}


void load_hdf5_vector(frowvec& x, const std::string name, const std::string dsname) {
    // These may be necessary to store the error handler (if we need to).
    herr_t (*old_func)(hid_t, void*);
    void *old_client_data;


    bool load_okay = false;

    H5File file( name, H5F_ACC_RDONLY );

    DataSet dataset = file.openDataSet(dsname);
    H5T_class_t type_class = dataset.getTypeClass();
    DataSpace dataspace = dataset.getSpace();
    int rank = dataspace.getSimpleExtentNdims();
    
    assert(rank == 1);
    
    hsize_t dims_out[1];
    int ndims = dataspace.getSimpleExtentDims( dims_out, NULL);

    x.resize(dims_out[0]);
    
    dataset.read(void_ptr(x.memptr()), PredType::NATIVE_FLOAT, DataSpace::ALL, DataSpace::ALL);
    
    printf("DONE...load_hdf5_vector. [%d] elements\n",dims_out[0]);
}
/**

void load_hdf5_vector_uint(urowvec& x, const std::string name, const std::string dsname) {
    // These may be necessary to store the error handler (if we need to).
    herr_t (*old_func)(hid_t, void*);
    void *old_client_data;


    bool load_okay = false;

    H5File file( name, H5F_ACC_RDONLY );

    DataSet dataset = file.openDataSet(dsname);
    H5T_class_t type_class = dataset.getTypeClass();
    DataSpace dataspace = dataset.getSpace();
    int rank = dataspace.getSimpleExtentNdims();
    
    assert(rank == 1);
    
    hsize_t dims_out[1];
    int ndims = dataspace.getSimpleExtentDims( dims_out, NULL);

    x.resize(dims_out[0]);
    
    dataset.read(void_ptr(x.memptr()), PredType::NATIVE_INT, DataSpace::ALL, DataSpace::ALL);
}
void load_hdf5_matrix_uint(umat& x, const string name, const string dsname) {
    // These may be necessary to store the error handler (if we need to).
    herr_t (*old_func)(hid_t, void*);
    void *old_client_data;


    bool load_okay = false;

    H5File file( name, H5F_ACC_RDONLY );

    DataSet dataset = file.openDataSet(dsname);
    H5T_class_t type_class = dataset.getTypeClass();
    DataSpace dataspace = dataset.getSpace();
    int rank = dataspace.getSimpleExtentNdims();

    //cout << "rank: " << rank << endl;
    
    assert(rank == 2);
    
    hsize_t dims_out[2];
    int ndims = dataspace.getSimpleExtentDims( dims_out, NULL);

    //cout << "dim1: " << dims_out[0] << endl;
    //cout << "dim2: " << dims_out[1] << endl;

    //assert(n == dims_out[0]);
    //assert(m == dims_out[1]);
    
    x.resize(dims_out[0],dims_out[1]);

    dataset.read( (void *) x.memptr(), PredType::NATIVE_INT, DataSpace::ALL, DataSpace::ALL);
}

void load_hdf5_matrix_sint(smat& x, const string name, const string dsname) {
    // These may be necessary to store the error handler (if we need to).
    herr_t (*old_func)(hid_t, void*);
    void *old_client_data;


    bool load_okay = false;

    H5File file( name, H5F_ACC_RDONLY );

    DataSet dataset = file.openDataSet(dsname);
    H5T_class_t type_class = dataset.getTypeClass();
    DataSpace dataspace = dataset.getSpace();
    int rank = dataspace.getSimpleExtentNdims();

    //cout << "rank: " << rank << endl;
    
    assert(rank == 2);
    
    hsize_t dims_out[2];
    int ndims = dataspace.getSimpleExtentDims( dims_out, NULL);

    //cout << "dim1: " << dims_out[0] << endl;
    //cout << "dim2: " << dims_out[1] << endl;

    //assert(n == dims_out[0]);
    //assert(m == dims_out[1]);
    
    x.resize(dims_out[0],dims_out[1]);
    cout << "x resized"  << endl;

    cout << "reading to buffer..."  << endl;
    dataset.read( (void *) x.memptr(), PredType::STD_I32LE, DataSpace::ALL, DataSpace::ALL);
    cout << "DONE"  << endl;
}
**/

void load_hdf5_matrix(fmat& x, const string name, const string dsname) {
    // These may be necessary to store the error handler (if we need to).
    herr_t (*old_func)(hid_t, void*);
    void *old_client_data;


    bool load_okay = false;

    H5File file( name, H5F_ACC_RDONLY );

    DataSet dataset = file.openDataSet(dsname);
    H5T_class_t type_class = dataset.getTypeClass();
    DataSpace dataspace = dataset.getSpace();
    int rank = dataspace.getSimpleExtentNdims();
    
    DataType h5set_datatype = dataset.getDataType();

    //cout << "rank: " << rank << endl;
    
    assert(rank == 2);
    
    hsize_t dims_out[2];
    int ndims = dataspace.getSimpleExtentDims( dims_out, NULL);

    //cout << "dim1: " << dims_out[0] << endl;
    //cout << "dim2: " << dims_out[1] << endl;

    //assert(n == dims_out[0]);
    //assert(m == dims_out[1]);
    
    x.resize(dims_out[1],dims_out[0]);
    
    dataset.read( (void *) x.memptr(), PredType::NATIVE_FLOAT, DataSpace::ALL, DataSpace::ALL);
    
    //inplace_trans(x, "lowmem");
    inplace_trans(x);

    //cout << "loading DONE" << endl;
    
    dataset.close();
    file.close();
}
/**
void load_hdf5_matrix_double(dmat& x, const string name, const string dsname) {
    // These may be necessary to store the error handler (if we need to).
    herr_t (*old_func)(hid_t, void*);
    void *old_client_data;


    bool load_okay = false;

    H5File file( name, H5F_ACC_RDONLY );

    DataSet dataset = file.openDataSet(dsname);
    H5T_class_t type_class = dataset.getTypeClass();
    DataSpace dataspace = dataset.getSpace();
    int rank = dataspace.getSimpleExtentNdims();

    //cout << "rank: " << rank << endl;
    
    assert(rank == 2);
    
    hsize_t dims_out[2];
    int ndims = dataspace.getSimpleExtentDims( dims_out, NULL);

    cout << "dim1: " << dims_out[0] << endl;
    cout << "dim2: " << dims_out[1] << endl;

    //assert(n == dims_out[0]);
    //assert(m == dims_out[1]);
    
    x.resize(dims_out[0],dims_out[1]);

    dataset.read( (void *) x.memptr(), PredType::NATIVE_DOUBLE, DataSpace::ALL, DataSpace::ALL);
}


void load_hdf5_double_matrix(mat& x, const std::string name, const std::string dsname)
{
    // These may be necessary to store the error handler (if we need to).
    bool load_okay = false;

    hid_t fid = H5Fopen(name, H5F_ACC_RDONLY, H5P_DEFAULT);

    if(fid >= 0)
      {
      // MATLAB HDF5 dataset names are user-specified;
      // Octave tends to store the datasets in a group, with the actual dataset being referred to as "value".
      // So we will search for "dataset" and "value", and if those are not found we will take the first dataset we do find.
      hid_t dataset = H5Dopen(fid, dsname, H5P_DEFAULT);

      if(dataset >= 0)
      {
        hid_t filespace = H5Dget_space(dataset);

        // This must be <= 2 due to our search rules.
        const int ndims = H5Sget_simple_extent_ndims(filespace);

        hsize_t dims[2];
        const herr_t query_status = H5Sget_simple_extent_dims(filespace, dims, NULL);


        if(ndims == 1) { dims[1] = 1; }  // Vector case; fake second dimension (one column).

        x.set_size(dims[1], dims[0]);

        // Now we have to see what type is stored to figure out how to load it.
        hid_t datatype = H5Dget_type(dataset);

        hid_t read_status = arma_H5Dread(dataset, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, void_ptr(x.memptr()));


        // Now clean up.
        H5Tclose(datatype);
        H5Tclose(mat_type);
        H5Sclose(filespace);
        }

      H5Dclose(dataset);

      H5Fclose(fid);

}
**/

//void load_hdf5_matrix_int_eigen(MatrixXi & mat, const string name, const string dsname) {
    //// These may be necessary to store the error handler (if we need to).
    //herr_t (*old_func)(hid_t, void*);
    //void *old_client_data;


    //bool load_okay = false;

    //H5File file( name, H5F_ACC_RDONLY );

    //DataSet dataset = file.openDataSet(dsname);
    //H5T_class_t type_class = dataset.getTypeClass();
    //DataSpace dataspace = dataset.getSpace();
    //int rank = dataspace.getSimpleExtentNdims();

    ////cout << "rank: " << rank << endl;
    
    //assert(rank == 2);
    
    //hsize_t dims_out[2];
    //int ndims = dataspace.getSimpleExtentDims( dims_out, NULL);

    //cout << "dim1: " << dims_out[0] << endl;
    //cout << "dim2: " << dims_out[1] << endl;

    ////assert(n == dims_out[0]);
    ////assert(m == dims_out[1]);
    ////MatrixXi mat(dims_out[1],dims_out[0]);
    //cout << "matrix resized" << endl;
    //mat.resize(dims_out[1],dims_out[0]);

    //dataset.read(mat.derived().data(), PredType::NATIVE_INT, DataSpace::ALL, DataSpace::ALL);
    //cout << "matrix loaded" << endl;

    //mat.transposeInPlace();
    
    ////for (int i=0;i<dims_out[0];i++) {
    ////    for (int j=0;j<dims_out[1];j++) {
    ////        mat(i,j) = tmp(j,i);
    ////    }
    ////}
//}
