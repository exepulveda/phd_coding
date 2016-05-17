#include <H5Cpp.h>
#include <cstring>

#include <armadillo>


using namespace std;
using namespace H5;
using namespace arma;

int dimensions_hdf5(const std::string name, const std::string dsname);
irowvec hdf5_shape(const std::string name, const std::string dsname);

void load_hdf5_vector(frowvec & x, const string name, const string dsname);
//void load_hdf5_vector_uint(urowvec& x, const string name, const string dsname);
//void load_hdf5_vector_int(irowvec& x, const string name, const string dsname);

//void load_hdf5_matrix(void *ptr,int n, int m,const string name, const string dsname, const PredType ptype);
//void load_hdf5_matrix_double(dmat& x, const string name, const string dsname);
//void load_hdf5_matrix_uint(umat& x, const string name, const string dsname);
void load_hdf5_matrix(Mat<int>& x, const string name, const string dsname);
void load_hdf5_matrix(fmat& x, const string name, const string dsname);
//void load_hdf5_matrix_sint(smat& x, const string name, const string dsname);

// void load_hdf5_matrix_int_eigen(MatrixXi & x, const string name, const string dsname);

template<class T>
void load_hdf5_matrix(T & x, const string name, const string dsname) {
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
    
    dataset.read( (void *) x.memptr(), PredType::NATIVE_INT, DataSpace::ALL, DataSpace::ALL);
    
    inplace_trans(x);
}
