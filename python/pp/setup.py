from numpy.distutils.core import Extension

ext1 = Extension(name = 'ppr', sources = ['ppr.f'])


# setup fortran 90 extension
#---------------------------------------------------------------------------  


# call setup
#--------------------------------------------------------------------------
if __name__ == "__main__":
    from numpy.distutils.core import setup


    setup( 
        name = 'ppreg',
        version = '0.1',        
        description='Python interface of Projection Pursuit Regression',
        author='Exequiel Sepulveda',
        author_email='exequiel.sepulveda@gmail.com',
        py_modules = ["ppreg_test","ppreg"],
        ext_modules = [ext1]
    )  
