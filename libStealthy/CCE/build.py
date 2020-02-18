def build():
    '''
    This is called by buildExt.py to build the CCE submodule.
    '''
    from distutils.core import setup
    from distutils.extension import Extension

    ext_modules = [
        Extension(
            "libCCE"
            ,sources = [
                #'libCCE.pyx'
                'libCCE.c'
                ,'cCCE/CCE.c'
                ,'cCCE/entropy.c'
                ,'cCCE/linkedList.c'
                ,'cCCE/llQueue.c'
                ,'cCCE/ntree.c'
                ,'cCCE/arenaAlloc.c'
            ]
            ,include_dirs=['cCCE']
            ,libraries=['m']
            # un-comment for address sanatizer, requires LD_PRELOAD of
            # the libasan shared lib
            ,extra_compile_args=[
                # Default arena size is 1 GB
                '-DARENA_SIZE=1073741824'
                #,'-DDEBUG_ALLOC'
                #,'-DDEBUG_CCE'
                #,'-pg'
                #,'-fsanitize=address'
                #,'-fno-omit-frame-pointer'
            ]
        )
    ]

    # from Cython.Build import cythonize
    # ext_modules = cythonize(ext_modules)

    setup(
        ext_modules = ext_modules
    )

if __name__ == '__main__':
    build()
