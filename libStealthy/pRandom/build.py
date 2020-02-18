#
# U.S. Naval Research Laboratory
# Center for High Assurance Computer Systems
#
def main():
    '''
    This is called by buildExt.py to build the CCE submodule.
    '''
    from distutils.core import setup
    from distutils.extension import Extension

    ext_modules = [
        Extension(
            "PRNG"
            ,sources = [
                # 'PRNG.pyx'
                'PRNG.c'
                ,'SFMT-src-1.5.1/SFMT.c'
                ,'pcg-c-basic-0.9/pcg_basic.c'
                ,'xoshiro/xoshiro256.c'
            ]
            ,include_dirs=[
                'SFMT-src-1.5.1'
                ,'pcg-c-basic-0.9'
                ,'xoshiro'
            ]
            ,extra_compile_args=[
                '-O3'
                ,'-finline-functions'
                ,'-fomit-frame-pointer'
                ,'-DNDEBUG'
                ,'-fno-strict-aliasing'
                ,'-Wmissing-prototypes'
                ,'-Wall'
                ,'-DSFMT_MEXP=216091'
            ]
        )
    ]

    # from Cython.Build import cythonize
    # ext_modules = cythonize(ext_modules)

    setup(
        ext_modules = ext_modules
    )
if __name__ == '__main__':
    main()