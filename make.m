function make()
    if (~verLessThan('matlab','9.4'))
        mex -R2018a -output as_zsl as_zsl.cpp as_zsl_matlab.cpp
    else
        mex -largeArrayDims -output as_zsl as_zsl.cpp as_zsl_matlab.cpp
    end
end