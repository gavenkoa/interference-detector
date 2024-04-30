
# Configuring VS Code

On Windows we need to set the prefix where rtlsdr header and lib can be found. Update your workspace env vars:

    cat .vscode/settings.json
    {
        "cmake.environment": {
            "CONDA_LIB_ROOT": "C:\\opt\\radioconda\\Library"
        },
    }

https://code.visualstudio.com/docs/editor/variables-reference
:   Predefined variables.

https://github.com/microsoft/vscode-cmake-tools/blob/main/docs/how-to.md
:   Set up include paths for C++ IntelliSense

https://github.com/microsoft/vscode-cmake-tools/blob/main/docs/configure.md
:   The CMake Tools configure step

https://github.com/microsoft/vscode-cmake-tools/blob/main/docs/cmake-settings.md
:   About `cmake.environment`.

To run from Mingw64 environment launch VS Code from Mingw64 configured console.

https://code.visualstudio.com/docs/cpp/config-mingw
:   Using GCC with MinGW

# 

https://wiki.gnuradio.org/index.php/Sample_Rate_Change
:   Sample Rate Change. About "Interpolating FIR Filter" +  "Low-Pass Filter Taps" block to increase sample rate for wider bandwidth.


# libusb

https://libusb.sourceforge.io/api-1.0/group__libusb__syncio.html
https://github.com/libusb/libusb
https://sourceforge.net/projects/libusb/
https://github.com/libusb/libusb/wiki/FAQ
https://github.com/libusb/libusb/wiki

# librtlsdr

There are several competing:

* https://github.com/librtlsdr/librtlsdr/issues
* https://github.com/steve-m/librtlsdr/issues
* 

https://pyrtlsdr.readthedocs.io/en/latest/rtlsdr.html
:   Python binding.

# FFT

https://www.fftw.org/
:    FFTW is a C subroutine library for computing the discrete Fourier transform (DFT) in one or more dimensions.


