# wdenoise
Wavelet Denoising in ANSI C using empirical bayes thresholding and a host of other thresholding methods.

|**[wdenoise object](https://github.com/rafat/wdenoise/wiki/wdenoise)**| WDenoise Object, Parameters and Functions|
|:--------------------------------------------------------|:-------------------------------------|
|**[Wdenoise 1D Example using EBayesThresh](https://github.com/rafat/wdenoise/wiki/Example-Code-1-:-wdenoise-(EBayesThresh))**| Example Code 1 : wdenoise (EBayesThresh)|
|**[Wdenoise 1D Example using fdrshrink](https://github.com/rafat/wdenoise/wiki/Example-Code-2-:-wdenoise)**| Example Code 2 : wdenoise|
|**[Image Denoising](https://github.com/rafat/wdenoise/wiki/Example-Code-3-:-Image-Denoising)**| Example Code 3 : Image Denoising using EBayesThresh and Visushrink|

_[Full Documentation is available here](https://github.com/rafat/wdenoise/wiki)_

## Dependencies

Git and CMake

## Getting Started
```
Clone the project.
cd to directory
cmake .
make
```
## Credits

The EbayesThresh software package was originally developed by Bernard W. Silverman and Ludger Evers, with extensions introduced by Kan Xu, Peter Carbonetto and Matthew Stephens in the Department of Statistics at the University of Chicago.

_[R code is available here](https://github.com/stephenslab/EbayesThresh)_

Matlab Version of the code is written by A. ANTONIADIS, M. JENSEN, I. JOHNSTONE & B. W. SILVERMAN.

_[Matlab code is available here](http://www-ljk.imag.fr/membres/Anestis.Antoniadis/EBayesThresh/)_

## License

All the codes available in the src folder of this repository are licensed under GNU General Public License 3.0

_[Full Text of the License is available here](https://github.com/rafat/wdenoise/blob/master/LICENSE)_
