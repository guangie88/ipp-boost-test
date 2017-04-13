#include "ipp.h"

#include <cstdio>
#include <cstdlib>

int main()
{
    // Set the size
    const int N=128;

    // Spec and working buffers
    IppsDFTSpec_C_32fc *pDFTSpec=0;
    Ipp8u  *pDFTInitBuf, *pDFTWorkBuf;

    // Allocate complex buffers
    Ipp32fc *pSrc=ippsMalloc_32fc(N);
    Ipp32fc *pDst=ippsMalloc_32fc(N); 

    // Query to get buffer sizes
    int sizeDFTSpec,sizeDFTInitBuf,sizeDFTWorkBuf;
    ippsDFTGetSize_C_32fc(N, IPP_FFT_NODIV_BY_ANY, 
        ippAlgHintAccurate, &sizeDFTSpec, &sizeDFTInitBuf, &sizeDFTWorkBuf);

    // Alloc DFT buffers
    pDFTSpec    = (IppsDFTSpec_C_32fc*)ippsMalloc_8u(sizeDFTSpec);
    pDFTInitBuf = ippsMalloc_8u(sizeDFTInitBuf);
    pDFTWorkBuf = ippsMalloc_8u(sizeDFTWorkBuf);

    // Initialize DFT
    ippsDFTInit_C_32fc(N, IPP_FFT_NODIV_BY_ANY, 
        ippAlgHintAccurate, pDFTSpec, pDFTInitBuf);
    if (pDFTInitBuf) ippFree(pDFTInitBuf);

    // Do the DFT
    ippsDFTFwd_CToC_32fc(pSrc,pDst,pDFTSpec,pDFTWorkBuf);

    //check results
    ippsDFTInv_CToC_32fc(pDst,pDst,pDFTSpec,pDFTWorkBuf);
    int OK=1;
    for (int i=0;i<N;i++) {
        pDst[i].re/=(Ipp32f)N;
        pDst[i].im/=(Ipp32f)N;
        if ((abs(pSrc[i].re-pDst[i].re)>.001) || 
            (abs(pSrc[i].im-pDst[i].im)>.001) ) 
        {
            OK=0;break;
        }
    }

    printf(OK==1?"DFT OK\n":"DFT Fail\n");

    if (pDFTWorkBuf) ippFree(pDFTWorkBuf);
    if (pDFTSpec) ippFree(pDFTSpec);

    ippFree(pSrc);
    ippFree(pDst);
}
