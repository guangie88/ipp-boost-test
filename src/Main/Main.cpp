#include "ipp.h"

#include "boost/filesystem.hpp"
#include "gtest/gtest.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>

// boost
using boost::filesystem::exists;

int dft()
{
    // Set the size
    static constexpr auto FILE_PATH = "data/data.bin";
    static constexpr int N = 128;

    if (!exists(FILE_PATH))
    {
        printf("'%s' does not exist!\n", FILE_PATH);
        return 1;
    }

    auto filePtr = fopen(FILE_PATH, "rb");

    // Allocate complex buffers
    Ipp32fc *pSrc = ippsMalloc_32fc(N);
    Ipp32fc *pDst = ippsMalloc_32fc(N); 

    if (filePtr)
    {
        fclose(filePtr);
    }

    // Spec and working buffers
    IppsDFTSpec_C_32fc *pDFTSpec = 0;
    Ipp8u *pDFTInitBuf, *pDFTWorkBuf;

    // Query to get buffer sizes
    int sizeDFTSpec, sizeDFTInitBuf, sizeDFTWorkBuf;

    ippsDFTGetSize_C_32fc(N, IPP_FFT_NODIV_BY_ANY, 
        ippAlgHintAccurate, &sizeDFTSpec, &sizeDFTInitBuf, &sizeDFTWorkBuf);

    // Alloc DFT buffers
    pDFTSpec = (IppsDFTSpec_C_32fc*)ippsMalloc_8u(sizeDFTSpec);
    pDFTInitBuf = ippsMalloc_8u(sizeDFTInitBuf);
    pDFTWorkBuf = ippsMalloc_8u(sizeDFTWorkBuf);

    // Initialize DFT
    ippsDFTInit_C_32fc(N, IPP_FFT_NODIV_BY_ANY, 
        ippAlgHintAccurate, pDFTSpec, pDFTInitBuf);

    if (pDFTInitBuf) ippFree(pDFTInitBuf);

    // Do the DFT
    ippsDFTFwd_CToC_32fc(pSrc, pDst, pDFTSpec, pDFTWorkBuf);

    // check results
    ippsDFTInv_CToC_32fc(pDst, pDst, pDFTSpec, pDFTWorkBuf);
    int OK = 1;

    for (int i = 0; i < N; i++)
    {
        pDst[i].re /= (Ipp32f)N;
        pDst[i].im /= (Ipp32f)N;

        if ((::fabs(pSrc[i].re - pDst[i].re) > .001) || 
            (::fabs(pSrc[i].im - pDst[i].im) > .001)) 
        {
            OK=0;
            break;
        }
    }

    printf(OK == 1 ? "DFT OK\n" : "DFT Fail\n");

    if (pDFTWorkBuf) ippFree(pDFTWorkBuf);
    if (pDFTSpec) ippFree(pDFTSpec);

    ippFree(pSrc);
    ippFree(pDst);

    return 0;
}

TEST(dft_test, dft_test_description)
{
    EXPECT_EQ(0, dft());
}

int main(int argc, char * argv[])
{
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
