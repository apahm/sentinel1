#pragma once
#include "sentinel1_packet_decode.h"
#include <complex>
#include <algorithm>
#include <numeric>
#include <vector>
#include "ipp.h"

struct ComplexMatrix {
    std::vector<std::vector<double>> re;
    std::vector<std::vector<double>> im;
    uint64_t rows = 0;
    uint64_t cols = 0;
};

class Sentinel
{
public:
    Sentinel();

    ~Sentinel();
    void polyfit(const std::vector<double>& t, const std::vector<double>& v, std::vector<double>& coeff, int order);
    void interp(const std::vector<double>& time,
        const std::vector<double>& velocity,
        const std::vector<double>& position,
        std::vector<double>& coeffVelocity,
        std::vector<double>& coeffPosition);

    float getNormSateliteVelocity(Sentinel1PacketDecode& sentinel1PacketDecode, uint64_t numberOfNavigPacket);
    float getNormGroundVelocity(Sentinel1PacketDecode& sentinel1PacketDecode);
    float getNormSatelitePosition(Sentinel1PacketDecode& sentinel1PacketDecode, uint64_t numberOfNavigPacket);
    float getEffectiveVelocity(Sentinel1PacketDecode& sentinel1PacketDecode);

    // Calculate the mean of the raw data.
    void MeanOfRawData(ComplexMatrix& rawData, RawDataAnalysis& rawDataAnalysis);
    
    // Calculate the standard deviations of the raw data.
    void StandardDeviationsOfRawData(ComplexMatrix& rawData, RawDataAnalysis& rawDataAnalysis);

    // Calculate the IQ gain imbalance
    void IQGainImbalance(ComplexMatrix& rawData, RawDataAnalysis& rawDataAnalysis);

    double getRowSumI(ComplexMatrix& rawData, uint64_t number);
    double getRow2SumI(ComplexMatrix& rawData, uint64_t number);
    double getRowSumQ(ComplexMatrix& rawData, uint64_t number);
    double getRow2SumQ(ComplexMatrix& rawData, uint64_t number);
    double getRowSumMultIQ(ComplexMatrix& rawData, uint64_t number);

    // IQ quadrature departure
    void IQQuadratureDeparture(ComplexMatrix& rawData, RawDataAnalysis& rawDataAnalysis);

    // Set the statistics significance flags
    void SetStatisticsSignificanceFlags(ComplexMatrix& rawData, RawDataAnalysis& rawDataAnalysis);

    // Remove constant biases from the I and Q channels. This correction enforces a zero mean to the
    // I and Q components of the signal.
    void CorrectBiases(ComplexMatrix& rawData, RawDataAnalysis& rawDataAnalysis);

    // Correct for gain imbalance between the I and Q channels.This correction
    // imposes the same standard deviation for the I and Q channels.
    void CorrectGainImbalance(ComplexMatrix& rawData, RawDataAnalysis& rawDataAnalysis);

    // Correct the Q channel for non - orthogonality
    void CorrectQChannelForNonOrthogonality(ComplexMatrix & rawData, RawDataAnalysis & rawDataAnalysis);

    template<typename T>
    std::complex<double> NominalImageReplicaGeneration(T  TXPL, T TXPSF, T TXPRR, double i);

    bool checkRefFunc();

    double MigrationFactor(double carrierFrequency, double freqAzimuth, double effectiveRadarVelocity);
    
    template<typename T>
    T AzimuthFMRate(T carrierFrequency, T frequencyDopplerCentroid, T effectiveRadarVelocity, T slantRange);

    int rangeCompression(ComplexMatrix& rawData, bool calcIfft = true);

    float RangeDecimationToSampleRate(RangeDecimation& rangeDecimation) {
        float fRef = 37.53472224 * 1e6;
        float sampleRate = 0.0;
        switch (rangeDecimation) {
        case RangeDecimation::Filter0:
            sampleRate = 3.0 * fRef;
            break;
        case RangeDecimation::Filter1:
            sampleRate = (8.0 / 3.0) * fRef;
            break;
        case RangeDecimation::Filter3:
            sampleRate = (20.0 / 9.0) * fRef;
            break;
        case RangeDecimation::Filter4:
            sampleRate = (16.0 / 9.0) * fRef;
            break;
        case RangeDecimation::Filter5:
            sampleRate = (3.0 / 2.0) * fRef;
            break;
        case RangeDecimation::Filter6:
            sampleRate = (4.0 / 3.0) * fRef;
            break;
        case RangeDecimation::Filter7:
            sampleRate = (2.0 / 3.0) * fRef;
            break;
        case RangeDecimation::Filter8:
            sampleRate = (12.0 / 7.0) * fRef;
            break;
        case RangeDecimation::Filter9:
            sampleRate = (5.0 / 4.0) * fRef;
            break;
        case RangeDecimation::Filter10:
            sampleRate = (6.0 / 13.0) * fRef;
            break;
        case RangeDecimation::Filter11:
            sampleRate = (16.0 / 11.0) * fRef;
            break;
        default:
            break;
        }
        return sampleRate;
    }

    int calcParams();

    int getRangeFilter();

    void dopplerCentroidEstimation();
    
    Sentinel1PacketDecode sentinel1PacketDecode;
    const double speedOfLight = 299792458.0;
    double rangeStartTime = 0.0;
    double wavelength = speedOfLight / 5.405e9;
    double rangeSampleFreq = 0.0;
    double rangeSamplePeriod = 0.0;
    double azSampleFreq = 0.0;
    double azSamplePeriod = 0.0;

    std::vector<double> fastTime;
    std::vector<double> slantRange;
    std::vector<double> azimuthFreq;
    std::vector<std::vector<double>> effectiveVelocity;
    std::vector<double> timeS;
    std::vector<double> velocity;
    std::vector<double> position;
    std::vector<std::complex<double>> dopplerCentroid;

    uint32_t rangeLengthBuffer = 0;
    uint32_t azimuthLengthBuffer = 0;
    std::vector<Ipp32fc> refFunc;

    int sizeRangeFFTSpec = 0;
    int sizeRangeFFTInitBuf = 0;
    int sizeRangeFFTWorkBuf = 0;
    Ipp8u* pRangeFFTSpec = nullptr;
    Ipp8u* pRangeFFTInitBuf = nullptr;
    Ipp8u* pRangeFFTWorkBuf = nullptr;
    IppsDFTSpec_C_32fc* pRangeSpec = nullptr;

    int sizeAzimuthFFTSpec = 0;
    int sizeAzimuthFFTInitBuf = 0;
    int sizeAzimuthFFTWorkBuf = 0;
    Ipp8u* pAzimuthFFTSpec = nullptr;
    Ipp8u* pAzimuthFFTInitBuf = nullptr;
    Ipp8u* pAzimuthFFTWorkBuf = nullptr;
    IppsDFTSpec_C_32fc* pAzimuthSpec = nullptr;

    uint32_t fftLengthRange = 0.0;
    uint32_t fftSizeRange = 0.0;

    uint32_t fftLengthAzimuth = 0.0;
    uint32_t fftSizeAzimuth = 0.0;
};



