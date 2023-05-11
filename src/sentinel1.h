#pragma once
#include "sentinel1_packet_decode.h"
#include <complex>
#include <algorithm>
#include <numeric>
#include <vector>

struct ComplexMatrix {
    std::vector<std::vector<double>> re;
    std::vector<std::vector<double>> im;
    uint64_t rows = 0;
    uint64_t cols = 0;
};

class Sentinel
{
public:
    Sentinel::Sentinel();

    Sentinel::~Sentinel();

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
    void NominalImageReplicaGeneration(T  TXPL, T TXPSF, T TXPRR);

    template<typename T>
    T getRangeReferenceFunction(T  TXPL, T TXPSF, T TXPRR);

    template<typename T>
    T MigrationFactor(T carrierFrequency, T freqAzimuth, T effectiveRadarVelocity);

    template<typename T>
    T AzimuthFMRate(T carrierFrequency, T frequencyDopplerCentroid, T effectiveRadarVelocity, T slantRange);

    int rangeCompression(ComplexMatrix& rawData, bool calcIfft = true);
private:
    const double speedOfLight = 299792458.0;
    double rangeStartTime = 0.0;
    double suppressedDataTime = 0.0;
    double wavelength = speedOfLight / 5.405e9;
    double rangeSampleFreq = 0.0;
    double rangeSamplePeriod = 0.0;
    double azSampleFreq = 0.0;
    double azSamplePeriod = 0.0;
};



