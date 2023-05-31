#include "sentinel1.h"
#include "ipp.h"
#include <Dense>
#include <QR>

Sentinel::Sentinel()
{
    sentinel1PacketDecode.readRawPacket("C:/S1A_S3_RAW__0SDH_20220710T213600_20220710T213625_044043_0541DB_56CE/S1A_S3_RAW__0SDH_20220710T213600_20220710T213625_044043_0541DB_56CE.SAFE/s1a-s3-raw-s-hh-20220710t213600-20220710t213625-044043-0541db.dat");
    sentinel1PacketDecode.getAuxData();
    calcParams();

    // Range direction
    IppStatus st = ippsDFTGetSize_C_32fc(fftLengthRange, IPP_DIV_FWD_BY_N, ippAlgHintAccurate, &sizeRangeFFTSpec, &sizeRangeFFTInitBuf, &sizeRangeFFTWorkBuf);

    pRangeSpec = (IppsDFTSpec_C_32fc*)ippsMalloc_8u(sizeRangeFFTSpec);
    pRangeFFTInitBuf = ippsMalloc_8u(sizeRangeFFTInitBuf);
    pRangeFFTWorkBuf = ippsMalloc_8u(sizeRangeFFTWorkBuf);

    st = ippsDFTInit_C_32fc(fftLengthRange, IPP_DIV_FWD_BY_N, ippAlgHintAccurate, pRangeSpec, pRangeFFTInitBuf);

    // Azimuth direction
    ippsDFTGetSize_C_32fc(fftLengthAzimuth, IPP_DIV_FWD_BY_N, ippAlgHintAccurate, &sizeAzimuthFFTSpec, &sizeAzimuthFFTInitBuf, &sizeAzimuthFFTWorkBuf);

    pAzimuthSpec = (IppsDFTSpec_C_32fc*)ippsMalloc_8u(sizeAzimuthFFTSpec);
    pAzimuthFFTInitBuf = ippsMalloc_8u(sizeAzimuthFFTInitBuf);
    pAzimuthFFTWorkBuf = ippsMalloc_8u(sizeAzimuthFFTWorkBuf);

    ippsDFTInit_C_32fc(fftLengthAzimuth, IPP_DIV_FWD_BY_N, ippAlgHintAccurate, pAzimuthSpec, pAzimuthFFTInitBuf);

    std::vector<double> coeffVelocity;
    std::vector<double> coeffPosition;

    interp(sentinel1PacketDecode.time, sentinel1PacketDecode.normOfVelocity, sentinel1PacketDecode.normOfPosition, coeffVelocity, coeffPosition);

    for (int p = 0; p < sentinel1PacketDecode.header.size(); p++)
    {
        double time = sentinel1PacketDecode.header.at(p).FineTime + static_cast<double>(sentinel1PacketDecode.header.at(p).Time);
        timeS.push_back(time);

        double vfitted = coeffVelocity[0] +
            coeffVelocity[1] * time +
            coeffVelocity[2] * (pow(time, 2)) +
            coeffVelocity[3] * (pow(time, 3)) +
            coeffVelocity[4] * (pow(time, 4)) +
            coeffVelocity[5] * (pow(time, 5)) +
            coeffVelocity[6] * (pow(time, 6)) +
            coeffVelocity[7] * (pow(time, 7)) +
            coeffVelocity[8] * (pow(time, 8)) +
            coeffVelocity[9] * (pow(time, 9)) +
            coeffVelocity[10] * (pow(time, 10));
        velocity.push_back(vfitted);

        vfitted = coeffPosition[0] +
            coeffPosition[1] * time +
            coeffPosition[2] * (pow(time, 2)) +
            coeffPosition[3] * (pow(time, 3)) +
            coeffPosition[4] * (pow(time, 4)) +
            coeffPosition[5] * (pow(time, 5)) +
            coeffPosition[6] * (pow(time, 6)) +
            coeffPosition[7] * (pow(time, 7)) +
            coeffPosition[8] * (pow(time, 8)) +
            coeffPosition[9] * (pow(time, 9)) +
            coeffPosition[10] * (pow(time, 10));
        position.push_back(vfitted);
    }

    //getEffectiveVelocity(sentinel1PacketDecode);
    
    //getRangeFilter();
    //
    //dopplerCentroidEstimation();
    //
    //for (size_t i = 0; i < sentinel1PacketDecode.out.size(); i++) {
    //    ippsDFTFwd_CToC_32fc(sentinel1PacketDecode.out[i].data(), sentinel1PacketDecode.out[i].data(), pRangeSpec, pRangeFFTWorkBuf);
    //    ippsMul_32fc(sentinel1PacketDecode.out[i].data(), refFunc.data(), sentinel1PacketDecode.out[i].data(), fftLengthRange);
    //}
    //
    //for (size_t j = 0; j < fftLengthRange; j++) {
    //    std::vector<Ipp32fc> az;
    //    for (size_t i = 0; i < fftLengthAzimuth; i++) {
    //        az.push_back(sentinel1PacketDecode.out[i][j]);
    //    }
    //    ippsDFTFwd_CToC_32fc(az.data(), az.data(), pAzimuthSpec, pAzimuthFFTWorkBuf);
    //    for (size_t i = 0; i < fftLengthAzimuth; i++) {
    //        sentinel1PacketDecode.out[i][j] = az[i];
    //    }
    //}
}

Sentinel::~Sentinel()
{

}

void Sentinel::dopplerCentroidEstimation()
{
    for (size_t j = 0; j < fftLengthRange; j++) {
        
        std::vector<Ipp32fc> az;
        
        for (size_t i = 0; i < fftLengthAzimuth; i++) {
            az.push_back(sentinel1PacketDecode.out[i][j]);
        }

        std::complex<double> res(0.0, 0.0);

        for (size_t i = 0; i < fftLengthAzimuth - 1; i++) {
            std::complex<double> oneTmp(az.at(i).re, az.at(i).im);
            std::complex<double> twoTmp(az.at(i + 1).re, -az.at(i + 1).im);

            res += oneTmp * twoTmp;
        }

        dopplerCentroid.push_back( - std::arg(res) / (2 * M_PI * sentinel1PacketDecode.header.at(0).PulseRepetitionInterval));
    }
    
}

void Sentinel::polyfit(const std::vector<double>& t, const std::vector<double>& v, std::vector<double>& coeff, int order) {
    // Create Matrix Placeholder of size n x k, n= number of datapoints, k = order of polynomial, for exame k = 3 for cubic polynomial
    Eigen::MatrixXd T(t.size(), order + 1);
    Eigen::VectorXd V = Eigen::VectorXd::Map(&v.front(), v.size());
    Eigen::VectorXd result;

    // check to make sure inputs are correct
    //assert(t.size() == v.size());
    //assert(t.size() >= order + 1);
    // Populate the matrix
    for (size_t i = 0; i < t.size(); ++i)
    {
        for (size_t j = 0; j < order + 1; ++j)
        {
            T(i, j) = pow(t.at(i), j);
        }
    }
    std::cout << T << std::endl;

    // Solve for linear least square fit
    result = T.householderQr().solve(V);
    coeff.resize(order + 1);
    for (int k = 0; k < order + 1; k++)
    {
        coeff[k] = result[k];
    }

}

void Sentinel::interp(const std::vector<double>& time, const std::vector<double>& velocity, const std::vector<double>& position,
                      std::vector<double>& coeffVelocity, std::vector<double>& coeffPosition) {
    polyfit(time, velocity, coeffVelocity, 10);
    polyfit(time, position, coeffPosition, 10);
}

void Sentinel::MeanOfRawData(ComplexMatrix& rawData, RawDataAnalysis& rawDataAnalysis) {
    double sumI = 0.0;
    double sumQ = 0.0;

    for (size_t i = 0; i < rawData.re.size(); i++) {
        sumI += static_cast<double>(std::accumulate(rawData.re[i].begin(), rawData.re[i].end(), 0));
    }

    for (size_t i = 0; i < rawData.im.size(); i++) {
        sumQ += static_cast<double>(std::accumulate(rawData.im[i].begin(), rawData.im[i].end(), 0));
    }

    rawDataAnalysis.meanOfRawDataI = sumI;
    rawDataAnalysis.meanOfRawDataQ = sumQ;

    sumI = sumI / (static_cast<double>(rawData.rows) * static_cast<double>(rawData.cols));
    sumQ = sumQ / (static_cast<double>(rawData.rows) * static_cast<double>(rawData.cols));
}



float Sentinel::getNormSateliteVelocity(Sentinel1PacketDecode& sentinel1PacketDecode, uint64_t numberOfNavigPacket) {
    return std::sqrt(std::pow(sentinel1PacketDecode.positionVelocityTime.at(numberOfNavigPacket).XvelocityECEF, 2) +
        std::pow(sentinel1PacketDecode.positionVelocityTime.at(numberOfNavigPacket).YvelocityECEF, 2) +
        std::pow(sentinel1PacketDecode.positionVelocityTime.at(numberOfNavigPacket).ZvelocityECEF, 2));
}

/*

Vr = sqrt(Vg * Vs)

Vg - ground velocity, which is a positive scalar representing the antenna
beam footprint velocity when projected onto the ground.

Vs - satellite velocity, which is the norm of the satellite velocity vector
expressed in ECR coordinates.

Note that, as Vg is range dependent, Ve will also be range dependent. 
The Range-Doppler azimuth focusing can fully take into account this variation.

The effective radar velocity is also azimuth dependent, as it depends on the satellite
position in orbit and its height above the Earth. This dependency will be taken into
account in IPF by updating V for each azimuth block.
*/

float Sentinel::getEffectiveVelocity(Sentinel1PacketDecode& sentinel1PacketDecode) {
    for (size_t i = 0; i < fftLengthAzimuth; i++) {
        float V = velocity.at(i);
        float H = position.at(i);

        float groundVelocity = 0.0;

        double a = 6378137.0; // WGS84 semi major axis
        double b = 6356752.3142; // WGS84 semi minor axis

        double lat = std::atan(sentinel1PacketDecode.positionVelocityTime.at(0).ZAxisPositionECEF /
            sentinel1PacketDecode.positionVelocityTime.at(0).XAxisPositionECEF);

        double local_earth_rad =
            std::sqrt(
                (std::pow(a * a * cos(lat), 2) + std::pow(b * b * sin(lat), 2)) /
                (std::pow(b * sin(lat), 2) + std::pow(a * cos(lat), 2)));

        double cos_beta = 0.0;
        double this_ground_velocity = 0.0;
        double effective_velocity = 0.0;

        std::vector<float> effectiveVelocityTemp;

        for (size_t j = 0; j < fftLengthRange; j++) {
            cos_beta = (local_earth_rad * local_earth_rad + H * H - slantRange[j] * slantRange[j]) / (2 * local_earth_rad * H);
            this_ground_velocity = local_earth_rad * V * cos_beta / H;
            effective_velocity = std::sqrt(V * this_ground_velocity);
            effectiveVelocityTemp.push_back(effective_velocity);
        }

        effectiveVelocity.emplace_back(effectiveVelocityTemp);
    }
    
    return 0.0;
}

float Sentinel::getNormSatelitePosition(Sentinel1PacketDecode& sentinel1PacketDecode, uint64_t numberOfNavigPacket) {
    return std::sqrt(std::pow(sentinel1PacketDecode.positionVelocityTime.at(numberOfNavigPacket).XAxisPositionECEF, 2) +
        std::pow(sentinel1PacketDecode.positionVelocityTime.at(numberOfNavigPacket).YAxisPositionECEF, 2) +
        std::pow(sentinel1PacketDecode.positionVelocityTime.at(numberOfNavigPacket).ZAxisPositionECEF, 2));
}

int Sentinel::rangeCompression(ComplexMatrix& rawData, bool calcIfft) {
    Ipp32fc* correl = nullptr;

    int fft_len = 0;

    correl = ippsMalloc_32fc(fft_len);
    if (correl == nullptr)
        return -1;

    ippsZero_32fc(correl, fft_len);

    //for (size_t i = 0; i < rawData.rows; i++)
    //{
    //    //ippsFFTFwd_CToC_32fc(, correl, pSpec, pDFTWorkBuf);
    //
    //    ippsMul_32fc(correl, refFunc.data(), correl, fft_len);
    //    if (calcIfft) {
    //        ippsFFTInv_CToC_32fc(correl, correl, pSpec, pDFTWorkBuf);
    //    }
    //}
    //
    //ippsFree(pDFTSpec);
    //ippsFree(pDFTInitBuf);
    //ippsFree(pDFTWorkBuf);
    //ippsFree(correl);

    return 0;
}

void Sentinel::StandardDeviationsOfRawData(ComplexMatrix& rawData, RawDataAnalysis& rawDataAnalysis) {
    double k = 1 / static_cast<double>(rawData.rows) * static_cast<double>(rawData.cols);

    double sumIMinusMean = 0.0;
    double sumQMinusMean = 0.0;

    for (size_t i = 0; i < rawData.rows; i++) {
        for (size_t j = 0; j < rawData.cols; j++) {
            sumIMinusMean += std::pow(rawData.re[i][j] - rawDataAnalysis.meanOfRawDataI, 2);
            sumQMinusMean += std::pow(rawData.im[i][j] - rawDataAnalysis.meanOfRawDataQ, 2);
        }
    }

    rawDataAnalysis.standardDeviationsOfRawDataI = std::sqrt(k * sumIMinusMean);
    rawDataAnalysis.standardDeviationsOfRawDataQ = std::sqrt(k * sumQMinusMean);
}

void Sentinel::IQGainImbalance(ComplexMatrix& rawData, RawDataAnalysis& rawDataAnalysis) {
    double k = 1 / static_cast<double>(rawData.rows) * static_cast<double>(rawData.cols);

    rawDataAnalysis.IQGainImbalance = rawDataAnalysis.standardDeviationsOfRawDataI / rawDataAnalysis.standardDeviationsOfRawDataQ;
    rawDataAnalysis.IQGainImbalanceLowerBounds = 1 - 3 / std::sqrt(k);
    rawDataAnalysis.IQGainImbalanceUpperBounds = 1 + 3 / std::sqrt(k);
}
double Sentinel::getRowSumI(ComplexMatrix& rawData, uint64_t number) {
    return std::accumulate(rawData.re[number].begin(), rawData.re[number].end(), 0);
}

double Sentinel::getRow2SumI(ComplexMatrix& rawData, uint64_t number) {
    double sum = 0.0;
    for (size_t i = 0; i < rawData.cols; i++) {
        sum += std::pow(rawData.re[number].at(i), 2);
    }
    return sum;
}

double Sentinel::getRowSumQ(ComplexMatrix& rawData, uint64_t number) {
    return std::accumulate(rawData.im[number].begin(), rawData.im[number].end(), 0);
}

double Sentinel::getRow2SumQ(ComplexMatrix& rawData, uint64_t number) {
    double sum = 0.0;
    for (size_t i = 0; i < rawData.cols; i++) {
        sum += std::pow(rawData.im[number].at(i), 2);
    }
    return sum;
}

double Sentinel::getRowSumMultIQ(ComplexMatrix& rawData, uint64_t number) {
    double sum = 0.0;
    for (size_t i = 0; i < rawData.cols; i++) {
        sum += rawData.re[number].at(i) * rawData.im[number].at(i);
    }
    return sum;
}

void Sentinel::IQQuadratureDeparture(ComplexMatrix& rawData, RawDataAnalysis& rawDataAnalysis) {
    std::vector<double> C;

    for (size_t i = 0; i < rawData.rows; i++) {
        C.push_back((getRowSumMultIQ(rawData, i) - getRow2SumI(rawData, i) * getRow2SumI(rawData, i) / rawData.cols) /
            std::sqrt((getRow2SumI(rawData, i) - getRow2SumI(rawData, i) * getRow2SumI(rawData, i) / rawData.cols) *
                (getRow2SumQ(rawData, i) - getRow2SumQ(rawData, i) * getRow2SumQ(rawData, i) / rawData.cols)));
    }

    for (size_t i = 0; i < rawData.rows; i++)
    {
        C[i] = 0.5 * std::log((1 + C[i]) / (1 - C[i]));
    }

    double sum = std::accumulate(C.begin(), C.end(), 0) / C.size();

    rawDataAnalysis.IQQuadratureDeparture = std::asin(sum);
    rawDataAnalysis.IQQuadratureDepartureLowerBounds = 0;
    rawDataAnalysis.IQQuadratureDepartureUpperBounds = 0;
}

void Sentinel::SetStatisticsSignificanceFlags(ComplexMatrix& rawData, RawDataAnalysis& rawDataAnalysis) {
    double k = 1 / static_cast<double>(rawData.rows) * static_cast<double>(rawData.cols);

    rawDataAnalysis.IBiasSignificanceFlag = 0.0;
    rawDataAnalysis.IQGainSignificanceFlag = 0.0;
    rawDataAnalysis.QBiasSignificanceFlag = 0.0;
}

void Sentinel::CorrectBiases(ComplexMatrix& rawData, RawDataAnalysis& rawDataAnalysis) {
    for (size_t i = 0; i < rawData.rows; i++) {
        std::transform(rawData.re[i].begin(), rawData.re[i].end(), rawData.re[i].begin(),
            [rawDataAnalysis](double d) -> double { return d + rawDataAnalysis.meanOfRawDataI; });
    }

    for (size_t i = 0; i < rawData.rows; i++) {
        std::transform(rawData.im[i].begin(), rawData.im[i].end(), rawData.re[i].begin(),
            [rawDataAnalysis](double d) -> double { return d + rawDataAnalysis.meanOfRawDataQ; });
    }
}

void Sentinel::CorrectGainImbalance(ComplexMatrix& rawData, RawDataAnalysis& rawDataAnalysis) {
    for (size_t i = 0; i < rawData.rows; i++) {
        std::transform(rawData.im[i].begin(), rawData.im[i].end(), rawData.re[i].begin(),
            [rawDataAnalysis](double d) -> double { return d * rawDataAnalysis.IQGainImbalance; });
    }
}

void Sentinel::CorrectQChannelForNonOrthogonality(ComplexMatrix& rawData, RawDataAnalysis& rawDataAnalysis) {
    for (size_t i = 0; i < rawData.rows; i++) {
        for (size_t j = 0; j < rawData.cols; j++)
        {
            rawData.im[i][j] = rawData.im[i][j] / cos(rawDataAnalysis.IQQuadratureDeparture) - rawData.re[i][j] * std::tan(rawDataAnalysis.IQQuadratureDeparture);
        }
    }
}

bool Sentinel::checkRefFunc() {
    double TXPRR = sentinel1PacketDecode.header[0].TxRampRate;
    double TXPSF = sentinel1PacketDecode.header[0].TxPulseStartFreq;
    double TXPL = sentinel1PacketDecode.header[0].TxPulseLength;

    bool f = false;

    for (size_t i = 0; i < sentinel1PacketDecode.out.size(); i++) {
        if (TXPRR != sentinel1PacketDecode.header[i].TxRampRate) {
            f = true;
        }
        if (TXPSF != sentinel1PacketDecode.header[i].TxPulseStartFreq) {
            f = true;
        }
        if (TXPL != sentinel1PacketDecode.header[i].TxPulseLength) {
            f = true;
        }
    }

    return f;
}

template<typename T>
std::complex<double> Sentinel::NominalImageReplicaGeneration(T  TXPL, T TXPSF, T TXPRR, double i) {
    return std::exp(std::complex<double>(2 * M_PI * (((TXPSF - TXPRR * (-TXPL / 2)) * i + (TXPRR / 2) * std::pow(i, 2)))));
}

double Sentinel::MigrationFactor(double carrierFrequency, double freqAzimuth, double effectiveRadarVelocity) {
    return std::sqrt(1 - (std::pow(speedOfLight / carrierFrequency, 2) * std::pow(freqAzimuth, 2) / 4 * std::pow(effectiveRadarVelocity, 2)));
}

template<typename T>
T Sentinel::AzimuthFMRate(T carrierFrequency, T frequencyDopplerCentroid, T effectiveRadarVelocity, T slantRange) {
    return 2 * std::pow(effectiveRadarVelocity, 2) * MigrationFactor<T>(carrierFrequency, frequencyDopplerCentroid, effectiveRadarVelocity) / ((speedOfLight / carrierFrequency) * slantRange);
}

int Sentinel::getRangeFilter() {
    double TXPRR = 1e12 * sentinel1PacketDecode.header[0].TxRampRate;
    double TXPSF = 1e6 * sentinel1PacketDecode.header[0].TxPulseStartFreq;
    double TXPL = 1e-6 * sentinel1PacketDecode.header[0].TxPulseLength;

    Ipp32fc tmp;

    for (size_t i = 0; i < static_cast<int>(TXPL * rangeSampleFreq); i++) {
        double samplingT = (double)i / (double)rangeSampleFreq - (double)TXPL / 2.0;
        double sum1 = (TXPSF + TXPRR * (TXPL / 2)) * samplingT;
        double sum2 = TXPRR * std::pow(samplingT, 2) / 2;

        double arg = 2 * M_PI * (sum1 + sum2);
        double re = cos(arg);
        double im = sin(arg);

        tmp.re = static_cast<float>(re);
        tmp.im = static_cast<float>(im);

        refFunc.emplace_back(tmp);
    }

    for (size_t i = static_cast<int>(TXPL * rangeSampleFreq); i < fftLengthRange; i++) {
        tmp.re = 0.0;
        tmp.im = 0.0;
        refFunc.push_back(tmp);
    }

    ippsDFTFwd_CToC_32fc(refFunc.data(), refFunc.data(), pRangeSpec, pRangeFFTWorkBuf);

    ippsConj_32fc(refFunc.data(), refFunc.data(), fftLengthRange);
    return 0;
}

int Sentinel::calcParams() {
    fftLengthRange = sentinel1PacketDecode.out[0].size();

    fftLengthAzimuth = sentinel1PacketDecode.out.size();

    rangeSampleFreq = RangeDecimationToSampleRate(sentinel1PacketDecode.header[0].RangeDecimation);
    rangeSamplePeriod = 1.0 / rangeSampleFreq;
    rangeStartTime = sentinel1PacketDecode.header[0].SamplingWindowStartTime + 320.0 * 1e-6 / (8 * 37.53472224);
    azSampleFreq = 1 / sentinel1PacketDecode.header[0].PulseRepetitionInterval;
    azSamplePeriod = sentinel1PacketDecode.header[0].PulseRepetitionInterval;

    for (size_t i = 0; i < sentinel1PacketDecode.out[0].size(); i++) {
        fastTime.push_back(rangeStartTime + i * rangeSamplePeriod);
    }

    for (size_t i = 0; i < sentinel1PacketDecode.out[0].size(); i++) {
        slantRange.push_back((sentinel1PacketDecode.header[0].Rank * sentinel1PacketDecode.header[0].PulseRepetitionInterval + fastTime.at(i)) * speedOfLight / 2);
    }
    float SWL = sentinel1PacketDecode.header[0].NumberOfQuads * 2 / rangeSampleFreq;

    for (size_t i = 0; i < sentinel1PacketDecode.out.size(); i++) {
        azimuthFreq.push_back((-azSampleFreq / 2) + i / (sentinel1PacketDecode.header[0].PulseRepetitionInterval * sentinel1PacketDecode.out.size()));
    }
    return 0;
}