#include <vector>
#include <iostream>

// Calculate the mean of the raw data.
template<typename T>
double MeanOfRawData(std::vector<std::vector<T>> rawData) {
    double sum = 0.0;

    for (size_t i = 0; i < rawData.size(); i++) {
        sum += static_cast<double>(accumulate(rawData[i].begin(), rawData[i].end(), 0));
    }

    sum = sum / (static_cast<double>(rawData.size()) * static_cast<double>(rawData.at(0).size()));
    return sum;
}

// Calculate the standard deviations of the raw data.
template<typename T>
T StandardDeviationsOfRawData(std::vector<std::vector<T>> rawData) {
    return static_cast <typeid(T).name()>(0);
}

// Calculate the IQ gain imbalance
template<typename T>
T IQGainImbalance(std::vector<std::vector<T>> rawData) {
    return static_cast <typeid(T).name()>(0);
}

// IQ quadrature departure
template<typename T>
T IQQuadratureDeparture(std::vector<std::vector<T>> rawData) {
    return static_cast <typeid(T).name()>(0);
}

// Set the statistics significance flags


int main(int argc, char *argv[])
{
    return 0;
}
