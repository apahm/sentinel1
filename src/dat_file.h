#pragma once
#include <iostream>
#include <string>

struct imop_desc_rec {

};

struct sdr_data_rec {

};

struct pdr_data_rec
{

};

struct sar_dat_file {
	imop_desc_rec _imop_desc_rec;
	sdr_data_rec _sdr_data_rec;
	pdr_data_rec _pdr_data_rec;
};
