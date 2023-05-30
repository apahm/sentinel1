#pragma once
#include <iostream>
#include <string>

struct trl_sar_desc_rec {

};

struct trl_dataset_sum_rec {

};

struct trl_qual_sum_rec {

};

struct trl_sdr_hist_rec {

};

struct trl_pdr8_hist_rec {

};

struct trl_proc_parm_rec {

};

struct trl_att_data_rec {

};

struct trl_radi_data_rec {

};

struct trl_radi_comp_rec {

};

struct sar_trl_file {
	trl_sar_desc_rec _sar_desc_rec;
	trl_dataset_sum_rec _dataset_sum_rec;
	trl_qual_sum_rec _qual_sum_rec;
	trl_sdr_hist_rec _sdr_hist_rec;
	trl_pdr8_hist_rec _pdr8_hist_rec;
	trl_proc_parm_rec _proc_parm_rec;
	trl_att_data_rec _att_data_rec;
	trl_radi_data_rec _radi_data_rec;
	trl_radi_comp_rec _radi_comp_rec;
};