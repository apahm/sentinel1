#pragma once
#include <iostream>
#include <string>
#include <vector>

struct ldr_sar_desc_rec {
    uint32_t rec_seq;
    uint8_t rec_sub1;
    uint8_t rec_type;
    uint8_t rec_sub2;
    uint8_t rec_sub3;
    uint32_t length;
    std::string ascii_flag;
    std::string spare1;
    std::string format_doc;
    std::string format_rev;
    std::string design_rev;
    std::string software_id;
    uint32_t file_num;
    std::string file_name;
    std::string rec_seq;
    uint64_t seq_loc;
    uint32_t seq_len;
    std::string rec_code;
    uint64_t code_loc;
    uint32_t code_len;
    std::string rec_len;
    uint64_t rlen_loc;
    uint32_t rlen_len;
    std::string spare2;
    std::string spare3;
    uint64_t n_dataset;
    uint64_t l_dataset;
    uint64_t n_map_proj;
    uint64_t l_map_proj;
    uint64_t n_plat_pos;
    uint64_t l_plat_pos;
    uint64_t n_att_data;
    uint64_t l_att_data;
    uint64_t n_radi_data;
    uint64_t l_radi_data;
    uint64_t n_radi_comp;
    uint64_t l_radi_comp;
    uint64_t n_qual_sum;
    uint64_t l_qual_sum;
    uint64_t n_data_hist;
    uint64_t l_data_hist;
    uint64_t n_rang_spec;
    uint64_t l_rang_spec;
    uint64_t n_dem_desc;
    uint64_t l_dem_desc;
    uint64_t n_radar_par;
    uint64_t l_radar_par;
    uint64_t n_anno_data;
    uint64_t l_anno_data;
    uint64_t n_det_proc;
    uint64_t l_det_proc;
    uint64_t n_cal;
    uint64_t l_cal;
    uint64_t n_gcp;
    uint64_t l_gcp;
    std::vector<uint64_t> spare4;
    uint64_t n_fac_data;
    uint64_t l_fac_data;
    std::string spare5;
};

struct ldr_dataset_sum_rec {
    uint32_t rec_seq;
    uint8_t rec_sub1;
    uint8_t rec_type;
    uint8_t rec_sub2;
    uint8_t rec_sub3;
    uint32_t Length;
    uint32_t seq_num;
    uint32_t sar_chn;
    std::string scene_id;
    std::string scene_des;
    std::string inp_sctim;
    std::string asc_des;
    double pro_lat;
    double pro_long;
    double pro_head;
    std::string ellip_des;
    double ellip_maj;
    double ellip_min;
    double earth_mass;
    double grav_const;
    double ellip_j;
    std::string spare2;
    double terrain_h;
    uint64_t sc_lin;
    uint64_t sc_pix;
    double scene_len;
    double scene_wid;
    std::string spare3;
    uint32_t nchn;
    std::string spare5;
    std::string mission_id;
    std::string sensor_id;
    std::string orbit_num;
    double plat_lat;
    double plat_long;
    double plat_head;
    double clock_ang;
    double incident_ang;
    std::string spare15;
    double wave_length;
    std::string motion_comp;
    std::string pulse_code;
    double ampl_coef;
    double phas_coef;
    uint64_t chirp_ext_ind;
    std::string spare6;
    double fr;
    double rng_gate;
    double rng_length;
    std::string baseband_f;
    std::string rngcmp_f;
    double gn_polar;
    double gn_cross;
    uint64_t chn_bits;
    std::string quant_desc;
    double i_bias;
    double q_bias;
    double iq_ratio;
    double spare7;
    double spare8;
    double ele_sight;
    double mech_sight;
    std::string echo_track;
    double fa;
    double elev_beam;
    double azim_beam;
    std::vector<uint8_t> sat_bintim;
    std::vector<uint8_t> sat_clktim;
    uint64_t sat_clkinc;
    std::string spare9;
    std::string fac_id;
    std::string sys_id;
    std::string ver_id;
    std::string fac_code;
    std::string lev_code;
    std::string prod_type;
    std::string algor_id;
    double n_azilok;
    double n_rnglok;
    double bnd_azilok;
    double bnd_rnglok;
    double bnd_azi;
    double bnd_rng;
    std::string azi_weight;
    std::string rng_weight;
    std::string data_inpsrc;
    double rng_res;
    double azi_res;
    double radi_stretch;
    double alt_dopcen;
    std::string spare10;
    double crt_dopcen;
    std::string time_dir_pix;
    std::string time_dir_lin;
    double alt_rate;
    std::string spare12;
    double crt_rate;
    std::string spare13;
    std::string line_cont;
    std::string clutter_lock;
    std::string auto_focus;
    double line_spacing;
    double pix_spacing;
    std::string rngcmp_desg;
    std::string spare14;
};

struct ldr_qual_sum_rec {
    uint32_t rec_seq;
    uint8_t rec_sub1;
    uint8_t rec_type;
    uint8_t rec_sub2;
    uint8_t rec_sub3;
    uint32_t length;
    uint32_t rec_seq;
    std::string sar_chn;
    std::string cali_date;
    uint32_t nchn;
    double islr;
    double pslr;
    double azi_ambig;
    double rng_ambig;
    double snr;
    double ber;
    double rng_res;
    double azi_res;
    double rad_res;
    double dyn_rng;
    double rad_unc_db;
    double rad_unc_deg;
    double db;
    double deg;
    double alt_locerr;
    double crt_locerr;
    double alt_scale;
    double crt_scale;
    double dis_skew;
    double ori_err;
    double alt_m;
    double crt_m;
    double nesz;
    double enl;
    std::string tb_update;
    std::string spare;
};

struct ldr_sdr_hist_rec {
    uint32_t rec_seq;
    uint8_t rec_sub1;
    uint8_t rec_type;
    uint8_t rec_sub2;
    uint8_t rec_sub3;
    uint32_t length;
    uint32_t rec_seq;
    uint32_t sar_chn;
    uint64_t ntab;
    uint64_t ltab;
    std::string hist_desc;
    uint32_t nrec;
    uint32_t tab_seq;
    uint64_t nbin;
    uint64_t ns_lin;
    uint64_t ns_pix;
    uint64_t ngrp_lin;
    uint64_t ngrp_pix;
    uint64_t nsamp_lin;
    uint64_t nsamp_pix;
    double min_smp;
    double max_smp;
    double mean_smp;
    double std_smp;
    double smp_inc;
    double min_hist;
    double max_hist;
    double mean_hist;
    double std_hist;
    uint64_t nhist;
    std::vector<uint64_t> hist;
    std::string spare;
};

struct ldr_pdr16_hist_rec {
    uint32_t rec_seq;
    uint8_t rec_sub1;
    uint8_t rec_type;
    uint8_t rec_sub2;
    uint8_t rec_sub3;
    uint32_t length;
    uint32_t rec_seq;
    uint32_t sar_chn;
    uint64_t ntab;
    uint64_t ltab;
    std::string hist_desc;
    uint32_t nrec;
    uint32_t tab_seq;
    uint64_t nbin;
    uint64_t ns_lin;
    uint64_t ns_pix;
    uint64_t ngrp_lin;
    uint64_t ngrp_pix;
    uint64_t nsamp_lin;
    uint64_t nsamp_pix;
    double min_smp;
    double max_smp;
    double mean_smp;
    double std_smp;
    double smp_inc;
    double min_hist;
    double max_hist;
    double mean_hist;
    double std_hist;
    uint64_t nhist;
    std::vector<uint64_t> hist;
    std::string hist_desc;
    uint32_t nrec;
    uint32_t tab_seq;
    uint64_t nbin;
    uint64_t ns_lin;
    uint64_t ns_pix;
    uint64_t ngrp_lin;
    uint64_t ngrp_pix;
    uint64_t nsamp_lin;
    uint64_t nsamp_pix;
    double min_smp;
    double max_smp;
    double mean_smp;
    double std_smp;
    double smp_inc;
    double min_hist;
    double max_hist;
    double mean_hist;
    double std_hist;
    uint64_t nhist;
    std::vector<uint64_t> hist;
    std::string spare;
};

struct ldr_proc_parm_rec {

};

struct ldr_map_proj_rec {

};

struct ldr_pos_data_rec {

};

struct ldr_att_data_rec {

};

struct ldr_radi_data_rec {

};

struct ldr_radi_comp_rec {

};

struct sar_ldr_file {
    ldr_sar_desc_rec _sar_desc_rec;
    ldr_dataset_sum_rec _dataset_sum_rec;
    ldr_qual_sum_rec _qual_sum_rec;
    ldr_sdr_hist_rec _sdr_hist_rec;
    ldr_pdr16_hist_rec _pdr16_hist_rec;
    ldr_proc_parm_rec _proc_parm_rec;
    ldr_map_proj_rec _map_proj_rec;
    ldr_pos_data_rec _pos_data_rec;
    ldr_att_data_rec _att_data_rec;
    ldr_radi_data_rec _radi_data_rec;
    ldr_radi_comp_rec _radi_comp_rec;
};