#pragma once

#include "data-load.h"
#include "baryon-pion-block-acc.h"
#include "two-point-function-init&correct.h"
#include "baryon-pion-coefficient.h"
#include "compute-utils.h"

namespace qlat
{

    inline std::string get_two_point_correlation_function_pro_all_psel_func_psel_path(const std::string &job_tag,
                                                                                      const int &traj)
    {
        return ssprintf("analysis/two_point_correlation_function_pro/%s/results=%d", job_tag.c_str(),
                        traj);
    }

    inline void compute_meson_two_point_correlation_function_pro_all_psel_func_psel_type(
        const std::string &job_tag, const int &traj,
        const std::vector<std::vector<int>> &Momentum_Targets,
        const std::vector<ComplexD> &coefficent_of_t2pf,
        const std::vector<ComplexD> &coefficent_of_t2pf_remain_snk,
        const bool &is_meson_smear,
        const PointsSelection &psel_meson,
        const Long &psel_meson_num,
        const std::vector<int> &psel_meson_num_list,
        const std::vector<int> &list_n_from_idx_meson,
        const std::vector<std::vector<Long>> &idx_class_by_time_meson,
        const Geometry &geo,
        const Coordinate &total_site)
    {
        check_sigterm();
        check_time_limit();
        TIMER_VERBOSE("compute_meson_two_point_correlation_function_pro_all_psel_func_psel_type");
        Timer::autodisplay();

        const int &mom_num = Momentum_Targets.size();
        const std::string path = get_two_point_correlation_function_pro_all_psel_func_psel_path(job_tag, traj);

        std::string fn_meson_twop_list;
        std::string fn_meson_twop_list_rs;

        fn_meson_twop_list = path + ssprintf("/pi_%s_two_point_psel.lat", (is_meson_smear ? "smear" : "point"));
        fn_meson_twop_list_rs = path + ssprintf("/pi_%s_two_point_psel_snk_remain.lat", (is_meson_smear ? "smear" : "point"));

        std::vector<std::vector<std::vector<Meson_Block_Scalar_Psel_Psel>>> block_meson;
        contract_meson_psel_two_point_block(block_meson, psel_meson_num_list, list_n_from_idx_meson, job_tag, traj, geo, psel_meson, total_site[3] - 1, is_meson_smear);

        LatData meson_twop_data = mk_two_point_momentum_sum_table(total_site,Momentum_Targets);
        LatData meson_twop_data_rs = mk_two_point_momentum_remain_snk_table(total_site,Momentum_Targets);

        contract_meson_psel_two_point_block_all_psel_gys_omp_optimiztion_pro(
            meson_twop_data_rs,
            meson_twop_data,
            coefficent_of_t2pf,
            coefficent_of_t2pf_remain_snk,
            block_meson,
            job_tag, traj,
            Momentum_Targets,
            is_meson_smear,
            psel_meson,
            psel_meson_num,
            psel_meson_num_list,
            list_n_from_idx_meson,
            idx_class_by_time_meson,
            geo);

        lat_data_save_info(fn_meson_twop_list, meson_twop_data);
        lat_data_save_info(fn_meson_twop_list_rs, meson_twop_data_rs);
    }

    inline void compute_baryon_two_point_correlation_function_pro_all_psel_func_psel_type(
        const std::string &job_tag, const int &traj,
        const std::vector<std::vector<int>> &Momentum_Targets,
        const std::vector<ComplexD> &coefficent_of_t2pf,
        const std::vector<ComplexD> &coefficent_of_t2pf_remain_snk,
        const bool &is_baryon_smear,
        const PointsSelection &psel_baryon,
        const Long &psel_baryon_num,
        const std::vector<int> &psel_baryon_num_list,
        const std::vector<int> &list_n_from_idx_baryon,
        const std::vector<std::vector<Long>> &idx_class_by_time_baryon,
        const Geometry &geo,
        const Coordinate &total_site)
    {
        check_sigterm();
        check_time_limit();
        TIMER_VERBOSE("compute_baryon_two_point_correlation_function_pro_all_psel_func_psel_typ");
        Timer::autodisplay();

        const int &mom_num = Momentum_Targets.size();

        const std::string path = get_two_point_correlation_function_pro_all_psel_func_psel_path(job_tag, traj);

        std::string fn_baryon_twop_list;
        std::string fn_baryon_twop_list_rs;

        fn_baryon_twop_list = path + ssprintf("/baryon_%s_two_point_psel.lat", (is_baryon_smear ? "smear" : "point"));
        fn_baryon_twop_list_rs = path + ssprintf("/baryon_%s_two_point_psel_snk_remain.lat", (is_baryon_smear ? "smear" : "point"));


        std::vector<std::vector<std::vector<Baryon_Block_Scalar_Psel_Psel>>> block_baryon;
        contract_baryon_psel_two_point_block(block_baryon, psel_baryon_num_list, list_n_from_idx_baryon, job_tag, traj, geo, psel_baryon, total_site[3] - 1, is_baryon_smear);
        displayln_info("contract_two_point_or_pselpsel_block_finish;");

        LatData baryon_twop_data = mk_two_point_momentum_sum_table(total_site,Momentum_Targets);
        LatData baryon_twop_data_rs = mk_two_point_momentum_remain_snk_table(total_site,Momentum_Targets);

        contract_baryon_psel_two_point_block_all_psel_gys_omp_optimiztion_pro(
            baryon_twop_data_rs,
            baryon_twop_data,
            coefficent_of_t2pf,
            coefficent_of_t2pf_remain_snk,
            block_baryon,
            job_tag, traj,
            Momentum_Targets,
            is_baryon_smear,
            psel_baryon,
            psel_baryon_num,
            psel_baryon_num_list,
            list_n_from_idx_baryon,
            idx_class_by_time_baryon,
            geo);

        lat_data_save_info(fn_baryon_twop_list, baryon_twop_data);
        lat_data_save_info(fn_baryon_twop_list_rs, baryon_twop_data_rs);
    }

    inline void compute_baryon_no_projection_two_point_correlation_function_pro_all_psel_func_psel_type(
        const std::string &job_tag, const int &traj,
        const std::vector<std::vector<int>> &Momentum_Targets,
        const std::vector<ComplexD> &coefficent_of_t2pf,
        const std::vector<ComplexD> &coefficent_of_t2pf_remain_snk,
        const bool &is_baryon_smear,
        const PointsSelection &psel_baryon,
        const Long &psel_baryon_num,
        const std::vector<int> &psel_baryon_num_list,
        const std::vector<int> &list_n_from_idx_baryon,
        const std::vector<std::vector<Long>> &idx_class_by_time_baryon,
        const Geometry &geo,
        const Coordinate &total_site)
    {
        check_sigterm();
        check_time_limit();
        TIMER_VERBOSE("compute_baryon_no_projection_two_point_correlation_function_pro_all_psel_func_psel_type");
        Timer::autodisplay();

        const int &mom_num = Momentum_Targets.size();

        const std::string path = get_two_point_correlation_function_pro_all_psel_func_psel_path(job_tag, traj);

        std::string fn_baryon_twop_list_no_projection;
        std::string fn_baryon_twop_list_no_projection_rs;

        fn_baryon_twop_list_no_projection = path + ssprintf("/baryon_%s_two_point_no_projection_psel.lat", (is_baryon_smear ? "smear" : "point"));

        fn_baryon_twop_list_no_projection_rs = path + ssprintf("/baryon_%s_two_point_no_projection_psel_snk_remain.lat", (is_baryon_smear ? "smear" : "point"));


        std::vector<std::vector<std::vector<Baryon_Block_2PF_No_Projection_Psel_Psel>>> block_baryon_no_projection;
        contract_baryon_psel_two_point_no_projection_block(block_baryon_no_projection, psel_baryon_num_list, list_n_from_idx_baryon, job_tag, traj, geo, psel_baryon, total_site[3] - 1, is_baryon_smear);

        displayln_info("contract_two_point_or_pselpsel_block_finish;");

        LatData baryon_twop_data_no_projection = mk_two_point_momentum_no_projection_sum_table(total_site,Momentum_Targets);
        LatData baryon_twop_data_no_projection_rs = mk_two_point_momentum_no_projection_remain_snk_table(total_site,Momentum_Targets);

        contract_baryon_psel_two_point_block_no_projection_all_psel_gys_omp_optimiztion_pro(
            baryon_twop_data_no_projection_rs,
            baryon_twop_data_no_projection,
            coefficent_of_t2pf,
            coefficent_of_t2pf_remain_snk,
            block_baryon_no_projection,
            job_tag, traj,
            Momentum_Targets,
            is_baryon_smear,
            psel_baryon,
            psel_baryon_num,
            psel_baryon_num_list,
            list_n_from_idx_baryon,
            idx_class_by_time_baryon,
            geo);

        lat_data_save_info(fn_baryon_twop_list_no_projection, baryon_twop_data_no_projection);
        lat_data_save_info(fn_baryon_twop_list_no_projection_rs, baryon_twop_data_no_projection_rs);
    }

    inline void compute_two_point_correlation_function_pro_all_psel_func_psel_type(const std::string &job_tag, const int &traj, const std::vector<std::vector<int>> &Momentum_Targets, const bool is_baryon_smear, const bool is_meson_smear)
    {
        check_sigterm();
        check_time_limit();
        Timer::autodisplay();

        const int &mom_num = Momentum_Targets.size();

        TIMER_VERBOSE("compute_two_point_correlation_function_pro_all_psel_func_psel_type");
        const PointsSelection &psel = get_point_selection(job_tag, traj);
        const long n_points = psel.size();
        const PointsSelection &psel_smear = get_point_selection_smear(job_tag, traj);
        const long n_points_smear = psel_smear.size();

        const PointsSelection &psel_meson = (!is_meson_smear ? psel : psel_smear);
        const long &psel_meson_num = (!is_meson_smear ? n_points : n_points_smear);
        const PointsSelection &psel_baryon = (!is_baryon_smear ? psel : psel_smear);
        const long &psel_baryon_num = (!is_baryon_smear ? n_points : n_points_smear);

        const FieldSelection &fsel = get_field_selection(job_tag, traj);

        const Geometry &geo = fsel.f_rank.geo();
        const Coordinate &total_site = geo.total_site();

        std::vector<int> psel_baryon_num_list;
        std::vector<int> psel_meson_num_list;
        std::vector<int> list_n_from_idx_baryon(psel_baryon_num);
        std::vector<int> list_n_from_idx_meson(psel_meson_num);
        std::vector<std::vector<long>> idx_class_by_time_baryon(total_site[3]);
        std::vector<std::vector<long>> idx_class_by_time_meson(total_site[3]);

        baryon_pion_two_point_simple_standard(
            job_tag, traj,
            psel_baryon,
            psel_meson,
            psel_baryon_num_list,
            psel_meson_num_list,
            list_n_from_idx_baryon,
            list_n_from_idx_meson,
            idx_class_by_time_baryon,
            idx_class_by_time_meson,
            Momentum_Targets,
            geo,
            is_baryon_smear, is_meson_smear);

        std::vector<ComplexD> coefficent_of_t2pf;
        std::vector<ComplexD> coefficent_of_t2pf_remain_snk;

        contract_two_point_function_coefficent_pro(
            coefficent_of_t2pf,
            coefficent_of_t2pf_remain_snk,
            Momentum_Targets,
            is_baryon_smear,
            psel_baryon,
            psel_baryon_num,
            psel_baryon_num_list,
            list_n_from_idx_baryon,
            idx_class_by_time_baryon,
            is_meson_smear,
            psel_meson,
            psel_meson_num,
            psel_meson_num_list,
            list_n_from_idx_meson,
            idx_class_by_time_meson,
            total_site);
            
        displayln_info("Coefficent Calculation finished!!!!");

        compute_meson_two_point_correlation_function_pro_all_psel_func_psel_type(
            job_tag, traj,
            Momentum_Targets,
            coefficent_of_t2pf,
            coefficent_of_t2pf_remain_snk,
            is_meson_smear,
            psel_meson,
            psel_meson_num,
            psel_meson_num_list,
            list_n_from_idx_meson,
            idx_class_by_time_meson,
            geo,
            total_site);

        compute_baryon_two_point_correlation_function_pro_all_psel_func_psel_type(
            job_tag, traj,
            Momentum_Targets,
            coefficent_of_t2pf,
            coefficent_of_t2pf_remain_snk,
            is_baryon_smear,
            psel_baryon,
            psel_baryon_num,
            psel_baryon_num_list,
            list_n_from_idx_baryon,
            idx_class_by_time_baryon,
            geo,
            total_site);

        compute_baryon_no_projection_two_point_correlation_function_pro_all_psel_func_psel_type(
            job_tag, traj,
            Momentum_Targets,
            coefficent_of_t2pf,
            coefficent_of_t2pf_remain_snk,
            is_baryon_smear,
            psel_baryon,
            psel_baryon_num,
            psel_baryon_num_list,
            list_n_from_idx_baryon,
            idx_class_by_time_baryon,
            geo,
            total_site);
    }

    inline void compute_two_point_correlation_function_pro_all_psel_func_psel(const std::string &job_tag,
                                                                              const int &traj, const bool is_baryon_smear, const bool is_meson_smear)
    {
        check_sigterm();
        Timer::autodisplay();
        TIMER_VERBOSE("compute_two_point_correlation_function_pro_all_psel_func_psel");
        const std::string path = get_two_point_correlation_function_pro_all_psel_func_psel_path(job_tag, traj);
        const std::string path_checkpoint = path + "/checkpoint.txt";
        if (does_file_exist_sync_node(path_checkpoint))
        {
            return;
        }
        if (does_file_exist_sync_node(path + ssprintf("/baryon_%s_two_point_no_projection_psel.lat", (is_baryon_smear ? "smear" : "point"))))
        {
            return;
        }
        if (does_file_exist_sync_node(path + ssprintf("/baryon_%s_two_pointn_psel.lat", (is_baryon_smear ? "smear" : "point"))))
        {
            return;
        }
        if (does_file_exist_sync_node(path + ssprintf("/meson_%s_two_pointn_psel.lat", (is_meson_smear ? "smear" : "point"))))
        {
            return;
        }
        if (not(check_prop_smear(job_tag, traj, 0)))
        {
            displayln_info("err_check_prop_smear;");
            displayln_info(get_prop_smear_path(job_tag, traj, 0));
            displayln_info(get_psel_prop_smear_path(job_tag, traj, 0));
            return;
        }
        check_sigterm();
        check_time_limit();
        if (not obtain_lock(
                ssprintf("lock-baryon-pi-psel-%s-%d", job_tag.c_str(), traj)))
        {
            return;
        }
        setup(job_tag, traj);
        qmkdir_info("analysis/two_point_correlation_function_pro");
        qmkdir_info(ssprintf("analysis/two_point_correlation_function_pro/%s", job_tag.c_str()));
        qmkdir_info(path);

        std::vector<std::vector<int>> Momentum_Targets;

        Momentum_Targets.push_back({0, 0, 0});

        Momentum_Targets.push_back({0, 0, 1});
        Momentum_Targets.push_back({0, 0, -1});
        Momentum_Targets.push_back({0, 1, 0});
        Momentum_Targets.push_back({0, -1, 0});
        Momentum_Targets.push_back({1, 0, 0});
        Momentum_Targets.push_back({-1, 0, 0});

        Momentum_Targets.push_back({0, 1, 1});
        Momentum_Targets.push_back({0, 1, -1});
        Momentum_Targets.push_back({0, -1, 1});
        Momentum_Targets.push_back({0, -1, -1});
        Momentum_Targets.push_back({1, 0, 1});
        Momentum_Targets.push_back({1, 0, -1});
        Momentum_Targets.push_back({1, 1, 0});
        Momentum_Targets.push_back({1, -1, 0});
        Momentum_Targets.push_back({-1, 0, 1});
        Momentum_Targets.push_back({-1, 0, -1});
        Momentum_Targets.push_back({-1, 1, 0});
        Momentum_Targets.push_back({-1, -1, 0});

        Momentum_Targets.push_back({1, 1, 1});
        Momentum_Targets.push_back({1, 1, -1});
        Momentum_Targets.push_back({1, -1, 1});
        Momentum_Targets.push_back({1, -1, -1});
        Momentum_Targets.push_back({-1, 1, 1});
        Momentum_Targets.push_back({-1, 1, -1});
        Momentum_Targets.push_back({-1, -1, 1});
        Momentum_Targets.push_back({-1, -1, -1});

        Momentum_Targets.push_back({0, 0, 2});
        Momentum_Targets.push_back({0, 0, -2});
        Momentum_Targets.push_back({0, 2, 0});
        Momentum_Targets.push_back({0, -2, 0});
        Momentum_Targets.push_back({2, 0, 0});
        Momentum_Targets.push_back({-2, 0, 0});

        compute_two_point_correlation_function_pro_all_psel_func_psel_type(job_tag, traj, Momentum_Targets, is_baryon_smear, is_meson_smear);

        Timer::autodisplay();
        release_lock();
    }
}