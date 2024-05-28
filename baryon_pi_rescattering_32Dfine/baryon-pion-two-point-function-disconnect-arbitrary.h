#pragma once

#include "data-load.h"
#include "baryon-pion-block-acc.h"
#include "baryon-pion-coefficient.h"
#include "compute-utils.h"

namespace qlat
{

    inline std::string get_baryon_pi_rescattering_disconnect_all_psel_arbitrarty_func_psel_path(const std::string &job_tag,
                                                                                                const int traj)
    {
        return ssprintf("analysis/baryon_pi_rescattering_disconnect/%s/results=%d", job_tag.c_str(),
                        traj);
    }

    inline LatData mk_baryon_pi_rescattering_for_arbitrary_disconnect(const int &num_dtxy, const int &gram_num, const std::vector<std::vector<int>> &Momentum_Targets, const Coordinate &total_site)
    {
        LatData ld;
        ld.info.push_back(lat_dim_number("txy", 0, num_dtxy - 1));
        ld.info.push_back(lat_dim_number("momentum_mode_snk", 0, Momentum_Targets.size() - 1));
        ld.info.push_back(lat_dim_number("momentum_mode_src", 0, Momentum_Targets.size() - 1));
        ld.info.push_back(lat_dim_number("gram", 0, gram_num - 1));
        ld.info.push_back(lat_dim_number("mu", 0, 3));
        ld.info.push_back(lat_dim_number("nu", 0, 3));
        // ld.info.push_back(lat_dim_number("drawers", 0, total_site[0] - 1));
        ld.info.push_back(lat_dim_re_im());
        lat_data_alloc(ld);
        set_zero(ld);
        return ld;
    }

    inline void
    contract_baryon_pi_rescattering_disconnect_all_psel_acc_psel(
        std::vector<std::vector<ComplexD>> &disconnect_result,
        const std::vector<std::vector<std::vector<Meson_Block_Scalar_Psel_Psel>>> &block_meson,
        const std::vector<std::vector<std::vector<Baryon_Block_2PF_No_Projection_Psel_Psel>>> &block_baryon,
        const std::vector<std::vector<std::vector<std::vector<std::vector<ComplexD>>>>> &baryon_pi_two_func_coeff,
        const std::vector<std::vector<int>> &Momentum_Targets,
        const bool &is_baryon_smear,
        const PointsSelection &psel_baryon,
        const long &psel_baryon_num,
        const std::vector<int> &psel_baryon_num_list,
        const std::vector<int> &list_n_from_idx_baryon,
        const std::vector<std::vector<long>> &idx_class_by_time_baryon,
        const bool &is_meson_smear,
        const PointsSelection &psel_meson,
        const long &psel_meson_num,
        const std::vector<int> &psel_meson_num_list,
        const std::vector<int> &list_n_from_idx_meson,
        const std::vector<std::vector<long>> &idx_class_by_time_meson,
        const Geometry &geo, const int &ti, const int &dtmax, const int &dtmin,
        const int &t_baryon_snk, const int &xt, const int &t_meson, const int &t_baryon_src, const int &dtxy,
        const double &anti_period, const int &gram_num)
    {
        TIMER_VERBOSE("contract_baryon_pi_rescattering_disconnect_all_psel_acc_psel");
        const Coordinate &total_site = geo.total_site();

        const int &mom_num = Momentum_Targets.size();

#pragma omp parallel for collapse(4)
        for (int idx_baryon_src = 0; idx_baryon_src < psel_baryon_num_list[t_baryon_src]; idx_baryon_src++)
        {
            for (int idx_baryon_snk = 0; idx_baryon_snk < psel_baryon_num_list[t_baryon_snk]; idx_baryon_snk++)
            {
                for (int idx_meson_ty = 0; idx_meson_ty < psel_meson_num_list[t_meson]; idx_meson_ty++)
                {
                    for (int idx_x = 0; idx_x < psel_meson_num_list[xt]; idx_x++)
                    {
                        const long &xg_src_psel_idx = idx_class_by_time_baryon[t_baryon_src][idx_baryon_src];
                        const Coordinate &xg_src = psel_baryon[xg_src_psel_idx];

                        const long &xg_snk_psel_idx = idx_class_by_time_baryon[t_baryon_snk][idx_baryon_snk];
                        const Coordinate &xg_snk = psel_baryon[xg_snk_psel_idx];

                        const long &xg_y_psel_idx = idx_class_by_time_meson[t_meson][idx_meson_ty];
                        const Coordinate &xg_y = psel_meson[xg_y_psel_idx];

                        const long &xg_x_psel_idx = idx_class_by_time_meson[xt][idx_x];
                        const Coordinate &xg_x = psel_meson[xg_x_psel_idx];

                        if (is_baryon_smear == is_meson_smear && (xg_snk_psel_idx == xg_y_psel_idx || xg_src_psel_idx == xg_x_psel_idx))
                        {
                            continue;
                        }

                        if (is_baryon_smear == is_meson_smear && (xg_snk_psel_idx == xg_x_psel_idx || xg_src_psel_idx == xg_y_psel_idx))
                        {
                            continue;
                        }

                        const int del_t_snk_src = mod(t_baryon_snk - t_baryon_src, total_site[3]);
                        const int del_t_x_y = mod(xt - t_meson, total_site[3]);
                        qassert(del_t_snk_src == dtxy);
                        qassert(del_t_x_y == dtxy);

                        const Meson_Block_Scalar_Psel_Psel &block_meson_temp = block_meson[xg_x_psel_idx][del_t_snk_src][idx_meson_ty];
                        const Baryon_Block_2PF_No_Projection_Psel_Psel &block_baryon_temp = block_baryon[xg_snk_psel_idx][del_t_x_y][idx_baryon_src];

                        qassert(block_meson_temp.is_build);
                        qassert(block_baryon_temp.is_build);

                        qassert(list_n_from_idx_baryon[xg_src_psel_idx] == idx_baryon_src);
                        qassert(list_n_from_idx_baryon[xg_snk_psel_idx] == idx_baryon_snk);
                        qassert(list_n_from_idx_meson[xg_y_psel_idx] == idx_meson_ty);
                        qassert(list_n_from_idx_meson[xg_x_psel_idx] == idx_x);

                        qassert(block_baryon_temp.xg_snk_psel_idx == xg_snk_psel_idx);
                        qassert(block_baryon_temp.xg_src_psel_idx == xg_src_psel_idx);
                        qassert(block_meson_temp.xg_snk_psel_idx == xg_x_psel_idx);
                        qassert(block_meson_temp.xg_src_psel_idx == xg_y_psel_idx);
                        
                        qassert(block_baryon_temp.n_snk == idx_baryon_snk);
                        qassert(block_baryon_temp.n_src == idx_baryon_src);
                        qassert(block_meson_temp.n_snk == idx_x);
                        qassert(block_meson_temp.n_src == idx_meson_ty);


                        std::vector<SpinMatrix> spin_result_list_temp;
                        spin_result_list_temp.resize(2);

                        spin_result_list_temp[0] = block_baryon_temp.block_pp[0] * block_meson_temp.block_pp;
                        spin_result_list_temp[1] = block_baryon_temp.block_pp[1] * block_meson_temp.block_pp;

                        for (int mom_src_id = 0; mom_src_id < mom_num; ++mom_src_id)
                        {
                            for (int mom_snk_id = 0; mom_snk_id < mom_num; ++mom_snk_id)
                            {
                                long phase_temp = space_dot(Momentum_Targets[mom_src_id], xg_src) - space_dot(Momentum_Targets[mom_src_id], xg_y) + space_dot(Momentum_Targets[mom_snk_id], xg_snk) - space_dot(Momentum_Targets[mom_snk_id], xg_x);
                                int phase = mod(phase_temp, total_site[0]);

                                const ComplexD normal_coeff = baryon_pi_two_func_coeff[ti][mom_src_id][mom_snk_id][dtxy][phase] * (ComplexD)anti_period;

                                for (int topology_id = 0; topology_id < gram_num; topology_id++)
                                {
                                    SpinMatrix spin_result_list = spin_result_list_temp[topology_id] * normal_coeff;
                                    // std::cout<<spin_result_list(0, 0)<<std::endl;
                                    for (int mu1 = 0; mu1 < 4; mu1++)
                                    {
                                        for (int mu2 = 0; mu2 < 4; mu2++)
                                        {
                                            if (spin_result_list(mu1, mu2) == 0.0)
                                            {
                                                continue;
                                            }
                                            disconnect_result[omp_get_thread_num()][(((mom_src_id * mom_num + mom_snk_id) * gram_num + topology_id) * 4 + mu1) * 4 + mu2] += spin_result_list(mu1, mu2);
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    inline void acc_data_to_disconnect_LatData(
        std::vector<std::vector<ComplexD>> &disconnect_result,
        LatData &baryon_pi_data_temp,
        const std::vector<std::vector<std::vector<std::vector<std::vector<ComplexD>>>>> &baryon_pi_two_func_coeff,
        const std::vector<std::vector<int>> &Momentum_Targets,
        const bool &is_baryon_smear,
        const PointsSelection &psel_baryon,
        const long &psel_baryon_num,
        const std::vector<int> &psel_baryon_num_list,
        const std::vector<int> &list_n_from_idx_baryon,
        const std::vector<std::vector<long>> &idx_class_by_time_baryon,
        const bool &is_meson_smear,
        const PointsSelection &psel_meson,
        const long &psel_meson_num,
        const std::vector<int> &psel_meson_num_list,
        const std::vector<int> &list_n_from_idx_meson,
        const std::vector<std::vector<long>> &idx_class_by_time_meson,
        const Geometry &geo, const int &ti, const int &dtmax, const int &dtmin,
        const int &t_baryon_snk, const int &xt, const int &t_meson, const int &t_baryon_src, const int &dtxy, const int &gram_num)
    {
        TIMER_VERBOSE("acc_data_to_disconnect_LatData");

        const Coordinate &total_site = geo.total_site();
        const int &mom_num = Momentum_Targets.size();

#pragma omp parallel for collapse(2)
        for (int mom_src_id = 0; mom_src_id < mom_num; ++mom_src_id)
        {
            for (int mom_snk_id = 0; mom_snk_id < mom_num; ++mom_snk_id)
            {
                Vector<ComplexD> baryon_pi_data_vec_temp = lat_data_complex_get(baryon_pi_data_temp, make_array(dtxy, mom_snk_id, mom_src_id));
                for (int thread_id = 0; thread_id < omp_get_max_threads(); thread_id++)
                {
                    for (int topology_id = 0; topology_id < gram_num; topology_id++)
                    {
                        for (int mu1 = 0; mu1 < 4; mu1++)
                        {
                            for (int mu2 = 0; mu2 < 4; mu2++)
                            {
                                baryon_pi_data_vec_temp[(topology_id * 4 + mu1) * 4 + mu2] += disconnect_result[thread_id][(((mom_src_id * mom_num + mom_snk_id) * gram_num + topology_id) * 4 + mu1) * 4 + mu2];
                            }
                        }
                    }
                }
            }
        }
    }

    inline void compute_baryon_pi_rescattering_disconnect_all_psel_arbitrarty_func_psel_type(const std::string &job_tag, const int &traj, const std::vector<std::vector<int>> &Momentum_Targets, const std::vector<int> &tsep_src_list, const std::vector<int> &tsep_snk_list, const bool is_baryon_smear, const bool is_meson_smear)
    {
        check_sigterm();
        check_time_limit();
        Timer::autodisplay();
        const int gram_num = 2;

        qassert(tsep_src_list.size() == tsep_snk_list.size());
        const int tsep_num = (int)tsep_src_list.size();

        const std::string path = get_baryon_pi_rescattering_disconnect_all_psel_arbitrarty_func_psel_path(job_tag, traj);
        std::vector<std::string> fn_disconnect_gram_list(tsep_num);
        for (int ti = 0; ti < tsep_num; ++ti)
        {
            const int &tsep_src = tsep_src_list[ti];
            const int &tsep_snk = tsep_snk_list[ti];
            fn_disconnect_gram_list[ti] = path + ssprintf("/baryon_%s_pi_%s_rescattering_disconnect_tsep_%d_%d.field", (is_baryon_smear ? "smear" : "point"), (is_meson_smear ? "smear" : "point"), tsep_src, tsep_snk);
        }

        bool is_complete = true;
        for (int ti = 0; ti < tsep_num; ++ti)
        {
            if (fn_disconnect_gram_list[ti] != "")
            {
                if (not is_d_field(fn_disconnect_gram_list[ti]))
                {
                    is_complete = false;
                    break;
                }
            }
        }

        if (is_complete)
        {
            return;
        }
        TIMER_VERBOSE("compute_baryon_pi_rescattering_all_psel_arbitrarty_func_psel_type");
        const PointsSelection &psel = get_point_selection(job_tag, traj);
        const long n_points = psel.size();
        const PointsSelection &psel_smear = get_point_selection_smear(job_tag, traj);
        const long n_points_smear = psel_smear.size();

        const PointsSelection &psel_meson = (!is_meson_smear ? psel : psel_smear);
        const long &psel_meson_num = (!is_meson_smear ? n_points : n_points_smear);
        const PointsSelection &psel_baryon = (!is_baryon_smear ? psel : psel_smear);
        const long &psel_baryon_num = (!is_baryon_smear ? n_points : n_points_smear);

        const FieldSelection &fsel = get_field_selection(job_tag, traj);
        std::vector<LatData> baryon_pi_disconnect_data_result(tsep_num);

        const Geometry &geo = fsel.f_rank.geo();
        const Coordinate total_site = geo.total_site();
        const int num_dtxy = get_baryon_pi_half_sep(total_site[3]);
        const int &dtmax = num_dtxy + (tsep_snk_list.back() > 0 ? tsep_snk_list.back() : 0) + (tsep_src_list.back() > 0 ? tsep_src_list.back() : 0);
        const int &dtmin = +(tsep_snk_list.front() > 0 ? tsep_snk_list.front() : 0) + (tsep_src_list.front() > 0 ? tsep_src_list.front() : 0); // block 所需时间片长度

        for (int ti = 0; ti < tsep_num; ++ti)
        {
            baryon_pi_disconnect_data_result[ti] = mk_baryon_pi_rescattering_for_arbitrary_disconnect(num_dtxy, gram_num, Momentum_Targets, total_site);
        }

        std::vector<std::vector<std::vector<std::vector<std::vector<ComplexD>>>>> baryon_pi_two_func_coeff;
        std::vector<int> psel_baryon_num_list;
        std::vector<int> psel_meson_num_list;
        std::vector<int> list_n_from_idx_baryon(psel_baryon_num);
        std::vector<int> list_n_from_idx_meson(psel_meson_num);
        std::vector<std::vector<long>> idx_class_by_time_baryon(psel_baryon_num);
        std::vector<std::vector<long>> idx_class_by_time_meson(psel_meson_num);

        baryon_pion_two_point_arbitrary_coeffience(
            job_tag, traj, baryon_pi_two_func_coeff,
            psel_baryon_num_list,
            psel_meson_num_list,
            list_n_from_idx_baryon,
            list_n_from_idx_meson,
            idx_class_by_time_baryon,
            idx_class_by_time_meson,
            Momentum_Targets,
            geo,
            psel_baryon, psel_meson,
            tsep_src_list, tsep_snk_list,
            is_baryon_smear, is_meson_smear);

        displayln_info("baryon_pion_two_point_coeffience_finish;");

        std::vector<std::vector<std::vector<Meson_Block_Scalar_Psel_Psel>>> block_meson;
        std::vector<std::vector<std::vector<Baryon_Block_2PF_No_Projection_Psel_Psel>>> block_baryon;

        contract_meson_psel_two_point_block(block_meson, psel_meson_num_list, list_n_from_idx_meson, job_tag, traj, geo, psel_meson, total_site[3] - 1, is_meson_smear);
        contract_baryon_psel_two_point_no_projection_block(block_baryon, psel_baryon_num_list, list_n_from_idx_baryon, job_tag, traj, geo, psel_baryon, total_site[3] - 1, is_baryon_smear);

        displayln_info("contract_two_point_or_pselpsel_block_finish;");

        Timer::autodisplay();

        // 获取通信域中的进程总数（节点数）
        int MPI_SIZE;
        MPI_Comm_size(MPI_COMM_WORLD, &MPI_SIZE);

        // 获取当前进程的进程号
        int MPI_RANK;
        MPI_Comm_rank(MPI_COMM_WORLD, &MPI_RANK);

        const int y_idx_sparsity = 1;
        for (int ti = 0; ti < tsep_num; ++ti)
        {
            const int &tsep_src = tsep_src_list[ti];
            const int &tsep_snk = tsep_snk_list[ti];
            long iter = 0;
            for (int t_meson = 0; t_meson < total_site[3]; ++t_meson)
            {
                if (t_meson % MPI_SIZE != MPI_RANK)
                {
                    continue;
                }
                for (int dt = 0; dt < num_dtxy; ++dt)
                {
                    const int &dtxy = dt;
                    int xt, t_baryon_src, t_baryon_snk;
                    double anti_period;
                    four_point_time_slice(dtxy, num_dtxy, tsep_src, tsep_snk, total_site, t_meson, t_baryon_src, t_baryon_snk, xt, anti_period);

                    Timer::autodisplay();
                    TIMER_VERBOSE("compute_baryon_pi_rescattering_disconnect_all_psel_func_noama_type-iter");
                    iter += 1;
                    displayln_info(fname + ssprintf(":iter=%ld ti=%ld tmeson=%ld dt=%ld", iter, ti, t_meson, dt));

                    std::vector<std::vector<ComplexD>> disconnect_result;
                    const int &mom_num = Momentum_Targets.size();

                    disconnect_result.resize(omp_get_max_threads());
                    for (int thread_id = 0; thread_id < omp_get_max_threads(); thread_id++)
                    {
                        disconnect_result[thread_id].resize(mom_num * mom_num * gram_num * 4 * 4);
                        set_zero(disconnect_result[thread_id]);
                    }

                    contract_baryon_pi_rescattering_disconnect_all_psel_acc_psel(disconnect_result,
                                                                                 block_meson,
                                                                                 block_baryon,
                                                                                 baryon_pi_two_func_coeff,
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
                                                                                 geo, ti, dtmax, dtmin,
                                                                                 t_baryon_snk, xt, t_meson, t_baryon_src, dtxy, anti_period, gram_num);
                    acc_data_to_disconnect_LatData(disconnect_result,
                                                   baryon_pi_disconnect_data_result[ti],
                                                   baryon_pi_two_func_coeff,
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
                                                   geo, ti, dtmax, dtmin,
                                                   t_baryon_snk, xt, t_meson, t_baryon_src, dtxy, gram_num);
                }
            }
            Timer::autodisplay();
            glb_sum(baryon_pi_disconnect_data_result[ti]);
            lat_data_save_info(fn_disconnect_gram_list[ti], baryon_pi_disconnect_data_result[ti]);
        }
    }

    inline void compute_baryon_pi_rescattering_disconnect_all_psel_arbitrarty_func_psel(const std::string &job_tag,
                                                                                        const int traj, const bool is_baryon_smear, const bool is_meson_smear)
    {
        check_sigterm();
        Timer::autodisplay();
        TIMER_VERBOSE("compute_baryon_pi_rescattering_disconnect_all_psel_arbitrarty_func_psel");
        const std::string path = get_baryon_pi_rescattering_disconnect_all_psel_arbitrarty_func_psel_path(job_tag, traj);
        const std::string path_checkpoint = path + "/checkpoint.txt";
        if (does_file_exist_sync_node(path_checkpoint))
        {
            return;
        }
        if (does_file_exist_sync_node(path + ssprintf("/baryon_%s_pi_%s_rescattering_disconnect_tsep_%d_%d.field", (is_baryon_smear ? "smear" : "point"), (is_meson_smear ? "smear" : "point"), 0, 0)))
        {
            return;
        }
        if (not(check_prop_smear_psel(job_tag, traj, 0)))
        {
            displayln_info("err_check_prop_smear;");
            displayln_info(get_prop_smear_path(job_tag, traj, 0));
            displayln_info(get_psel_prop_smear_path(job_tag, traj, 0));
            return;
        }
        check_sigterm();
        check_time_limit();
        if (not obtain_lock(
                ssprintf("lock-baryon-pi-disconnect-%s-%d", job_tag.c_str(), traj)))
        {
            return;
        }
        setup(job_tag, traj);
        qmkdir_info("analysis/baryon_pi_rescattering_disconnect");
        qmkdir_info(ssprintf("analysis/baryon_pi_rescattering_disconnect/%s", job_tag.c_str()));
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

        std::vector<int> tsep_src_list;
        tsep_src_list.push_back(0);

        std::vector<int> tsep_snk_list;
        tsep_snk_list.push_back(0);

        compute_baryon_pi_rescattering_disconnect_all_psel_arbitrarty_func_psel_type(job_tag, traj, Momentum_Targets, tsep_src_list, tsep_snk_list, is_baryon_smear, is_meson_smear);

        release_lock();
    }
}
