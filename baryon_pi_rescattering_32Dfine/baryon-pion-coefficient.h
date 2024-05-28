#pragma once

#include "baryon-pion-block-acc.h"

namespace qlat
{
    inline void baryon_pion_two_point_simple(
        const std::string &job_tag, const int &traj,
        std::vector<int> &psel_baryon_num_list,
        std::vector<int> &psel_meson_num_list,
        std::vector<int> &list_n_from_idx_baryon,
        std::vector<int> &list_n_from_idx_meson,
        const Geometry &geo, const PointsSelection &psel_baryon, const PointsSelection &psel_meson,
        const bool is_baryon_smear, const bool is_meson_smear)
    {
        TIMER_VERBOSE("baryon_pion_two_point_coeffience");
        const Coordinate &total_site = geo.total_site();
        const int num_dtxy = get_baryon_pi_max_sep(total_site[3]);

        const long &psel_meson_num = psel_meson.size();
        const long &psel_baryon_num = psel_baryon.size();

        psel_baryon_num_list.resize(total_site[3]);
        psel_meson_num_list.resize(total_site[3]);
        set_zero(psel_baryon_num_list);
        set_zero(psel_meson_num_list);

        list_n_from_idx_baryon.resize(psel_baryon_num);
        list_n_from_idx_meson.resize(psel_meson_num);
        set_zero(list_n_from_idx_baryon);
        set_zero(list_n_from_idx_meson);

        for (long n = 0; n < psel_meson_num; ++n)
        {
            const long &xg_psel_idx = n;
            const Coordinate xg = psel_meson[xg_psel_idx];
            const int t = xg[3];
            list_n_from_idx_meson[xg_psel_idx] = psel_meson_num_list[t];
            psel_meson_num_list[t] += 1;
        }

        for (long n = 0; n < psel_baryon_num; ++n)
        {
            const long &xg_psel_idx = n;
            const Coordinate xg = psel_baryon[xg_psel_idx];
            const int t = xg[3];
            list_n_from_idx_baryon[xg_psel_idx] = psel_baryon_num_list[t];
            psel_baryon_num_list[t] += 1;
        }
    }

    inline void baryon_pion_two_point_coeffience(
        const std::string &job_tag, const int &traj,
        std::vector<std::vector<long>> &baryon_pi_two_func_coeff,
        std::vector<int> &psel_baryon_num_list,
        std::vector<int> &psel_meson_num_list,
        std::vector<int> &list_n_from_idx_baryon,
        std::vector<int> &list_n_from_idx_meson,
        const Geometry &geo, const PointsSelection &psel_baryon, const PointsSelection &psel_meson,
        const std::vector<int> &tsep_list, const bool is_baryon_smear, const bool is_meson_smear)
    {
        TIMER_VERBOSE("baryon_pion_two_point_coeffience");
        const int tsep_num = (int)tsep_list.size();
        const Coordinate &total_site = geo.total_site();
        const int num_dtxy = get_baryon_pi_max_sep(total_site[3]);

        const long &psel_meson_num = psel_meson.size();
        const long &psel_baryon_num = psel_baryon.size();

        psel_baryon_num_list.resize(total_site[3]);
        psel_meson_num_list.resize(total_site[3]);
        set_zero(psel_baryon_num_list);
        set_zero(psel_meson_num_list);

        list_n_from_idx_baryon.resize(psel_baryon_num);
        list_n_from_idx_meson.resize(psel_meson_num);
        set_zero(list_n_from_idx_baryon);
        set_zero(list_n_from_idx_meson);

        for (long n = 0; n < psel_meson_num; ++n)
        {
            const long &xg_psel_idx = n;
            const Coordinate xg = psel_meson[xg_psel_idx];
            const int t = xg[3];
            list_n_from_idx_meson[xg_psel_idx] = psel_meson_num_list[t];
            psel_meson_num_list[t] += 1;
        }

        for (long n = 0; n < psel_baryon_num; ++n)
        {
            const long &xg_psel_idx = n;
            const Coordinate xg = psel_baryon[xg_psel_idx];
            const int t = xg[3];
            list_n_from_idx_baryon[xg_psel_idx] = psel_baryon_num_list[t];
            psel_baryon_num_list[t] += 1;
        }

        baryon_pi_two_func_coeff.resize(tsep_num);
        for (int ti = 0; ti < tsep_num; ++ti)
        {
            baryon_pi_two_func_coeff[ti].resize(num_dtxy);
            set_zero(baryon_pi_two_func_coeff[ti]);
        }
        for (int ti = 0; ti < tsep_num; ++ti)
        {
            const int tsep = tsep_list[ti];

            for (int t_meson = 0; t_meson < total_site[3]; t_meson++)
            {
                for (int dt = 0; dt < num_dtxy; ++dt)
                {
                    const int dtxy = dt - num_dtxy / 2;
                    int xt = (tsep >= 0) ? mod(t_meson + dtxy, total_site[3]) : (dtxy >= 0 ? mod(t_meson - 2 * tsep + dtxy, total_site[3]) : mod(t_meson + 2 * tsep + dtxy, total_site[3]));
                    int t_baryon_src, t_baryon_snk;
                    if (dtxy >= 0)
                    {
                        t_baryon_src = mod(t_meson - tsep, total_site[3]);
                        t_baryon_snk = mod(xt + tsep, total_site[3]);
                    }
                    else
                    {
                        t_baryon_src = mod(xt - tsep, total_site[3]);
                        t_baryon_snk = mod(t_meson + tsep, total_site[3]);
                    }
                    if (tsep == 0 && is_baryon_smear == is_meson_smear)
                    {
                        baryon_pi_two_func_coeff[ti][dt] += psel_baryon_num_list[t_baryon_src] * psel_baryon_num_list[t_baryon_snk] * (psel_meson_num_list[xt] - 1) * (psel_meson_num_list[t_meson] - 1);
                    }
                    else
                    {
                        baryon_pi_two_func_coeff[ti][dt] += psel_baryon_num_list[t_baryon_src] * psel_baryon_num_list[t_baryon_snk] * psel_meson_num_list[xt] * psel_meson_num_list[t_meson];
                    }
                }
            }
        }
    }

    inline void baryon_pion_two_point_fsel4meson_coeffience(
        const std::string &job_tag, const int &traj,
        std::vector<std::vector<long>> &baryon_pi_two_func_coeff,
        std::vector<int> &psel_baryon_num_list,
        std::vector<int> &psel_meson_num_list,
        std::vector<int> &list_n_from_idx_baryon,
        std::vector<int> &list_n_from_idx_meson,
        const Geometry &geo, const PointsSelection &psel_baryon, const PointsSelection &psel_meson,
        const std::vector<int> &tsep_list, const bool is_baryon_smear, const bool is_meson_smear)
    {
        TIMER_VERBOSE("baryon_pion_two_point_coeffience");
        const int tsep_num = (int)tsep_list.size();
        const Coordinate &total_site = geo.total_site();
        const int num_dtxy = get_baryon_pi_max_sep(total_site[3]);

        const long &psel_meson_num = psel_meson.size();
        const long &psel_baryon_num = psel_baryon.size();

        psel_baryon_num_list.resize(total_site[3]);
        psel_meson_num_list.resize(total_site[3]);
        set_zero(psel_baryon_num_list);
        set_zero(psel_meson_num_list);

        list_n_from_idx_baryon.resize(psel_baryon_num);
        list_n_from_idx_meson.resize(psel_meson_num);
        set_zero(list_n_from_idx_baryon);
        set_zero(list_n_from_idx_meson);

        for (long n = 0; n < psel_meson_num; ++n)
        {
            const long &xg_psel_idx = n;
            const Coordinate xg = psel_meson[xg_psel_idx];
            const int t = xg[3];
            list_n_from_idx_meson[xg_psel_idx] = psel_meson_num_list[t];
            psel_meson_num_list[t] += 1;
        }

        for (long n = 0; n < psel_baryon_num; ++n)
        {
            const long &xg_psel_idx = n;
            const Coordinate xg = psel_baryon[xg_psel_idx];
            const int t = xg[3];
            list_n_from_idx_baryon[xg_psel_idx] = psel_baryon_num_list[t];
            psel_baryon_num_list[t] += 1;
        }

        baryon_pi_two_func_coeff.resize(tsep_num);
        for (int ti = 0; ti < tsep_num; ++ti)
        {
            baryon_pi_two_func_coeff[ti].resize(num_dtxy);
            set_zero(baryon_pi_two_func_coeff[ti]);
        }
        for (int ti = 0; ti < tsep_num; ++ti)
        {
            const int tsep = tsep_list[ti];

            for (int t_meson = 0; t_meson < total_site[3]; t_meson++)
            {

                for (int dt = 0; dt < num_dtxy; ++dt)
                {
                    const int dtxy = dt - num_dtxy / 2;
                    int xt = (tsep >= 0) ? mod(t_meson + dtxy, total_site[3]) : (dtxy >= 0 ? mod(t_meson - 2 * tsep + dtxy, total_site[3]) : mod(t_meson + 2 * tsep + dtxy, total_site[3]));
                    int t_baryon_src, t_baryon_snk;
                    if (dtxy >= 0)
                    {
                        t_baryon_src = mod(t_meson - tsep, total_site[3]);
                        t_baryon_snk = mod(xt + tsep, total_site[3]);
                    }
                    else
                    {
                        t_baryon_src = mod(xt - tsep, total_site[3]);
                        t_baryon_snk = mod(t_meson + tsep, total_site[3]);
                    }
                    if (tsep == 0 && is_baryon_smear == is_meson_smear)
                    {
                        baryon_pi_two_func_coeff[ti][dt] += psel_baryon_num_list[t_baryon_src] * psel_baryon_num_list[t_baryon_snk] * (psel_meson_num_list[t_meson] - 1);
                    }
                    else
                    {
                        baryon_pi_two_func_coeff[ti][dt] += psel_baryon_num_list[t_baryon_src] * psel_baryon_num_list[t_baryon_snk] * psel_meson_num_list[t_meson];
                    }
                }
            }
        }
    }

    inline void baryon_pion_two_point_simple_standard(
        const std::string &job_tag, const int &traj,
        const PointsSelection &psel_baryon,
        const PointsSelection &psel_meson,
        std::vector<int> &psel_baryon_num_list,
        std::vector<int> &psel_meson_num_list,
        std::vector<int> &list_n_from_idx_baryon,
        std::vector<int> &list_n_from_idx_meson,
        std::vector<std::vector<long>> &idx_class_by_time_baryon,
        std::vector<std::vector<long>> &idx_class_by_time_meson,
        const std::vector<std::vector<int>> &Momentum_Targets,
        const Geometry &geo,
        const bool is_baryon_smear, const bool is_meson_smear)
    {
        TIMER_VERBOSE("baryon_pion_two_point_coeffience");

        const Coordinate &total_site = geo.total_site();
        const int num_dtxy = get_baryon_pi_half_sep(total_site[3]);

        const long &psel_meson_num = psel_meson.size();
        const long &psel_baryon_num = psel_baryon.size();

        const unsigned mom_num = Momentum_Targets.size();

        psel_baryon_num_list.resize(total_site[3]);
        psel_meson_num_list.resize(total_site[3]);
        set_zero(psel_baryon_num_list);
        set_zero(psel_meson_num_list);

        list_n_from_idx_baryon.resize(psel_baryon_num);
        list_n_from_idx_meson.resize(psel_meson_num);
        set_zero(list_n_from_idx_baryon);
        set_zero(list_n_from_idx_meson);

        idx_class_by_time_baryon.resize(total_site[3]);
        idx_class_by_time_meson.resize(total_site[3]);

        for (int t_idx = 0; t_idx < total_site[3]; t_idx++)
        {
            std::vector<long>().swap(idx_class_by_time_baryon[t_idx]);
            std::vector<long>().swap(idx_class_by_time_meson[t_idx]);

            // qassert(idx_class_by_time_baryon[t_idx].size() == 0);
            // qassert(idx_class_by_time_meson[t_idx].size() == 0);
        }

        for (long n = 0; n < psel_meson_num; ++n)
        {
            const long &xg_psel_idx = n;
            const Coordinate xg = psel_meson[xg_psel_idx];
            const int t = xg[3];
            list_n_from_idx_meson[xg_psel_idx] = psel_meson_num_list[t];
            psel_meson_num_list[t] += 1;
            idx_class_by_time_meson[t].push_back(xg_psel_idx);
        }

        for (long n = 0; n < psel_baryon_num; ++n)
        {
            const long &xg_psel_idx = n;
            const Coordinate xg = psel_baryon[xg_psel_idx];
            const int t = xg[3];
            list_n_from_idx_baryon[xg_psel_idx] = psel_baryon_num_list[t];
            psel_baryon_num_list[t] += 1;
            idx_class_by_time_baryon[t].push_back(xg_psel_idx);
        }
    }

    inline void baryon_pion_two_point_arbitrary_coeffience(
        const std::string &job_tag, const int &traj,
        std::vector<std::vector<std::vector<std::vector<std::vector<ComplexD>>>>> &baryon_pi_two_func_coeff,
        // std::vector<std::vector<std::vector<std::vector<std::vector<long>>>>> &baryon_pi_two_func_coeff,
        std::vector<int> &psel_baryon_num_list,
        std::vector<int> &psel_meson_num_list,
        std::vector<int> &list_n_from_idx_baryon,
        std::vector<int> &list_n_from_idx_meson,
        std::vector<std::vector<long>> &idx_class_by_time_baryon,
        std::vector<std::vector<long>> &idx_class_by_time_meson,
        const std::vector<std::vector<int>> &Momentum_Targets,
        const Geometry &geo, const PointsSelection &psel_baryon, const PointsSelection &psel_meson,
        const std::vector<int> &tsep_src_list, const std::vector<int> &tsep_snk_list, const bool is_baryon_smear, const bool is_meson_smear)
    {
        TIMER_VERBOSE("baryon_pion_two_point_coeffience");

        // qassert(tsep_src_list.size() == tsep_snk_list.size());
        const int tsep_num = (int)tsep_src_list.size();

        const Coordinate &total_site = geo.total_site();
        const int num_dtxy = get_baryon_pi_half_sep(total_site[3]);

        const long &psel_meson_num = psel_meson.size();
        const long &psel_baryon_num = psel_baryon.size();

        const unsigned mom_num = Momentum_Targets.size();

        psel_baryon_num_list.resize(total_site[3]);
        psel_meson_num_list.resize(total_site[3]);
        set_zero(psel_baryon_num_list);
        set_zero(psel_meson_num_list);

        list_n_from_idx_baryon.resize(psel_baryon_num);
        list_n_from_idx_meson.resize(psel_meson_num);
        set_zero(list_n_from_idx_baryon);
        set_zero(list_n_from_idx_meson);

        idx_class_by_time_baryon.resize(total_site[3]);
        idx_class_by_time_meson.resize(total_site[3]);

        for (int t_idx = 0; t_idx < total_site[3]; t_idx++)
        {
            std::vector<long>().swap(idx_class_by_time_baryon[t_idx]);
            std::vector<long>().swap(idx_class_by_time_meson[t_idx]);

            // qassert(idx_class_by_time_baryon[t_idx].size() == 0);
            // qassert(idx_class_by_time_meson[t_idx].size() == 0);
        }

        for (long n = 0; n < psel_meson_num; ++n)
        {
            const long &xg_psel_idx = n;
            const Coordinate xg = psel_meson[xg_psel_idx];
            const int t = xg[3];
            list_n_from_idx_meson[xg_psel_idx] = psel_meson_num_list[t];
            psel_meson_num_list[t] += 1;
            idx_class_by_time_meson[t].push_back(xg_psel_idx);
        }

        for (long n = 0; n < psel_baryon_num; ++n)
        {
            const long &xg_psel_idx = n;
            const Coordinate xg = psel_baryon[xg_psel_idx];
            const int t = xg[3];
            list_n_from_idx_baryon[xg_psel_idx] = psel_baryon_num_list[t];
            psel_baryon_num_list[t] += 1;
            idx_class_by_time_baryon[t].push_back(xg_psel_idx);
        }

        baryon_pi_two_func_coeff.resize(tsep_num);
        for (int ti = 0; ti < tsep_num; ++ti)
        {
            baryon_pi_two_func_coeff[ti].resize(mom_num);
            for (unsigned int mom_src_id = 0; mom_src_id < mom_num; ++mom_src_id)
            {
                baryon_pi_two_func_coeff[ti][mom_src_id].resize(mom_num);
                for (unsigned int mom_snk_id = 0; mom_snk_id < mom_num; ++mom_snk_id)
                {
                    baryon_pi_two_func_coeff[ti][mom_src_id][mom_snk_id].resize(num_dtxy);
                    for (int dt = 0; dt < num_dtxy; ++dt)
                    {
                        baryon_pi_two_func_coeff[ti][mom_src_id][mom_snk_id][dt].resize(total_site[0]);
                        set_zero(baryon_pi_two_func_coeff[ti][mom_src_id][mom_snk_id][dt]);
                    }
                }
            }
        }

        const int OMP_MAX_THREADS_NUM = omp_get_max_threads();

        // TIMER_VERBOSE("baryon_pion_two_point_coeffience_counting");

        std::vector<std::vector<long>> baryon_pi_two_func_coeff_temp;
        // std::vector<long> baryon_pi_two_func_coeff_temp;

        baryon_pi_two_func_coeff_temp.resize(OMP_MAX_THREADS_NUM);
        for (int thread_id = 0; thread_id < OMP_MAX_THREADS_NUM; thread_id++)
        {
            baryon_pi_two_func_coeff_temp[thread_id].resize(mom_num * mom_num * num_dtxy * total_site[0]);
        }

        baryon_pi_two_func_coeff_temp.resize(mom_num * mom_num * num_dtxy * total_site[0]);

        for (int ti = 0; ti < tsep_num; ++ti)
        {
            const int &tsep_src = tsep_src_list[ti];
            const int &tsep_snk = tsep_snk_list[ti];
            for (int thread_id = 0; thread_id < OMP_MAX_THREADS_NUM; thread_id++)
            {
                set_zero(baryon_pi_two_func_coeff_temp[thread_id]);
            }
// set_zero(baryon_pi_two_func_coeff_temp);
// #pragma omp parallel for reduction(+ : baryon_pi_two_func_coeff_temp)
// shared(is_baryon_smear, psel_baryon, psel_baryon_num, psel_baryon_num_list, list_n_from_idx_baryon, idx_class_by_time_baryon, is_meson_smear, psel_meson, psel_meson_num, psel_meson_num_list, list_n_from_idx_meson, idx_class_by_time_meson)
#pragma omp parallel for collapse(2)
            for (int t_meson = 0; t_meson < total_site[3]; ++t_meson)
            {
                for (int dt = 0; dt < num_dtxy; ++dt)
                {
                    // int dt = 1;
                    const int &dtxy = dt;
                    int xt, t_baryon_src, t_baryon_snk;
                    double anti_period;
                    four_point_time_slice(dtxy, num_dtxy, tsep_src, tsep_snk, total_site, t_meson, t_baryon_src, t_baryon_snk, xt, anti_period);
                    for (int idx_baryon_src = 0; idx_baryon_src < psel_baryon_num_list[t_baryon_src]; idx_baryon_src++)
                    {
                        // int idx_baryon_src = 0;
                        const long &xg_src_psel_idx = idx_class_by_time_baryon[t_baryon_src][idx_baryon_src];
                        const Coordinate &xg_src = psel_baryon[xg_src_psel_idx];
                        // qassert(xg_src[3] == t_baryon_src);
                        qassert(psel_baryon_num_list[t_baryon_src] == (int)idx_class_by_time_baryon[t_baryon_src].size());

                        for (int idx_baryon_snk = 0; idx_baryon_snk < psel_baryon_num_list[t_baryon_snk]; idx_baryon_snk++)
                        {
                            // int idx_baryon_snk = 0;
                            const long &xg_snk_psel_idx = idx_class_by_time_baryon[t_baryon_snk][idx_baryon_snk];
                            const Coordinate &xg_snk = psel_baryon[xg_snk_psel_idx];
                            // qassert(xg_snk[3] == t_baryon_snk);
                            for (int idx_meson_ty = 0; idx_meson_ty < psel_meson_num_list[t_meson]; idx_meson_ty++)
                            {
                                // int idx_meson_ty = 1;
                                const long &xg_y_psel_idx = idx_class_by_time_meson[t_meson][idx_meson_ty];
                                const Coordinate &xg_y = psel_meson[xg_y_psel_idx];
                                // qassert(xg_y[3] == t_meson);

                                if (is_baryon_smear == is_meson_smear && ((xg_y_psel_idx == xg_src_psel_idx) || (xg_y_psel_idx == xg_snk_psel_idx)))
                                {
                                    continue;
                                }

                                for (int idx_x = 0; idx_x < psel_meson_num_list[xt]; idx_x++)
                                {
                                    // int idx_x = 1;
                                    const long &xg_x_psel_idx = idx_class_by_time_meson[xt][idx_x];
                                    const Coordinate &xg_x = psel_meson[xg_x_psel_idx];
                                    // qassert(xg_x[3] == xt);

                                    // if (omp_get_thread_num() == 0)
                                    // {
                                    //     std::cout << "idx_baryon_src: " << idx_baryon_src << std::endl;
                                    //     std::cout << "xg_src_psel_idx: " << xg_src_psel_idx << std::endl;
                                    //     std::cout << "xg_src: (" << xg_src[0] << ", " << xg_src[1] << ", " << xg_src[2] << ", " << xg_src[3] << ")" << std::endl;

                                    //     std::cout << "idx_baryon_snk: " << idx_baryon_snk << std::endl;
                                    //     std::cout << "xg_snk_psel_idx: " << xg_snk_psel_idx << std::endl;
                                    //     std::cout << "xg_snk: (" << xg_snk[0] << ", " << xg_snk[1] << ", " << xg_snk[2] << ", " << xg_snk[3] << ")" << std::endl;

                                    //     std::cout << "idx_meson_ty: " << idx_meson_ty << std::endl;
                                    //     std::cout << "xg_y_psel_idx: " << xg_y_psel_idx << std::endl;
                                    //     std::cout << "xg_y: (" << xg_y[0] << ", " << xg_y[1] << ", " << xg_y[2] << ", " << xg_y[3] << ")" << std::endl;

                                    //     std::cout << "idx_x: " << idx_x << std::endl;
                                    //     std::cout << "xg_x_psel_idx: " << xg_x_psel_idx << std::endl;
                                    //     std::cout << "xg_x: (" << xg_x[0] << ", " << xg_x[1] << ", " << xg_x[2] << ", " << xg_x[3] << ")" << std::endl;
                                    // }
                                    if (is_baryon_smear == is_meson_smear && ((xg_x_psel_idx == xg_src_psel_idx) || (xg_x_psel_idx == xg_snk_psel_idx)))
                                    {
                                        continue;
                                    }

                                    for (unsigned int mom_src_id = 0; mom_src_id < mom_num; ++mom_src_id)
                                    {
                                        for (unsigned int mom_snk_id = 0; mom_snk_id < mom_num; ++mom_snk_id)
                                        {
                                            long phase_temp = space_dot(Momentum_Targets[mom_src_id], xg_src) - space_dot(Momentum_Targets[mom_src_id], xg_y) + space_dot(Momentum_Targets[mom_snk_id], xg_snk) - space_dot(Momentum_Targets[mom_snk_id], xg_x);

                                            // std::cout << "phase_temp: " << phase_temp << std::endl;
                                            int phase = mod(phase_temp, total_site[0]);

                                            // std::cout << "phase: " << phase << "omp_get_thread_num(): " << omp_get_thread_num() << std::endl;
                                            baryon_pi_two_func_coeff_temp[omp_get_thread_num()][((mom_src_id * mom_num + mom_snk_id) * num_dtxy + dt) * total_site[0] + phase] += 1;
                                            // std::cout << "baryon_pi_two_func_coeff_temp[omp_get_thread_num()][((mom_src_id * mom_num + mom_snk_id) * num_dtxy + dt) * total_site[0] + phase] " << baryon_pi_two_func_coeff_temp[omp_get_thread_num()][((mom_src_id * mom_num + mom_snk_id) * num_dtxy + dt) * total_site[0] + phase] << std::endl;
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
// 累加临时数组到最终数组
#pragma omp parallel for collapse(4)
            for (unsigned int mom_src_id = 0; mom_src_id < mom_num; ++mom_src_id)
            {
                for (unsigned int mom_snk_id = 0; mom_snk_id < mom_num; ++mom_snk_id)
                {
                    for (int dt = 0; dt < num_dtxy; ++dt)
                    {
                        for (int phase = 0; phase < total_site[0]; ++phase)
                        {
                            int thread_sum = 0;
                            // 累加所有线程的值
                            for (int thread_id = 0; thread_id < omp_get_max_threads(); ++thread_id)
                            {
                                thread_sum += baryon_pi_two_func_coeff_temp[thread_id][((mom_src_id * mom_num + mom_snk_id) * num_dtxy + dt) * total_site[0] + phase];
                            }
                            // 将累加值写入最终数组
                            // if (thread_sum != 0)
                            // {
                            baryon_pi_two_func_coeff[ti][mom_src_id][mom_snk_id][dt][phase] = exp(ii * 2.0 * PI * (double)phase / (double)total_site[0]) * (1 / (double)thread_sum);
                            // baryon_pi_two_func_coeff[ti][mom_src_id][mom_snk_id][dt][phase] = thread_sum;
                            // }
                            // else
                            // {
                            //     baryon_pi_two_func_coeff[ti][mom_src_id][mom_snk_id][dt][phase] = 0;
                            // }
                        }
                    }
                }
            }
        }
    }

    inline void baryon_pion_two_point_arbitrary_coeffience_fsel(
        const std::string &job_tag, const int &traj,
        std::vector<std::vector<std::vector<std::vector<std::vector<ComplexD>>>>> &baryon_pi_two_func_coeff,
        // std::vector<std::vector<std::vector<std::vector<std::vector<long>>>>> &baryon_pi_two_func_coeff,
        std::vector<int> &psel_baryon_num_list,
        std::vector<int> &psel_meson_num_list,
        std::vector<int> &list_n_from_idx_baryon,
        std::vector<int> &list_n_from_idx_meson,
        std::vector<std::vector<long>> &idx_class_by_time_baryon,
        std::vector<std::vector<long>> &idx_class_by_time_meson,
        std::vector<int> &fsel_num_list,
        std::vector<int> &list_n_from_idx_fsel,
        std::vector<std::vector<long>> &idx_class_by_time_fsel,
        const std::vector<std::vector<int>> &Momentum_Targets_src, const std::vector<std::vector<int>> &Momentum_Targets_snk,
        const FieldSelection &fsel, const Geometry &geo,
        const PointsSelection &psel_baryon, const PointsSelection &psel_meson,
        const std::vector<int> &tsep_src_list, const std::vector<int> &tsep_snk_list, const bool is_baryon_smear, const bool is_meson_smear)
    {
        TIMER_VERBOSE("baryon_pion_two_point_coeffience");

        // qassert(tsep_src_list.size() == tsep_snk_list.size());
        const int tsep_num = (int)tsep_src_list.size();

        const Coordinate &total_site = geo.total_site();
        const int num_dtxy = get_baryon_pi_half_sep(total_site[3]);

        const long &psel_meson_num = psel_meson.size();
        const long &psel_baryon_num = psel_baryon.size();

        const int &mom_num_src = Momentum_Targets_src.size();
        const int &mom_num_snk = Momentum_Targets_snk.size();

        psel_baryon_num_list.resize(total_site[3]);
        psel_meson_num_list.resize(total_site[3]);
        fsel_num_list.resize(total_site[3]);
        set_zero(psel_baryon_num_list);
        set_zero(psel_meson_num_list);
        set_zero(fsel_num_list);

        list_n_from_idx_baryon.resize(psel_baryon_num);
        list_n_from_idx_meson.resize(psel_meson_num);
        list_n_from_idx_fsel.resize(fsel.n_elems);
        set_zero(list_n_from_idx_baryon);
        set_zero(list_n_from_idx_meson);
        set_zero(list_n_from_idx_fsel);

        idx_class_by_time_baryon.resize(total_site[3]);
        idx_class_by_time_meson.resize(total_site[3]);
        idx_class_by_time_fsel.resize(total_site[3]);

        for (int t_idx = 0; t_idx < total_site[3]; t_idx++)
        {
            std::vector<long>().swap(idx_class_by_time_baryon[t_idx]);
            std::vector<long>().swap(idx_class_by_time_meson[t_idx]);
            std::vector<long>().swap(idx_class_by_time_fsel[t_idx]);
        }

        for (long n = 0; n < psel_meson_num; ++n)
        {
            const long &xg_psel_idx = n;
            const Coordinate xg = psel_meson[xg_psel_idx];
            const int t = xg[3];
            list_n_from_idx_meson[xg_psel_idx] = psel_meson_num_list[t];
            psel_meson_num_list[t] += 1;
            idx_class_by_time_meson[t].push_back(xg_psel_idx);
        }

        for (long n = 0; n < psel_baryon_num; ++n)
        {
            const long &xg_psel_idx = n;
            const Coordinate xg = psel_baryon[xg_psel_idx];
            const int t = xg[3];
            list_n_from_idx_baryon[xg_psel_idx] = psel_baryon_num_list[t];
            psel_baryon_num_list[t] += 1;
            idx_class_by_time_baryon[t].push_back(xg_psel_idx);
        }

        for (int n = 0; n < fsel.n_elems; n++)
        {
            long xg_fsel_id = n;
            const long index = fsel.indices[n];
            const Coordinate xg_l = geo.coordinate_from_index(index);
            const Coordinate xg_g = geo.coordinate_g_from_l(xg_l);
            const Coordinate &xg = xg_g;
            const int &t = xg[3];

            list_n_from_idx_fsel[xg_fsel_id] = fsel_num_list[t];
            fsel_num_list[t] += 1;
            idx_class_by_time_fsel[t].push_back(xg_fsel_id);
        }

        baryon_pi_two_func_coeff.resize(tsep_num);
        for (int ti = 0; ti < tsep_num; ++ti)
        {
            baryon_pi_two_func_coeff[ti].resize(mom_num_src);
            for (int mom_src_id = 0; mom_src_id < mom_num_src; ++mom_src_id)
            {
                baryon_pi_two_func_coeff[ti][mom_src_id].resize(mom_num_snk);
                for (int mom_snk_id = 0; mom_snk_id < mom_num_snk; ++mom_snk_id)
                {
                    baryon_pi_two_func_coeff[ti][mom_src_id][mom_snk_id].resize(num_dtxy);
                    for (int dt = 0; dt < num_dtxy; ++dt)
                    {
                        baryon_pi_two_func_coeff[ti][mom_src_id][mom_snk_id][dt].resize(total_site[0]);
                        set_zero(baryon_pi_two_func_coeff[ti][mom_src_id][mom_snk_id][dt]);
                    }
                }
            }
        }

        const int OMP_MAX_THREADS_NUM = omp_get_max_threads();

        // TIMER_VERBOSE("baryon_pion_two_point_coeffience_counting");

        std::vector<std::vector<long>> baryon_pi_two_func_coeff_temp;
        // std::vector<long> baryon_pi_two_func_coeff_temp;

        baryon_pi_two_func_coeff_temp.resize(OMP_MAX_THREADS_NUM);
        for (int thread_id = 0; thread_id < OMP_MAX_THREADS_NUM; thread_id++)
        {
            baryon_pi_two_func_coeff_temp[thread_id].resize(mom_num_src * mom_num_snk * num_dtxy * total_site[0]);
        }

        baryon_pi_two_func_coeff_temp.resize(mom_num_src * mom_num_snk * num_dtxy * total_site[0]);

        for (int ti = 0; ti < tsep_num; ++ti)
        {
            const int &tsep_src = tsep_src_list[ti];
            const int &tsep_snk = tsep_snk_list[ti];
            for (int thread_id = 0; thread_id < OMP_MAX_THREADS_NUM; thread_id++)
            {
                set_zero(baryon_pi_two_func_coeff_temp[thread_id]);
            }
// set_zero(baryon_pi_two_func_coeff_temp);
// #pragma omp parallel for reduction(+ : baryon_pi_two_func_coeff_temp)
// shared(is_baryon_smear, psel_baryon, psel_baryon_num, psel_baryon_num_list, list_n_from_idx_baryon, idx_class_by_time_baryon, is_meson_smear, psel_meson, psel_meson_num, psel_meson_num_list, list_n_from_idx_meson, idx_class_by_time_meson)
#pragma omp parallel for collapse(2)
            for (int t_meson = 0; t_meson < total_site[3]; ++t_meson)
            {
                for (int dt = 0; dt < num_dtxy; ++dt)
                {
                    // int dt = 1;
                    const int &dtxy = dt;
                    int xt, t_baryon_src, t_baryon_snk;
                    double anti_period;
                    four_point_time_slice(dtxy, num_dtxy, tsep_src, tsep_snk, total_site, t_meson, t_baryon_src, t_baryon_snk, xt, anti_period);
                    for (int idx_baryon_src = 0; idx_baryon_src < psel_baryon_num_list[t_baryon_src]; idx_baryon_src++)
                    {
                        // int idx_baryon_src = 0;
                        const long &xg_src_psel_idx = idx_class_by_time_baryon[t_baryon_src][idx_baryon_src];
                        const Coordinate &xg_src = psel_baryon[xg_src_psel_idx];
                        // qassert(xg_src[3] == t_baryon_src);
                        qassert(psel_baryon_num_list[t_baryon_src] == (int)idx_class_by_time_baryon[t_baryon_src].size());

                        for (int idx_baryon_snk = 0; idx_baryon_snk < psel_baryon_num_list[t_baryon_snk]; idx_baryon_snk++)
                        {
                            // int idx_baryon_snk = 0;
                            const long &xg_snk_psel_idx = idx_class_by_time_baryon[t_baryon_snk][idx_baryon_snk];
                            const Coordinate &xg_snk = psel_baryon[xg_snk_psel_idx];
                            // qassert(xg_snk[3] == t_baryon_snk);
                            for (int idx_meson_ty = 0; idx_meson_ty < psel_meson_num_list[t_meson]; idx_meson_ty++)
                            {
                                // int idx_meson_ty = 1;
                                const long &xg_y_psel_idx = idx_class_by_time_meson[t_meson][idx_meson_ty];
                                const Coordinate &xg_y = psel_meson[xg_y_psel_idx];
                                // qassert(xg_y[3] == t_meson);

                                if (is_baryon_smear == is_meson_smear && ((xg_y_psel_idx == xg_src_psel_idx) || (xg_y_psel_idx == xg_snk_psel_idx)))
                                {
                                    continue;
                                }

                                for (int idx_x = 0; idx_x < fsel_num_list[xt]; idx_x++)
                                {
                                    long xg_x_fsel_idx = idx_class_by_time_fsel[xt][idx_x];
                                    const long index = fsel.indices[xg_x_fsel_idx];
                                    const Coordinate xg_x_l = geo.coordinate_from_index(index);
                                    const Coordinate xg_x_g = geo.coordinate_g_from_l(xg_x_l);
                                    const Coordinate &xg_x = xg_x_g;
                                    const int &t_x = xg_x[3];
                                    qassert(t_x == xt);
                                    qassert(idx_x == list_n_from_idx_fsel[xg_x_fsel_idx]);
                                    for (int mom_src_id = 0; mom_src_id < mom_num_src; ++mom_src_id)
                                    {
                                        for (int mom_snk_id = 0; mom_snk_id < mom_num_snk; ++mom_snk_id)
                                        {
                                            long phase_temp = space_dot(Momentum_Targets_src[mom_src_id], xg_src) - space_dot(Momentum_Targets_src[mom_src_id], xg_y) + space_dot(Momentum_Targets_snk[mom_snk_id], xg_snk);
                                            //  - space_dot(Momentum_Targets_snk[mom_snk_id], xg_x);

                                            int phase = mod(phase_temp, total_site[0]);

                                            // std::cout << "phase: " << phase << "omp_get_thread_num(): " << omp_get_thread_num() << std::endl;
                                            baryon_pi_two_func_coeff_temp[omp_get_thread_num()][((mom_src_id * mom_num_snk + mom_snk_id) * num_dtxy + dt) * total_site[0] + phase] += 1;
                                            // std::cout << "baryon_pi_two_func_coeff_temp[omp_get_thread_num()][((mom_src_id * mom_num + mom_snk_id) * num_dtxy + dt) * total_site[0] + phase] " << baryon_pi_two_func_coeff_temp[omp_get_thread_num()][((mom_src_id * mom_num + mom_snk_id) * num_dtxy + dt) * total_site[0] + phase] << std::endl;
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }

            std::vector<Long> baryon_pi_two_func_phase_iteration;
            baryon_pi_two_func_phase_iteration.resize(mom_num_src * mom_num_snk * num_dtxy * total_site[0]);
// 累加临时数组到最终数组
#pragma omp parallel for collapse(4)
            for (int mom_src_id = 0; mom_src_id < mom_num_src; ++mom_src_id)
            {
                for (int mom_snk_id = 0; mom_snk_id < mom_num_snk; ++mom_snk_id)
                {
                    for (int dt = 0; dt < num_dtxy; ++dt)
                    {
                        for (int phase = 0; phase < total_site[0]; ++phase)
                        {
                            int thread_sum = 0;
                            // 累加所有线程的值
                            for (int thread_id = 0; thread_id < omp_get_max_threads(); ++thread_id)
                            {
                                thread_sum += baryon_pi_two_func_coeff_temp[thread_id][((mom_src_id * mom_num_snk + mom_snk_id) * num_dtxy + dt) * total_site[0] + phase];
                            }
                            baryon_pi_two_func_phase_iteration[((mom_src_id * mom_num_snk + mom_snk_id) * num_dtxy + dt) * total_site[0] + phase] = thread_sum;
                        }
                    }
                }
            }

            glb_sum(baryon_pi_two_func_phase_iteration);

// 累加临时数组到最终数组
#pragma omp parallel for collapse(4)
            for (int mom_src_id = 0; mom_src_id < mom_num_src; ++mom_src_id)
            {
                for (int mom_snk_id = 0; mom_snk_id < mom_num_snk; ++mom_snk_id)
                {
                    for (int dt = 0; dt < num_dtxy; ++dt)
                    {
                        for (int phase = 0; phase < total_site[0]; ++phase)
                        {
                            baryon_pi_two_func_coeff[ti][mom_src_id][mom_snk_id][dt][phase] = exp(ii * 2.0 * PI * (double)phase / (double)total_site[0]) * ((double)total_site[0] * (double)total_site[1] * (double)total_site[2] / 16 / (double)baryon_pi_two_func_phase_iteration[((mom_src_id * mom_num_snk + mom_snk_id) * num_dtxy + dt) * total_site[0] + phase]);
                        }
                    }
                }
            }
        }
    }

    inline void baryon_pion_two_point_arbitrary_coeffience_mate_fsel(
        const std::string &job_tag, const int &traj,
        std::vector<std::vector<std::vector<std::vector<std::vector<ComplexD>>>>> &baryon_pi_two_func_coeff,
        std::vector<int> &psel_baryon_num_list,
        std::vector<int> &psel_meson_num_list,
        std::vector<int> &list_n_from_idx_baryon,
        std::vector<int> &list_n_from_idx_meson,
        std::vector<std::vector<long>> &idx_class_by_time_baryon,
        std::vector<std::vector<long>> &idx_class_by_time_meson,
        std::vector<int> &fsel_num_list,
        std::vector<int> &list_n_from_idx_fsel,
        std::vector<std::vector<long>> &idx_class_by_time_fsel,
        const std::vector<std::vector<int>> &Momentum_Targets_src,
        const std::vector<std::vector<int>> &Momentum_Targets_curr,
        const std::vector<std::vector<int>> &Momentum_Targets_snk,
        const FieldSelection &fsel, const Geometry &geo,
        const PointsSelection &psel_baryon, const PointsSelection &psel_meson,
        const std::vector<int> &tsep_src_list, const std::vector<int> &tsep_snk_list, const bool is_baryon_smear, const bool is_meson_smear)
    {
        TIMER_VERBOSE("baryon_pion_two_point_coeffience");

        // qassert(tsep_src_list.size() == tsep_snk_list.size());
        const int tsep_num = (int)tsep_src_list.size();

        const Coordinate &total_site = geo.total_site();
        const int num_dtxy = get_baryon_pi_half_sep(total_site[3]);

        const long &psel_meson_num = psel_meson.size();
        const long &psel_baryon_num = psel_baryon.size();

        const int &mom_num_src = Momentum_Targets_src.size();
        const int &mom_num_curr = Momentum_Targets_curr.size();
        const int &mom_num_snk = Momentum_Targets_snk.size();

        psel_baryon_num_list.resize(total_site[3]);
        psel_meson_num_list.resize(total_site[3]);
        fsel_num_list.resize(total_site[3]);
        set_zero(psel_baryon_num_list);
        set_zero(psel_meson_num_list);
        set_zero(fsel_num_list);

        list_n_from_idx_baryon.resize(psel_baryon_num);
        list_n_from_idx_meson.resize(psel_meson_num);
        list_n_from_idx_fsel.resize(fsel.n_elems);
        set_zero(list_n_from_idx_baryon);
        set_zero(list_n_from_idx_meson);
        set_zero(list_n_from_idx_fsel);

        idx_class_by_time_baryon.resize(total_site[3]);
        idx_class_by_time_meson.resize(total_site[3]);
        idx_class_by_time_fsel.resize(total_site[3]);

        for (int t_idx = 0; t_idx < total_site[3]; t_idx++)
        {
            std::vector<long>().swap(idx_class_by_time_baryon[t_idx]);
            std::vector<long>().swap(idx_class_by_time_meson[t_idx]);
            std::vector<long>().swap(idx_class_by_time_fsel[t_idx]);
        }

        for (long n = 0; n < psel_meson_num; ++n)
        {
            const long &xg_psel_idx = n;
            const Coordinate xg = psel_meson[xg_psel_idx];
            const int t = xg[3];
            list_n_from_idx_meson[xg_psel_idx] = psel_meson_num_list[t];
            psel_meson_num_list[t] += 1;
            idx_class_by_time_meson[t].push_back(xg_psel_idx);
        }

        for (long n = 0; n < psel_baryon_num; ++n)
        {
            const long &xg_psel_idx = n;
            const Coordinate xg = psel_baryon[xg_psel_idx];
            const int t = xg[3];
            list_n_from_idx_baryon[xg_psel_idx] = psel_baryon_num_list[t];
            psel_baryon_num_list[t] += 1;
            idx_class_by_time_baryon[t].push_back(xg_psel_idx);
        }

        for (int n = 0; n < fsel.n_elems; n++)
        {
            long xg_fsel_id = n;
            const long index = fsel.indices[n];
            const Coordinate xg_l = geo.coordinate_from_index(index);
            const Coordinate xg_g = geo.coordinate_g_from_l(xg_l);
            const Coordinate &xg = xg_g;
            const int &t = xg[3];

            list_n_from_idx_fsel[xg_fsel_id] = fsel_num_list[t];
            fsel_num_list[t] += 1;
            idx_class_by_time_fsel[t].push_back(xg_fsel_id);
        }

        baryon_pi_two_func_coeff.resize(tsep_num);
        for (int ti = 0; ti < tsep_num; ++ti)
        {
            baryon_pi_two_func_coeff[ti].resize(mom_num_src);
            for (int mom_src_id = 0; mom_src_id < mom_num_src; ++mom_src_id)
            {
                baryon_pi_two_func_coeff[ti][mom_src_id].resize(mom_num_snk);
                for (int mom_snk_id = 0; mom_snk_id < mom_num_snk; ++mom_snk_id)
                {
                    baryon_pi_two_func_coeff[ti][mom_src_id][mom_snk_id].resize(num_dtxy);
                    for (int dt = 0; dt < num_dtxy; ++dt)
                    {
                        baryon_pi_two_func_coeff[ti][mom_src_id][mom_snk_id][dt].resize(total_site[0]);
                        set_zero(baryon_pi_two_func_coeff[ti][mom_src_id][mom_snk_id][dt]);
                    }
                }
            }
        }

        const int OMP_MAX_THREADS_NUM = omp_get_max_threads();

        // TIMER_VERBOSE("baryon_pion_two_point_coeffience_counting");

        std::vector<std::vector<long>> baryon_pi_two_func_coeff_temp;
        // std::vector<long> baryon_pi_two_func_coeff_temp;

        baryon_pi_two_func_coeff_temp.resize(OMP_MAX_THREADS_NUM);
        for (int thread_id = 0; thread_id < OMP_MAX_THREADS_NUM; thread_id++)
        {
            baryon_pi_two_func_coeff_temp[thread_id].resize(mom_num_src * mom_num_snk * num_dtxy * total_site[0]);
        }

        for (int ti = 0; ti < tsep_num; ++ti)
        {
            const int &tsep_src = tsep_src_list[ti];
            const int &tsep_snk = tsep_snk_list[ti];
            for (int thread_id = 0; thread_id < OMP_MAX_THREADS_NUM; thread_id++)
            {
                set_zero(baryon_pi_two_func_coeff_temp[thread_id]);
            }
// set_zero(baryon_pi_two_func_coeff_temp);
// #pragma omp parallel for reduction(+ : baryon_pi_two_func_coeff_temp)
// shared(is_baryon_smear, psel_baryon, psel_baryon_num, psel_baryon_num_list, list_n_from_idx_baryon, idx_class_by_time_baryon, is_meson_smear, psel_meson, psel_meson_num, psel_meson_num_list, list_n_from_idx_meson, idx_class_by_time_meson)
#pragma omp parallel for collapse(2)
            for (int t_meson = 0; t_meson < total_site[3]; ++t_meson)
            {
                for (int dt = 0; dt < num_dtxy; ++dt)
                {
                    // int dt = 1;
                    const int &dtxy = dt;
                    int xt, t_baryon_src, t_baryon_snk;
                    double anti_period;
                    four_point_time_slice(dtxy, num_dtxy, tsep_src, tsep_snk, total_site, t_meson, t_baryon_src, t_baryon_snk, xt, anti_period);
                    for (int idx_baryon_src = 0; idx_baryon_src < psel_baryon_num_list[t_baryon_src]; idx_baryon_src++)
                    {
                        // int idx_baryon_src = 0;
                        const long &xg_src_psel_idx = idx_class_by_time_baryon[t_baryon_src][idx_baryon_src];
                        const Coordinate &xg_src = psel_baryon[xg_src_psel_idx];
                        // qassert(xg_src[3] == t_baryon_src);
                        qassert(psel_baryon_num_list[t_baryon_src] == (int)idx_class_by_time_baryon[t_baryon_src].size());

                        for (int idx_baryon_snk = 0; idx_baryon_snk < psel_baryon_num_list[t_baryon_snk]; idx_baryon_snk++)
                        {
                            // int idx_baryon_snk = 0;
                            const long &xg_snk_psel_idx = idx_class_by_time_baryon[t_baryon_snk][idx_baryon_snk];
                            const Coordinate &xg_snk = psel_baryon[xg_snk_psel_idx];
                            // qassert(xg_snk[3] == t_baryon_snk);
                            for (int idx_meson_ty = 0; idx_meson_ty < psel_meson_num_list[t_meson]; idx_meson_ty++)
                            {
                                // int idx_meson_ty = 1;
                                const long &xg_y_psel_idx = idx_class_by_time_meson[t_meson][idx_meson_ty];
                                const Coordinate &xg_y = psel_meson[xg_y_psel_idx];
                                // qassert(xg_y[3] == t_meson);

                                if (is_baryon_smear == is_meson_smear && ((xg_y_psel_idx == xg_src_psel_idx) || (xg_y_psel_idx == xg_snk_psel_idx)))
                                {
                                    continue;
                                }

                                for (int mom_src_id = 0; mom_src_id < mom_num_src; ++mom_src_id)
                                {
                                    for (int mom_snk_id = 0; mom_snk_id < mom_num_snk; ++mom_snk_id)
                                    {
                                        long phase_temp = space_dot(Momentum_Targets_src[mom_src_id], xg_src) - space_dot(Momentum_Targets_src[mom_src_id], xg_y) + space_dot(Momentum_Targets_snk[mom_snk_id], xg_snk);

                                        int phase = mod(phase_temp, total_site[0]);

                                        // std::cout << "phase: " << phase << "omp_get_thread_num(): " << omp_get_thread_num() << std::endl;
                                        baryon_pi_two_func_coeff_temp[omp_get_thread_num()][((mom_src_id * mom_num_snk + mom_snk_id) * num_dtxy + dt) * total_site[0] + phase] += 1;
                                        // std::cout << "baryon_pi_two_func_coeff_temp[omp_get_thread_num()][((mom_src_id * mom_num + mom_snk_id) * num_dtxy + dt) * total_site[0] + phase] " << baryon_pi_two_func_coeff_temp[omp_get_thread_num()][((mom_src_id * mom_num + mom_snk_id) * num_dtxy + dt) * total_site[0] + phase] << std::endl;
                                    }
                                }
                            }
                        }
                    }
                }
            }
            std::vector<Long> baryon_pi_two_func_phase_iteration;
            baryon_pi_two_func_phase_iteration.resize(mom_num_src * mom_num_snk * num_dtxy * total_site[0]);
// 累加临时数组到最终数组
#pragma omp parallel for collapse(4)
            for (int mom_src_id = 0; mom_src_id < mom_num_src; ++mom_src_id)
            {
                for (int mom_snk_id = 0; mom_snk_id < mom_num_snk; ++mom_snk_id)
                {
                    for (int dt = 0; dt < num_dtxy; ++dt)
                    {
                        for (int phase = 0; phase < total_site[0]; ++phase)
                        {
                            int thread_sum = 0;
                            // 累加所有线程的值
                            for (int thread_id = 0; thread_id < omp_get_max_threads(); ++thread_id)
                            {
                                thread_sum += baryon_pi_two_func_coeff_temp[thread_id][((mom_src_id * mom_num_snk + mom_snk_id) * num_dtxy + dt) * total_site[0] + phase];
                            }
                            baryon_pi_two_func_phase_iteration[((mom_src_id * mom_num_snk + mom_snk_id) * num_dtxy + dt) * total_site[0] + phase] = thread_sum;
                        }
                    }
                }
            }

            // glb_sum(baryon_pi_two_func_phase_iteration);

// 累加临时数组到最终数组
#pragma omp parallel for collapse(4)
            for (int mom_src_id = 0; mom_src_id < mom_num_src; ++mom_src_id)
            {
                for (int mom_snk_id = 0; mom_snk_id < mom_num_snk; ++mom_snk_id)
                {
                    for (int dt = 0; dt < num_dtxy; ++dt)
                    {
                        for (int phase = 0; phase < total_site[0]; ++phase)
                        {
                            baryon_pi_two_func_coeff[ti][mom_src_id][mom_snk_id][dt][phase] = exp(ii * 2.0 * PI * (double)phase / (double)total_site[0]) / (double)baryon_pi_two_func_phase_iteration[((mom_src_id * mom_num_snk + mom_snk_id) * num_dtxy + dt) * total_site[0] + phase];
                        }
                    }
                }
            }
        }
    }

    inline void baryon_pion_two_point_arbitrary_coeffience_mate_magic_fsel(
        const std::string &job_tag, const int &traj,
        std::vector<std::vector<std::vector<std::vector<std::vector<std::vector<ComplexD>>>>>> &baryon_pi_two_func_coeff,
        std::vector<int> &psel_baryon_num_list,
        std::vector<int> &psel_meson_num_list,
        std::vector<int> &list_n_from_idx_baryon,
        std::vector<int> &list_n_from_idx_meson,
        std::vector<std::vector<long>> &idx_class_by_time_baryon,
        std::vector<std::vector<long>> &idx_class_by_time_meson,
        std::vector<int> &fsel_num_list,
        std::vector<int> &list_n_from_idx_fsel,
        std::vector<std::vector<long>> &idx_class_by_time_fsel,
        const std::vector<std::vector<int>> &Momentum_Targets_src,
        const std::vector<std::vector<int>> &Momentum_Targets_curr,
        const FieldSelection &fsel, const Geometry &geo,
        const PointsSelection &psel_baryon, const PointsSelection &psel_meson,
        const std::vector<int> &tsep_src_list, const std::vector<int> &tsep_snk_list, const bool is_baryon_smear, const bool is_meson_smear)
    {
        TIMER_VERBOSE("baryon_pion_two_point_coeffience");

        // qassert(tsep_src_list.size() == tsep_snk_list.size());
        const int tsep_num = (int)tsep_src_list.size();

        const Coordinate &total_site = geo.total_site();
        const int num_dtxy = get_baryon_pi_half_sep(total_site[3]);

        const long &psel_meson_num = psel_meson.size();
        const long &psel_baryon_num = psel_baryon.size();

        const int &mom_num_src = Momentum_Targets_src.size();
        const int &mom_num_curr = Momentum_Targets_curr.size();

        psel_baryon_num_list.resize(total_site[3]);
        psel_meson_num_list.resize(total_site[3]);
        fsel_num_list.resize(total_site[3]);
        set_zero(psel_baryon_num_list);
        set_zero(psel_meson_num_list);
        set_zero(fsel_num_list);

        list_n_from_idx_baryon.resize(psel_baryon_num);
        list_n_from_idx_meson.resize(psel_meson_num);
        list_n_from_idx_fsel.resize(fsel.n_elems);
        set_zero(list_n_from_idx_baryon);
        set_zero(list_n_from_idx_meson);
        set_zero(list_n_from_idx_fsel);

        idx_class_by_time_baryon.resize(total_site[3]);
        idx_class_by_time_meson.resize(total_site[3]);
        idx_class_by_time_fsel.resize(total_site[3]);

        for (int t_idx = 0; t_idx < total_site[3]; t_idx++)
        {
            std::vector<long>().swap(idx_class_by_time_baryon[t_idx]);
            std::vector<long>().swap(idx_class_by_time_meson[t_idx]);
            std::vector<long>().swap(idx_class_by_time_fsel[t_idx]);
        }

        for (long n = 0; n < psel_meson_num; ++n)
        {
            const long &xg_psel_idx = n;
            const Coordinate xg = psel_meson[xg_psel_idx];
            const int t = xg[3];
            list_n_from_idx_meson[xg_psel_idx] = psel_meson_num_list[t];
            psel_meson_num_list[t] += 1;
            idx_class_by_time_meson[t].push_back(xg_psel_idx);
        }

        for (long n = 0; n < psel_baryon_num; ++n)
        {
            const long &xg_psel_idx = n;
            const Coordinate xg = psel_baryon[xg_psel_idx];
            const int t = xg[3];
            list_n_from_idx_baryon[xg_psel_idx] = psel_baryon_num_list[t];
            psel_baryon_num_list[t] += 1;
            idx_class_by_time_baryon[t].push_back(xg_psel_idx);
        }

        for (int n = 0; n < fsel.n_elems; n++)
        {
            long xg_fsel_id = n;
            const long index = fsel.indices[n];
            const Coordinate xg_l = geo.coordinate_from_index(index);
            const Coordinate xg_g = geo.coordinate_g_from_l(xg_l);
            const Coordinate &xg = xg_g;
            const int &t = xg[3];

            list_n_from_idx_fsel[xg_fsel_id] = fsel_num_list[t];
            fsel_num_list[t] += 1;
            idx_class_by_time_fsel[t].push_back(xg_fsel_id);
        }

        baryon_pi_two_func_coeff.resize(tsep_num);
        for (int ti = 0; ti < tsep_num; ++ti)
        {
            const int &tsep_src = tsep_src_list[ti];
            const int &tsep_snk = tsep_snk_list[ti];

            baryon_pi_two_func_coeff[ti].resize(num_dtxy);
            for (int t_src_snk = 0; t_src_snk < num_dtxy; ++t_src_snk)
            {
                baryon_pi_two_func_coeff[ti][t_src_snk].resize(mom_num_src);
                for (int mom_src_id = 0; mom_src_id < mom_num_src; ++mom_src_id)
                {
                    baryon_pi_two_func_coeff[ti][t_src_snk][mom_src_id].resize(mom_num_curr);
                    for (int mom_curr_id = 0; mom_curr_id < mom_num_curr; ++mom_curr_id)
                    {
                        baryon_pi_two_func_coeff[ti][t_src_snk][mom_src_id][mom_curr_id].resize(t_src_snk);
                        for (int dt = 0; dt < t_src_snk; ++dt)
                        {
                            baryon_pi_two_func_coeff[ti][t_src_snk][mom_src_id][mom_curr_id][dt].resize(total_site[0]);
                            set_zero(baryon_pi_two_func_coeff[ti][t_src_snk][mom_src_id][mom_curr_id][dt]);
                        }
                    }
                }
            }
        }

        const int OMP_MAX_THREADS_NUM = omp_get_max_threads();

        // TIMER_VERBOSE("baryon_pion_two_point_coeffience_counting");

        std::vector<std::vector<long>> baryon_pi_two_func_coeff_temp;
        // std::vector<long> baryon_pi_two_func_coeff_temp;

        baryon_pi_two_func_coeff_temp.resize(OMP_MAX_THREADS_NUM);

        for (int ti = 0; ti < tsep_num; ++ti)
        {
            const int &tsep_src = tsep_src_list[ti];
            const int &tsep_snk = tsep_snk_list[ti];

            for (int t_src_snk = 0; t_src_snk < num_dtxy; ++t_src_snk)
            {
                for (int thread_id = 0; thread_id < OMP_MAX_THREADS_NUM; thread_id++)
                {
                    baryon_pi_two_func_coeff_temp[thread_id].resize(mom_num_src * mom_num_curr * t_src_snk * total_site[0]);
                    set_zero(baryon_pi_two_func_coeff_temp[thread_id]);
                }
#pragma omp parallel for collapse(2)
                for (int t_baryon_src = 0; t_baryon_src < total_site[3]; ++t_baryon_src)
                {
                    for (int dt = 0; dt < t_src_snk; ++dt)
                    {
                        const int &dtxy = dt;
                        int xt, t_meson, t_baryon_snk;
                        double anti_period;
                        three_point_time_slice(dtxy, t_src_snk, tsep_src, tsep_snk, total_site, t_meson, t_baryon_src, t_baryon_snk, xt, anti_period);

                        for (int idx_baryon_src = 0; idx_baryon_src < psel_baryon_num_list[t_baryon_src]; idx_baryon_src++)
                        {
                            const long &xg_src_psel_idx = idx_class_by_time_baryon[t_baryon_src][idx_baryon_src];
                            const Coordinate &xg_src = psel_baryon[xg_src_psel_idx];
                            qassert(xg_src[3] == t_baryon_src);
                            qassert(psel_baryon_num_list[t_baryon_src] == (int)idx_class_by_time_baryon[t_baryon_src].size());

                            for (int idx_baryon_snk = 0; idx_baryon_snk < psel_baryon_num_list[t_baryon_snk]; idx_baryon_snk++)
                            {
                                const long &xg_snk_psel_idx = idx_class_by_time_baryon[t_baryon_snk][idx_baryon_snk];
                                const Coordinate &xg_snk = psel_baryon[xg_snk_psel_idx];
                                qassert(xg_snk[3] == t_baryon_snk);
                                for (int idx_meson_ty = 0; idx_meson_ty < psel_meson_num_list[t_meson]; idx_meson_ty++)
                                {
                                    const long &xg_y_psel_idx = idx_class_by_time_meson[t_meson][idx_meson_ty];
                                    const Coordinate &xg_y = psel_meson[xg_y_psel_idx];
                                    qassert(xg_y[3] == t_meson);

                                    if (is_baryon_smear == is_meson_smear && ((xg_y_psel_idx == xg_src_psel_idx) || (xg_y_psel_idx == xg_snk_psel_idx)))
                                    {
                                        continue;
                                    }
                                    for (int mom_src_id = 0; mom_src_id < mom_num_src; ++mom_src_id)
                                    {
                                        for (int mom_curr_id = 0; mom_curr_id < mom_num_curr; ++mom_curr_id)
                                        {
                                            long phase_temp = space_dot(Momentum_Targets_src[mom_src_id], xg_src) - space_dot(Momentum_Targets_src[mom_src_id], xg_y) + space_dot(Momentum_Targets_curr[mom_curr_id], xg_snk);

                                            int phase = mod(phase_temp, total_site[0]);

                                            baryon_pi_two_func_coeff_temp[omp_get_thread_num()][((mom_src_id * mom_num_curr + mom_curr_id) * t_src_snk + dt) * total_site[0] + phase] += 1;
                                        }
                                    }
                                }
                            }
                        }
                    }
                }

                std::vector<Long> baryon_pi_two_func_phase_iteration;
                baryon_pi_two_func_phase_iteration.resize(mom_num_src * mom_num_curr * t_src_snk * total_site[0]);
// 累加临时数组到最终数组
#pragma omp parallel for collapse(4)
                for (int mom_src_id = 0; mom_src_id < mom_num_src; ++mom_src_id)
                {
                    for (int mom_curr_id = 0; mom_curr_id < mom_num_curr; ++mom_curr_id)
                    {
                        for (int dt = 0; dt < t_src_snk; ++dt)
                        {
                            for (int phase = 0; phase < total_site[0]; ++phase)
                            {
                                int thread_sum = 0;
                                // 累加所有线程的值
                                for (int thread_id = 0; thread_id < omp_get_max_threads(); ++thread_id)
                                {
                                    thread_sum += baryon_pi_two_func_coeff_temp[thread_id][((mom_src_id * mom_num_curr + mom_curr_id) * t_src_snk + dt) * total_site[0] + phase];
                                }
                                baryon_pi_two_func_phase_iteration[((mom_src_id * mom_num_curr + mom_curr_id) * t_src_snk + dt) * total_site[0] + phase] = thread_sum;
                            }
                        }
                    }
                }

                // glb_sum(baryon_pi_two_func_phase_iteration);

// 累加临时数组到最终数组
#pragma omp parallel for collapse(4)
                for (int mom_src_id = 0; mom_src_id < mom_num_src; ++mom_src_id)
                {
                    for (int mom_curr_id = 0; mom_curr_id < mom_num_curr; ++mom_curr_id)
                    {
                        for (int dt = 0; dt < t_src_snk; ++dt)
                        {
                            for (int phase = 0; phase < total_site[0]; ++phase)
                            {
                                baryon_pi_two_func_coeff[ti][t_src_snk][mom_src_id][mom_curr_id][dt][phase] = exp(ii * 2.0 * PI * (double)phase / (double)total_site[0]) / (double)baryon_pi_two_func_phase_iteration[((mom_src_id * mom_num_curr + mom_curr_id) * t_src_snk + dt) * total_site[0] + phase];
                            }
                        }
                    }
                }
            }
        }
    }
}
