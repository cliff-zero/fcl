#pragma once

#include "data-load.h"
#include "contract_block.h"
#include "baryon-pion-block-acc.h"
#include "compute-utils.h"

namespace qlat
{
    inline void contract_two_point_function_coefficent_pro(
        std::vector<ComplexD> &coefficent_of_t2pf,
        std::vector<ComplexD> &coefficent_of_t2pf_remain_snk,
        const std::vector<std::vector<int>> &Momentum_Targets,
        const bool &is_baryon_smear,
        const PointsSelection &psel_baryon,
        const Long &psel_baryon_num,
        const std::vector<int> &psel_baryon_num_list,
        const std::vector<int> &list_n_from_idx_baryon,
        const std::vector<std::vector<Long>> &idx_class_by_time_baryon,
        const bool &is_meson_smear,
        const PointsSelection &psel_meson,
        const Long &psel_meson_num,
        const std::vector<int> &psel_meson_num_list,
        const std::vector<int> &list_n_from_idx_meson,
        const std::vector<std::vector<Long>> &idx_class_by_time_meson,
        const Coordinate &total_site)
    {
        TIMER_VERBOSE("contract_two_point_function_coefficent_pro");
        const int &mom_num = Momentum_Targets.size();

        const Long data_length = total_site[3] * total_site[3] * mom_num * mom_num * total_site[0];
        std::vector<Long> coefficent_of_t2pf_omp(data_length * omp_get_max_threads());
        set_zero(coefficent_of_t2pf);

#pragma omp parallel for collapse(2)
        for (Long n_snk = 0; n_snk < psel_baryon_num; n_snk++)
        {
            for (int dt = 0; dt < total_site[3]; ++dt)
            {
                const Long &xg_snk_psel_idx = n_snk;
                const Coordinate &xg_snk = psel_baryon[xg_snk_psel_idx];
                const int &t_snk = xg_snk[3];

                const int &del_t = dt;
                const int t_src = mod(t_snk - del_t, total_site[3]);
                for (int idx_src = 0; idx_src < psel_baryon_num_list[t_src]; idx_src++)
                {
                    const Long &xg_src_psel_idx = idx_class_by_time_baryon[t_src][idx_src];
                    const Coordinate &xg_src = psel_baryon[xg_src_psel_idx];
                    for (int mom_src_id = 0; mom_src_id < mom_num; ++mom_src_id)
                    {
                        for (int mom_snk_id = 0; mom_snk_id < mom_num; ++mom_snk_id)
                        {
                            Long phase_temp = space_dot(Momentum_Targets[mom_src_id], xg_src) - space_dot(Momentum_Targets[mom_snk_id], xg_snk);
                            int phase = mod(phase_temp, total_site[0]);

                            coefficent_of_t2pf_omp[(((((Long)omp_get_thread_num() * total_site[3] + t_src) * total_site[3] + dt) * mom_num + mom_snk_id) * mom_num + mom_src_id) * total_site[0] + phase] += 1.0;
                        }
                    }
                }
            }
        }

        coefficent_of_t2pf_remain_snk.resize(data_length);
        set_zero(coefficent_of_t2pf_remain_snk);
#pragma omp parallel for collapse(4)
        for (int t_src = 0; t_src < total_site[3]; ++t_src)
        {
            for (int dt = 0; dt < total_site[3]; ++dt)
            {
                for (int mom_src_id = 0; mom_src_id < mom_num; ++mom_src_id)
                {
                    for (int mom_snk_id = 0; mom_snk_id < mom_num; ++mom_snk_id)
                    {
                        for (int phase = 0; phase < total_site[0]; ++phase)
                        {
                            for (int thread_id = 1; thread_id < omp_get_max_threads(); thread_id++)
                            {
                                coefficent_of_t2pf_omp[((((Long)t_src * total_site[3] + dt) * mom_num + mom_snk_id) * mom_num + mom_src_id) * total_site[0] + phase] += coefficent_of_t2pf_omp[(((((Long)thread_id * total_site[3] + t_src) * total_site[3] + dt) * mom_num + mom_snk_id) * mom_num + mom_src_id) * total_site[0] + phase];
                            }
                            coefficent_of_t2pf_remain_snk[((((Long)t_src * total_site[3] + dt) * mom_num + mom_snk_id) * mom_num + mom_src_id) * total_site[0] + phase] = exp(ii * 2.0 * PI * (double)phase / (double)total_site[0]) / (double)coefficent_of_t2pf_omp[((((Long)t_src * total_site[3] + dt) * mom_num + mom_snk_id) * mom_num + mom_src_id) * total_site[0] + phase];
                        }
                    }
                }
            }
        }

        Long data_length_sum = total_site[3] * mom_num * mom_num * total_site[0];
        coefficent_of_t2pf.resize(data_length_sum);
        set_zero(coefficent_of_t2pf);
#pragma omp parallel for collapse(3)
        for (int dt = 0; dt < total_site[3]; ++dt)
        {
            for (int mom_src_id = 0; mom_src_id < mom_num; ++mom_src_id)
            {
                for (int mom_snk_id = 0; mom_snk_id < mom_num; ++mom_snk_id)
                {
                    for (int phase = 0; phase < total_site[0]; ++phase)
                    {
                        for (int t_src = 1; t_src < total_site[3]; ++t_src)
                        {
                            coefficent_of_t2pf_omp[(((Long)dt * mom_num + mom_snk_id) * mom_num + mom_src_id) * total_site[0] + phase] += coefficent_of_t2pf_omp[((((Long)t_src * total_site[3] + dt) * mom_num + mom_snk_id) * mom_num + mom_src_id) * total_site[0] + phase];
                        }

                        coefficent_of_t2pf[(((Long)dt * mom_num + mom_snk_id) * mom_num + mom_src_id) * total_site[0] + phase] = exp(ii * 2.0 * PI * (double)phase / (double)total_site[0]) / (double)coefficent_of_t2pf_omp[(((Long)dt * mom_num + mom_snk_id) * mom_num + mom_src_id) * total_site[0] + phase];
                    }
                }
            }
        }
    }

    inline void contract_meson_psel_two_point_block_all_psel_gys_omp_optimiztion(
        std::vector<ComplexD> &meson_two_raw,
        const std::vector<std::vector<std::vector<Meson_Block_Scalar_Psel_Psel>>> &block_meson,
        const std::string &job_tag, const int &traj,
        const PointsSelection &psel_meson,
        const Long &psel_meson_num,
        const std::vector<int> &psel_meson_num_list,
        const std::vector<int> &list_n_from_idx_meson,
        const Geometry &geo,
        const bool &is_meson_smear)
    {
        const Coordinate &total_site = geo.total_site();

        std::vector<ComplexD> gc_ts(total_site[3] * omp_get_max_threads());
        set_zero(gc_ts);
#pragma omp parallel for
        for (Long n_snk = 0; n_snk < psel_meson_num; n_snk++)
        {
            const Long xg_snk_psel_idx = n_snk;
            const Coordinate &xg_snk = psel_meson[xg_snk_psel_idx];
            const int &t_meson_snk = xg_snk[3];

            std::vector<ComplexD> c_ts;
            c_ts.resize(total_site[3]);
            set_zero(c_ts);
            for (int dt = 0; dt < total_site[3]; ++dt)
            {
                const int del_t = dt;
                const int t_meson_src = mod(t_meson_snk - del_t, total_site[3]);

                for (int idx_src = 0; idx_src < psel_meson_num_list[t_meson_src]; idx_src++)
                {
                    qassert(block_meson[xg_snk_psel_idx][dt][idx_src].n_src == idx_src);
                    c_ts[dt] += block_meson[xg_snk_psel_idx][dt][idx_src].block_pp;
                }
            }

            for (int dt = 0; dt < total_site[3]; ++dt)
            {
                gc_ts[omp_get_thread_num() * total_site[3] + dt] += c_ts[dt];
            }
        }

        std::vector<ComplexD> &m_ts = meson_two_raw;
        m_ts.resize(total_site[3]);
        set_zero(m_ts);
#pragma omp parallel for
        for (int dt = 0; dt < total_site[3]; ++dt)
        {
            for (int thread_id = 0; thread_id < omp_get_max_threads(); thread_id++)
            {
                m_ts[dt] += gc_ts[thread_id * total_site[3] + dt];
            }
        }
    }

    inline void contract_meson_psel_two_point_block_all_psel_gys_omp_optimiztion_pro(
        LatData &meson_twop_data_rs_ld,
        LatData &meson_twop_data_ld,
        const std::vector<ComplexD> &coefficent_of_t2pf,
        const std::vector<ComplexD> &coefficent_of_t2pf_remain_snk,
        const std::vector<std::vector<std::vector<Meson_Block_Scalar_Psel_Psel>>> &block_meson,
        const std::string &job_tag, const int &traj,
        const std::vector<std::vector<int>> &Momentum_Targets,
        const bool &is_meson_smear,
        const PointsSelection &psel_meson,
        const Long &psel_meson_num,
        const std::vector<int> &psel_meson_num_list,
        const std::vector<int> &list_n_from_idx_meson,
        const std::vector<std::vector<Long>> &idx_class_by_time_meson,
        const Geometry &geo)
    {
        TIMER_VERBOSE("contract_meson_psel_two_point_block_all_psel_gys_omp_optimiztion_pro");
        const Coordinate &total_site = geo.total_site();
        const int &mom_num = Momentum_Targets.size();

        const Long data_length = total_site[3] * total_site[3] * mom_num * mom_num;
        const Long data_length_sum = total_site[3] * mom_num * mom_num;

        std::vector<ComplexD> meson_two_rs;
        std::vector<ComplexD> meson_two;
        meson_two_rs.resize(data_length);
        meson_two.resize(data_length_sum);

        set_zero(meson_two_rs);
        set_zero(meson_two);

#pragma omp parallel for collapse(2)
        for (int t_meson_snk = 0; t_meson_snk < total_site[3]; ++t_meson_snk)
        {
            for (int dt = 0; dt < total_site[3]; ++dt)
            {
                for (int idx_meson_snk = 0; idx_meson_snk < psel_meson_num_list[t_meson_snk]; idx_meson_snk++)
                {
                    const long &xg_snk_psel_idx = idx_class_by_time_meson[t_meson_snk][idx_meson_snk];
                    const Coordinate &xg_snk = psel_meson[xg_snk_psel_idx];

                    const int &del_t = dt;
                    const int t_meson_src = mod(t_meson_snk - del_t, total_site[3]);

                    for (int idx_meson_src = 0; idx_meson_src < psel_meson_num_list[t_meson_src]; idx_meson_src++)
                    {
                        const long &xg_src_psel_idx = idx_class_by_time_meson[t_meson_src][idx_meson_src];
                        const Coordinate &xg_src = psel_meson[xg_src_psel_idx];

                        const Meson_Block_Scalar_Psel_Psel &block_meson_temp = block_meson[xg_snk_psel_idx][dt][idx_meson_src];
                        qassert(block_meson_temp.is_build);
                        qassert(block_meson_temp.xg_snk_psel_idx == xg_snk_psel_idx);
                        qassert(block_meson_temp.xg_src_psel_idx == xg_src_psel_idx);
                        qassert(block_meson_temp.n_snk == idx_meson_snk);
                        qassert(block_meson_temp.n_src == idx_meson_src);

                        for (int mom_src_id = 0; mom_src_id < mom_num; ++mom_src_id)
                        {
                            for (int mom_snk_id = 0; mom_snk_id < mom_num; ++mom_snk_id)
                            {
                                Long phase_temp = space_dot(Momentum_Targets[mom_src_id], xg_src) - space_dot(Momentum_Targets[mom_snk_id], xg_snk);
                                int phase = mod(phase_temp, total_site[0]);

                                const ComplexD &normal_coeff = coefficent_of_t2pf_remain_snk[(((t_meson_src * total_site[3] + dt) * mom_num + mom_snk_id) * mom_num + mom_src_id) * total_site[0] + phase];

                                meson_two_rs[(((t_meson_src * total_site[3] + dt) * mom_num + mom_snk_id) * mom_num + mom_src_id)] += normal_coeff * block_meson_temp.block_pp;
                            }
                        }
                    }
                }
            }
        }

#pragma omp parallel for
        for (int dt = 0; dt < total_site[3]; ++dt)
        {
            for (int t_meson_snk = 0; t_meson_snk < total_site[3]; ++t_meson_snk)
            {
                for (int idx_meson_snk = 0; idx_meson_snk < psel_meson_num_list[t_meson_snk]; idx_meson_snk++)
                {
                    const long &xg_snk_psel_idx = idx_class_by_time_meson[t_meson_snk][idx_meson_snk];
                    const Coordinate &xg_snk = psel_meson[xg_snk_psel_idx];

                    const int &del_t = dt;
                    const int t_meson_src = mod(t_meson_snk - del_t, total_site[3]);

                    for (int idx_meson_src = 0; idx_meson_src < psel_meson_num_list[t_meson_src]; idx_meson_src++)
                    {
                        const long &xg_src_psel_idx = idx_class_by_time_meson[t_meson_src][idx_meson_src];
                        const Coordinate &xg_src = psel_meson[xg_src_psel_idx];

                        const Meson_Block_Scalar_Psel_Psel &block_meson_temp = block_meson[xg_snk_psel_idx][dt][idx_meson_src];
                        qassert(block_meson_temp.is_build);
                        qassert(block_meson_temp.xg_snk_psel_idx == xg_snk_psel_idx);
                        qassert(block_meson_temp.xg_src_psel_idx == xg_src_psel_idx);
                        qassert(block_meson_temp.n_snk == idx_meson_snk);
                        qassert(block_meson_temp.n_src == idx_meson_src);

                        for (int mom_src_id = 0; mom_src_id < mom_num; ++mom_src_id)
                        {
                            for (int mom_snk_id = 0; mom_snk_id < mom_num; ++mom_snk_id)
                            {
                                Long phase_temp = space_dot(Momentum_Targets[mom_src_id], xg_src) - space_dot(Momentum_Targets[mom_snk_id], xg_snk);
                                int phase = mod(phase_temp, total_site[0]);

                                const ComplexD &normal_coeff = coefficent_of_t2pf[((dt * mom_num + mom_snk_id) * mom_num + mom_src_id) * total_site[0] + phase];

                                meson_two[((dt * mom_num + mom_snk_id) * mom_num + mom_src_id)] += normal_coeff * block_meson_temp.block_pp;
                            }
                        }
                    }
                }
            }
        }
#pragma omp parallel for collapse(4)
        for (int t_src = 0; t_src < total_site[3]; ++t_src)
        {
            for (int dt = 0; dt < total_site[3]; ++dt)
            {
                for (int mom_src_id = 0; mom_src_id < mom_num; ++mom_src_id)
                {
                    for (int mom_snk_id = 0; mom_snk_id < mom_num; ++mom_snk_id)
                    {
                        Vector<ComplexD> meson_np_twop_sum_rs = lat_data_complex_get(meson_twop_data_rs_ld, make_array(t_src, dt, mom_src_id, mom_snk_id));
                        meson_np_twop_sum_rs[0] = meson_two_rs[(((t_src * total_site[3] + dt) * mom_num + mom_snk_id) * mom_num + mom_src_id)];
                    }
                }
            }
        }

#pragma omp parallel for collapse(3)
        for (int dt = 0; dt < total_site[3]; ++dt)
        {
            for (int mom_src_id = 0; mom_src_id < mom_num; ++mom_src_id)
            {
                for (int mom_snk_id = 0; mom_snk_id < mom_num; ++mom_snk_id)
                {
                    Vector<ComplexD> meson_np_twop_sum = lat_data_complex_get(meson_twop_data_ld, make_array(dt, mom_src_id, mom_snk_id));

                    meson_np_twop_sum[0] = meson_two[((dt * mom_num + mom_snk_id) * mom_num + mom_src_id)];
                }
            }
        }
    }

    inline void contract_baryon_psel_two_point_block_all_psel_gys_omp_optimiztion(
        std::vector<ComplexD> &baryon_two_raw,
        const std::vector<std::vector<std::vector<Baryon_Block_Scalar_Psel_Psel>>> &block_baryon,
        const std::string &job_tag, const int &traj,
        const PointsSelection &psel_baryon,
        const Long &psel_baryon_num,
        const std::vector<int> &psel_baryon_num_list,
        const std::vector<int> &list_n_from_idx_baryon,
        const Geometry &geo,
        const bool &is_baryon_smear)
    {
        const Coordinate &total_site = geo.total_site();

        std::vector<ComplexD> gc_ts(total_site[3] * omp_get_max_threads());
        set_zero(gc_ts);
#pragma omp parallel for
        for (Long n_snk = 0; n_snk < psel_baryon_num; n_snk++)
        {
            const Long xg_snk_psel_idx = n_snk;
            const Coordinate &xg_snk = psel_baryon[xg_snk_psel_idx];
            const int &t_baryon_snk = xg_snk[3];

            std::vector<ComplexD> c_ts;
            c_ts.resize(total_site[3]);
            set_zero(c_ts);

            double anti_period = 0;
            for (int dt = 0; dt < total_site[3]; ++dt)
            {
                const int del_t = dt;
                const int t_baryon_src = mod(t_baryon_snk - del_t, total_site[3]);

                if (t_baryon_src <= t_baryon_snk)
                {
                    anti_period = 1.0;
                }
                else
                {
                    anti_period = -1.0;
                }
                for (int idx_src = 0; idx_src < psel_baryon_num_list[t_baryon_src]; idx_src++)
                {
                    qassert(block_baryon[xg_snk_psel_idx][dt][idx_src].n_src == idx_src);
                    c_ts[dt] += anti_period * block_baryon[xg_snk_psel_idx][dt][idx_src].block_pp[0];
                    c_ts[dt] -= anti_period * block_baryon[xg_snk_psel_idx][dt][idx_src].block_pp[1];
                }
            }

            for (int dt = 0; dt < total_site[3]; ++dt)
            {
                gc_ts[omp_get_thread_num() * total_site[3] + dt] += c_ts[dt];
            }
        }

        std::vector<ComplexD> &m_ts = baryon_two_raw;
        m_ts.resize(total_site[3]);
        set_zero(m_ts);
#pragma omp parallel for
        for (int dt = 0; dt < total_site[3]; ++dt)
        {
            for (int thread_id = 0; thread_id < omp_get_max_threads(); thread_id++)
            {
                m_ts[dt] += gc_ts[thread_id * total_site[3] + dt];
            }
        }
    }

    inline void contract_baryon_psel_two_point_block_all_psel_gys_omp_optimiztion_pro(
        LatData &baryon_twop_data_rs_ld,
        LatData &baryon_twop_data_ld,
        const std::vector<ComplexD> &coefficent_of_t2pf,
        const std::vector<ComplexD> &coefficent_of_t2pf_remain_snk,
        const std::vector<std::vector<std::vector<Baryon_Block_Scalar_Psel_Psel>>> &block_baryon,
        const std::string &job_tag, const int &traj,
        const std::vector<std::vector<int>> &Momentum_Targets,
        const bool &is_baryon_smear,
        const PointsSelection &psel_baryon,
        const Long &psel_baryon_num,
        const std::vector<int> &psel_baryon_num_list,
        const std::vector<int> &list_n_from_idx_baryon,
        const std::vector<std::vector<Long>> &idx_class_by_time_baryon,
        const Geometry &geo)
    {
        TIMER_VERBOSE("contract_baryon_psel_two_point_block_all_psel_gys_omp_optimiztion_pro");
        const Coordinate &total_site = geo.total_site();
        const int &mom_num = Momentum_Targets.size();

        const Long data_length = total_site[3] * total_site[3] * mom_num * mom_num;
        const Long data_length_sum = total_site[3] * mom_num * mom_num;

        std::vector<ComplexD> baryon_two_rs;
        std::vector<ComplexD> baryon_two;
        baryon_two_rs.resize(data_length);
        baryon_two.resize(data_length_sum);

        set_zero(baryon_two_rs);
        set_zero(baryon_two);

#pragma omp parallel for collapse(2)
        for (int t_baryon_snk = 0; t_baryon_snk < total_site[3]; ++t_baryon_snk)
        {
            for (int dt = 0; dt < total_site[3]; ++dt)
            {
                for (int idx_baryon_snk = 0; idx_baryon_snk < psel_baryon_num_list[t_baryon_snk]; idx_baryon_snk++)
                {
                    const long &xg_snk_psel_idx = idx_class_by_time_baryon[t_baryon_snk][idx_baryon_snk];
                    const Coordinate &xg_snk = psel_baryon[xg_snk_psel_idx];

                    const int &del_t = dt;
                    const int t_baryon_src = mod(t_baryon_snk - del_t, total_site[3]);

                    ComplexD anti_period = 0.0;

                    if (t_baryon_src <= t_baryon_snk)
                    {
                        anti_period = 1.0;
                    }
                    else
                    {
                        anti_period = -1.0;
                    }

                    qassert(anti_period != 0.0);

                    for (int idx_baryon_src = 0; idx_baryon_src < psel_baryon_num_list[t_baryon_src]; idx_baryon_src++)
                    {
                        const long &xg_src_psel_idx = idx_class_by_time_baryon[t_baryon_src][idx_baryon_src];
                        const Coordinate &xg_src = psel_baryon[xg_src_psel_idx];

                        const Baryon_Block_Scalar_Psel_Psel &block_baryon_temp = block_baryon[xg_snk_psel_idx][dt][idx_baryon_src];
                        qassert(block_baryon_temp.is_build);
                        qassert(block_baryon_temp.xg_snk_psel_idx == xg_snk_psel_idx);
                        qassert(block_baryon_temp.xg_src_psel_idx == xg_src_psel_idx);
                        qassert(block_baryon_temp.n_snk == idx_baryon_snk);
                        qassert(block_baryon_temp.n_src == idx_baryon_src);

                        for (int mom_src_id = 0; mom_src_id < mom_num; ++mom_src_id)
                        {
                            for (int mom_snk_id = 0; mom_snk_id < mom_num; ++mom_snk_id)
                            {
                                Long phase_temp = space_dot(Momentum_Targets[mom_src_id], xg_src) - space_dot(Momentum_Targets[mom_snk_id], xg_snk);
                                int phase = mod(phase_temp, total_site[0]);

                                const ComplexD &normal_coeff = coefficent_of_t2pf_remain_snk[(((t_baryon_src * total_site[3] + dt) * mom_num + mom_snk_id) * mom_num + mom_src_id) * total_site[0] + phase];

                                baryon_two_rs[(((t_baryon_src * total_site[3] + dt) * mom_num + mom_snk_id) * mom_num + mom_src_id)] += anti_period * normal_coeff * block_baryon_temp.block_pp[0];
                                baryon_two_rs[(((t_baryon_src * total_site[3] + dt) * mom_num + mom_snk_id) * mom_num + mom_src_id)] -= anti_period * normal_coeff * block_baryon_temp.block_pp[1];
                            }
                        }
                    }
                }
            }
        }

#pragma omp parallel for
        for (int dt = 0; dt < total_site[3]; ++dt)
        {
            for (int t_baryon_snk = 0; t_baryon_snk < total_site[3]; ++t_baryon_snk)
            {
                for (int idx_baryon_snk = 0; idx_baryon_snk < psel_baryon_num_list[t_baryon_snk]; idx_baryon_snk++)
                {
                    const long &xg_snk_psel_idx = idx_class_by_time_baryon[t_baryon_snk][idx_baryon_snk];
                    const Coordinate &xg_snk = psel_baryon[xg_snk_psel_idx];

                    const int &del_t = dt;
                    const int t_baryon_src = mod(t_baryon_snk - del_t, total_site[3]);

                    ComplexD anti_period = 0.0;

                    if (t_baryon_src <= t_baryon_snk)
                    {
                        anti_period = 1.0;
                    }
                    else
                    {
                        anti_period = -1.0;
                    }

                    qassert(anti_period != 0.0);

                    for (int idx_baryon_src = 0; idx_baryon_src < psel_baryon_num_list[t_baryon_src]; idx_baryon_src++)
                    {
                        const long &xg_src_psel_idx = idx_class_by_time_baryon[t_baryon_src][idx_baryon_src];
                        const Coordinate &xg_src = psel_baryon[xg_src_psel_idx];

                        const Baryon_Block_Scalar_Psel_Psel &block_baryon_temp = block_baryon[xg_snk_psel_idx][dt][idx_baryon_src];
                        qassert(block_baryon_temp.is_build);
                        qassert(block_baryon_temp.xg_snk_psel_idx == xg_snk_psel_idx);
                        qassert(block_baryon_temp.xg_src_psel_idx == xg_src_psel_idx);
                        qassert(block_baryon_temp.n_snk == idx_baryon_snk);
                        qassert(block_baryon_temp.n_src == idx_baryon_src);

                        for (int mom_src_id = 0; mom_src_id < mom_num; ++mom_src_id)
                        {
                            for (int mom_snk_id = 0; mom_snk_id < mom_num; ++mom_snk_id)
                            {
                                Long phase_temp = space_dot(Momentum_Targets[mom_src_id], xg_src) - space_dot(Momentum_Targets[mom_snk_id], xg_snk);
                                int phase = mod(phase_temp, total_site[0]);

                                const ComplexD &normal_coeff = coefficent_of_t2pf[((dt * mom_num + mom_snk_id) * mom_num + mom_src_id) * total_site[0] + phase];

                                baryon_two[((dt * mom_num + mom_snk_id) * mom_num + mom_src_id)] += anti_period * normal_coeff * block_baryon_temp.block_pp[0];
                                baryon_two[((dt * mom_num + mom_snk_id) * mom_num + mom_src_id)] -= anti_period * normal_coeff * block_baryon_temp.block_pp[1];
                            }
                        }
                    }
                }
            }
        }
#pragma omp parallel for collapse(4)
        for (int t_src = 0; t_src < total_site[3]; ++t_src)
        {
            for (int dt = 0; dt < total_site[3]; ++dt)
            {
                for (int mom_src_id = 0; mom_src_id < mom_num; ++mom_src_id)
                {
                    for (int mom_snk_id = 0; mom_snk_id < mom_num; ++mom_snk_id)
                    {
                        Vector<ComplexD> baryon_np_twop_sum_rs = lat_data_complex_get(baryon_twop_data_rs_ld, make_array(t_src, dt, mom_src_id, mom_snk_id));
                        baryon_np_twop_sum_rs[0] = baryon_two_rs[(((t_src * total_site[3] + dt) * mom_num + mom_snk_id) * mom_num + mom_src_id)];
                    }
                }
            }
        }

#pragma omp parallel for collapse(3)
        for (int dt = 0; dt < total_site[3]; ++dt)
        {
            for (int mom_src_id = 0; mom_src_id < mom_num; ++mom_src_id)
            {
                for (int mom_snk_id = 0; mom_snk_id < mom_num; ++mom_snk_id)
                {
                    Vector<ComplexD> baryon_np_twop_sum = lat_data_complex_get(baryon_twop_data_ld, make_array(dt, mom_src_id, mom_snk_id));

                    baryon_np_twop_sum[0] = baryon_two[((dt * mom_num + mom_snk_id) * mom_num + mom_src_id)];
                }
            }
        }
    }

    inline void contract_baryon_psel_two_point_block_no_projection_all_psel_gys_omp_optimiztion(
        std::vector<SpinMatrix> &baryon_two_no_projection_raw,
        const std::vector<std::vector<std::vector<Baryon_Block_2PF_No_Projection_Psel_Psel>>> &block_baryon_no_projection,
        const std::string &job_tag, const int &traj,
        const PointsSelection &psel_baryon,
        const Long &psel_baryon_num,
        const std::vector<int> &psel_baryon_num_list,
        const std::vector<int> &list_n_from_idx_baryon,
        const Geometry &geo,
        const bool &is_baryon_smear)
    {
        const Coordinate &total_site = geo.total_site();

        std::vector<SpinMatrix> gc_ts(total_site[3] * omp_get_max_threads());
        for (size_t i = 0; i < gc_ts.size(); i++)
        {
            set_zero(gc_ts[i]);
        }
#pragma omp parallel for
        for (Long n_snk = 0; n_snk < psel_baryon_num; n_snk++)
        {
            const Long xg_snk_psel_idx = n_snk;
            const Coordinate &xg_snk = psel_baryon[xg_snk_psel_idx];
            const int &t_baryon_snk = xg_snk[3];

            std::vector<SpinMatrix> c_ts;
            c_ts.resize(total_site[3]);
            for (size_t i = 0; i < c_ts.size(); i++)
            {
                set_zero(c_ts[i]);
            }

            double anti_period = 0;
            for (int dt = 0; dt < total_site[3]; ++dt)
            {
                const int del_t = dt;
                const int t_baryon_src = mod(t_baryon_snk - del_t, total_site[3]);

                if (t_baryon_src <= t_baryon_snk)
                {
                    anti_period = 1.0;
                }
                else
                {
                    anti_period = -1.0;
                }

                for (int idx_src = 0; idx_src < psel_baryon_num_list[t_baryon_src]; idx_src++)
                {
                    qassert(block_baryon_no_projection[xg_snk_psel_idx][dt][idx_src].n_src == idx_src);
                    c_ts[dt] += (ComplexD)anti_period * block_baryon_no_projection[xg_snk_psel_idx][dt][idx_src].block_pp[0];
                    c_ts[dt] -= (ComplexD)anti_period * block_baryon_no_projection[xg_snk_psel_idx][dt][idx_src].block_pp[1];
                }
            }

            for (int dt = 0; dt < total_site[3]; ++dt)
            {
                gc_ts[omp_get_thread_num() * total_site[3] + dt] += c_ts[dt];
            }
        }

        std::vector<SpinMatrix> &m_ts = baryon_two_no_projection_raw;
        m_ts.resize(total_site[3]);
        for (size_t i = 0; i < m_ts.size(); i++)
        {
            set_zero(m_ts[i]);
        }

#pragma omp parallel for
        for (int dt = 0; dt < total_site[3]; ++dt)
        {
            for (int thread_id = 0; thread_id < omp_get_max_threads(); thread_id++)
            {
                m_ts[dt] += gc_ts[thread_id * total_site[3] + dt];
            }
        }
    }

    inline void contract_baryon_psel_two_point_block_no_projection_all_psel_gys_omp_optimiztion_pro(
        LatData &baryon_twop_data_no_projection_rs_ld,
        LatData &baryon_twop_data_no_projection_ld,
        const std::vector<ComplexD> &coefficent_of_t2pf,
        const std::vector<ComplexD> &coefficent_of_t2pf_remain_snk,
        const std::vector<std::vector<std::vector<Baryon_Block_2PF_No_Projection_Psel_Psel>>> &block_baryon_no_projection,
        const std::string &job_tag, const int &traj,
        const std::vector<std::vector<int>> &Momentum_Targets,
        const bool &is_baryon_smear,
        const PointsSelection &psel_baryon,
        const Long &psel_baryon_num,
        const std::vector<int> &psel_baryon_num_list,
        const std::vector<int> &list_n_from_idx_baryon,
        const std::vector<std::vector<Long>> &idx_class_by_time_baryon,
        const Geometry &geo)
    {
        TIMER_VERBOSE("contract_baryon_psel_two_point_block_no_projection_all_psel_gys_omp_optimiztion_pro");
        const Coordinate &total_site = geo.total_site();
        const int &mom_num = Momentum_Targets.size();

        const Long data_length = total_site[3] * total_site[3] * mom_num * mom_num;
        const Long data_length_sum = total_site[3] * mom_num * mom_num;

        std::vector<SpinMatrix> baryon_two_no_projection_rs;
        std::vector<SpinMatrix> baryon_two_no_projection;
        baryon_two_no_projection_rs.resize(data_length);
        baryon_two_no_projection.resize(data_length_sum);

        set_zero(baryon_two_no_projection_rs);
        set_zero(baryon_two_no_projection);

#pragma omp parallel for collapse(2)
        for (int t_baryon_snk = 0; t_baryon_snk < total_site[3]; ++t_baryon_snk)
        {
            for (int dt = 0; dt < total_site[3]; ++dt)
            {
                for (int idx_baryon_snk = 0; idx_baryon_snk < psel_baryon_num_list[t_baryon_snk]; idx_baryon_snk++)
                {
                    const long &xg_snk_psel_idx = idx_class_by_time_baryon[t_baryon_snk][idx_baryon_snk];
                    const Coordinate &xg_snk = psel_baryon[xg_snk_psel_idx];

                    const int &del_t = dt;
                    const int t_baryon_src = mod(t_baryon_snk - del_t, total_site[3]);

                    ComplexD anti_period = 0.0;

                    if (t_baryon_src <= t_baryon_snk)
                    {
                        anti_period = 1.0;
                    }
                    else
                    {
                        anti_period = -1.0;
                    }

                    for (int idx_baryon_src = 0; idx_baryon_src < psel_baryon_num_list[t_baryon_src]; idx_baryon_src++)
                    {
                        const long &xg_src_psel_idx = idx_class_by_time_baryon[t_baryon_src][idx_baryon_src];
                        const Coordinate &xg_src = psel_baryon[xg_src_psel_idx];

                        const Baryon_Block_2PF_No_Projection_Psel_Psel &block_baryon_no_projection_temp = block_baryon_no_projection[xg_snk_psel_idx][dt][idx_baryon_src];
                        qassert(block_baryon_no_projection_temp.is_build);
                        qassert(block_baryon_no_projection_temp.xg_snk_psel_idx == xg_snk_psel_idx);
                        qassert(block_baryon_no_projection_temp.xg_src_psel_idx == xg_src_psel_idx);
                        qassert(block_baryon_no_projection_temp.n_snk == idx_baryon_snk);
                        qassert(block_baryon_no_projection_temp.n_src == idx_baryon_src);

                        for (int mom_src_id = 0; mom_src_id < mom_num; ++mom_src_id)
                        {
                            for (int mom_snk_id = 0; mom_snk_id < mom_num; ++mom_snk_id)
                            {
                                Long phase_temp = space_dot(Momentum_Targets[mom_src_id], xg_src) - space_dot(Momentum_Targets[mom_snk_id], xg_snk);
                                int phase = mod(phase_temp, total_site[0]);

                                const ComplexD &normal_coeff = coefficent_of_t2pf_remain_snk[(((t_baryon_src * total_site[3] + dt) * mom_num + mom_snk_id) * mom_num + mom_src_id) * total_site[0] + phase];

                                baryon_two_no_projection_rs[(((t_baryon_src * total_site[3] + dt) * mom_num + mom_snk_id) * mom_num + mom_src_id)] += anti_period * normal_coeff * block_baryon_no_projection_temp.block_pp[0];
                                baryon_two_no_projection_rs[(((t_baryon_src * total_site[3] + dt) * mom_num + mom_snk_id) * mom_num + mom_src_id)] -= anti_period * normal_coeff * block_baryon_no_projection_temp.block_pp[1];
                            }
                        }
                    }
                }
            }
        }

#pragma omp parallel for
        for (int dt = 0; dt < total_site[3]; ++dt)
        {
            for (int t_baryon_snk = 0; t_baryon_snk < total_site[3]; ++t_baryon_snk)
            {
                for (int idx_baryon_snk = 0; idx_baryon_snk < psel_baryon_num_list[t_baryon_snk]; idx_baryon_snk++)
                {
                    const long &xg_snk_psel_idx = idx_class_by_time_baryon[t_baryon_snk][idx_baryon_snk];
                    const Coordinate &xg_snk = psel_baryon[xg_snk_psel_idx];

                    const int &del_t = dt;
                    const int t_baryon_src = mod(t_baryon_snk - del_t, total_site[3]);

                    ComplexD anti_period = 0.0;

                    if (t_baryon_src <= t_baryon_snk)
                    {
                        anti_period = 1.0;
                    }
                    else
                    {
                        anti_period = -1.0;
                    }

                    for (int idx_baryon_src = 0; idx_baryon_src < psel_baryon_num_list[t_baryon_src]; idx_baryon_src++)
                    {
                        const long &xg_src_psel_idx = idx_class_by_time_baryon[t_baryon_src][idx_baryon_src];
                        const Coordinate &xg_src = psel_baryon[xg_src_psel_idx];

                        const Baryon_Block_2PF_No_Projection_Psel_Psel &block_baryon_no_projection_temp = block_baryon_no_projection[xg_snk_psel_idx][dt][idx_baryon_src];
                        qassert(block_baryon_no_projection_temp.is_build);
                        qassert(block_baryon_no_projection_temp.xg_snk_psel_idx == xg_snk_psel_idx);
                        qassert(block_baryon_no_projection_temp.xg_src_psel_idx == xg_src_psel_idx);
                        qassert(block_baryon_no_projection_temp.n_snk == idx_baryon_snk);
                        qassert(block_baryon_no_projection_temp.n_src == idx_baryon_src);

                        for (int mom_src_id = 0; mom_src_id < mom_num; ++mom_src_id)
                        {
                            for (int mom_snk_id = 0; mom_snk_id < mom_num; ++mom_snk_id)
                            {
                                Long phase_temp = space_dot(Momentum_Targets[mom_src_id], xg_src) - space_dot(Momentum_Targets[mom_snk_id], xg_snk);
                                int phase = mod(phase_temp, total_site[0]);

                                const ComplexD &normal_coeff = coefficent_of_t2pf[((dt * mom_num + mom_snk_id) * mom_num + mom_src_id) * total_site[0] + phase];

                                baryon_two_no_projection[((dt * mom_num + mom_snk_id) * mom_num + mom_src_id)] += anti_period * normal_coeff * block_baryon_no_projection_temp.block_pp[0];
                                baryon_two_no_projection[((dt * mom_num + mom_snk_id) * mom_num + mom_src_id)] -= anti_period * normal_coeff * block_baryon_no_projection_temp.block_pp[1];
                            }
                        }
                    }
                }
            }
        }
#pragma omp parallel for collapse(4)
        for (int t_src = 0; t_src < total_site[3]; ++t_src)
        {
            for (int dt = 0; dt < total_site[3]; ++dt)
            {
                for (int mom_src_id = 0; mom_src_id < mom_num; ++mom_src_id)
                {
                    for (int mom_snk_id = 0; mom_snk_id < mom_num; ++mom_snk_id)
                    {
                        Vector<ComplexD> baryon_np_twop_sum_rs = lat_data_complex_get(baryon_twop_data_no_projection_rs_ld, make_array(t_src, dt, mom_src_id, mom_snk_id));
                        for (int mu1 = 0; mu1 < 4; mu1++)
                        {
                            for (int mu2 = 0; mu2 < 4; mu2++)
                            {
                                baryon_np_twop_sum_rs[mu1 * 4 + mu2] = baryon_two_no_projection_rs[(((t_src * total_site[3] + dt) * mom_num + mom_snk_id) * mom_num + mom_src_id)](mu1, mu2);
                            }
                        }
                    }
                }
            }
        }

#pragma omp parallel for collapse(3)
        for (int dt = 0; dt < total_site[3]; ++dt)
        {
            for (int mom_src_id = 0; mom_src_id < mom_num; ++mom_src_id)
            {
                for (int mom_snk_id = 0; mom_snk_id < mom_num; ++mom_snk_id)
                {
                    Vector<ComplexD> baryon_np_twop_sum = lat_data_complex_get(baryon_twop_data_no_projection_ld, make_array(dt, mom_src_id, mom_snk_id));

                    for (int mu1 = 0; mu1 < 4; mu1++)
                    {
                        for (int mu2 = 0; mu2 < 4; mu2++)
                        {
                            baryon_np_twop_sum[mu1 * 4 + mu2] = baryon_two_no_projection[((dt * mom_num + mom_snk_id) * mom_num + mom_src_id)](mu1, mu2);
                        }
                    }
                }
            }
        }
    }
    // inline LatData mk_two_point_table_omp(const Coordinate &total_site)
    // {
    //     LatData ld;
    //     ld.info.push_back(lat_dim_number("omp", 0, omp_get_max_threads() - 1));
    //     ld.info.push_back(lat_dim_number("tsrc", 0, total_site[3] - 1));
    //     ld.info.push_back(lat_dim_number("tsep", 0, total_site[3] - 1));
    //     ld.info.push_back(lat_dim_re_im());
    //     lat_data_alloc(ld);
    //     set_zero(ld);
    //     return ld;
    // }

    inline LatData mk_two_point_remain_snk_table(const Coordinate &total_site)
    {
        LatData ld;
        ld.info.push_back(lat_dim_number("tsrc", 0, total_site[3] - 1));
        ld.info.push_back(lat_dim_number("tsep", 0, total_site[3] - 1));
        ld.info.push_back(lat_dim_re_im());
        lat_data_alloc(ld);
        set_zero(ld);
        return ld;
    }

    inline LatData mk_two_point_sum_table(const Coordinate &total_site)
    {
        LatData ld;
        ld.info.push_back(lat_dim_number("tsep", 0, total_site[3] - 1));
        ld.info.push_back(lat_dim_re_im());
        lat_data_alloc(ld);
        set_zero(ld);
        return ld;
    }

    inline LatData mk_two_point_momentum_remain_snk_table(const Coordinate &total_site, const std::vector<std::vector<int>> &Momentum_Targets)
    {
        LatData ld;
        ld.info.push_back(lat_dim_number("tsrc", 0, total_site[3] - 1));
        ld.info.push_back(lat_dim_number("tsep", 0, total_site[3] - 1));
        ld.info.push_back(lat_dim_number("momentum_mode_snk", 0, Momentum_Targets.size() - 1));
        ld.info.push_back(lat_dim_number("momentum_mode_src", 0, Momentum_Targets.size() - 1));
        ld.info.push_back(lat_dim_re_im());
        lat_data_alloc(ld);
        set_zero(ld);
        return ld;
    }

    inline LatData mk_two_point_momentum_sum_table(const Coordinate &total_site, const std::vector<std::vector<int>> &Momentum_Targets)
    {
        LatData ld;
        ld.info.push_back(lat_dim_number("tsep", 0, total_site[3] - 1));
        ld.info.push_back(lat_dim_number("momentum_mode_snk", 0, Momentum_Targets.size() - 1));
        ld.info.push_back(lat_dim_number("momentum_mode_src", 0, Momentum_Targets.size() - 1));
        ld.info.push_back(lat_dim_re_im());
        lat_data_alloc(ld);
        set_zero(ld);
        return ld;
    }

    inline LatData mk_two_point_no_projection_remain_snk_table(const Coordinate &total_site)
    {
        LatData ld;
        ld.info.push_back(lat_dim_number("tsrc", 0, total_site[3] - 1));
        ld.info.push_back(lat_dim_number("tsep", 0, total_site[3] - 1));
        ld.info.push_back(lat_dim_number("mu", 0, 3));
        ld.info.push_back(lat_dim_number("nu", 0, 3));
        ld.info.push_back(lat_dim_re_im());
        lat_data_alloc(ld);
        set_zero(ld);
        return ld;
    }

    inline LatData mk_two_point_no_projection_sum_table(const Coordinate &total_site)
    {
        LatData ld;
        ld.info.push_back(lat_dim_number("tsep", 0, total_site[3] - 1));
        ld.info.push_back(lat_dim_number("mu", 0, 3));
        ld.info.push_back(lat_dim_number("nu", 0, 3));
        ld.info.push_back(lat_dim_re_im());
        lat_data_alloc(ld);
        set_zero(ld);
        return ld;
    }

    inline LatData mk_two_point_momentum_no_projection_remain_snk_table(const Coordinate &total_site, const std::vector<std::vector<int>> &Momentum_Targets)
    {
        LatData ld;
        ld.info.push_back(lat_dim_number("tsrc", 0, total_site[3] - 1));
        ld.info.push_back(lat_dim_number("tsep", 0, total_site[3] - 1));
        ld.info.push_back(lat_dim_number("momentum_mode_snk", 0, Momentum_Targets.size() - 1));
        ld.info.push_back(lat_dim_number("momentum_mode_src", 0, Momentum_Targets.size() - 1));
        ld.info.push_back(lat_dim_number("mu", 0, 3));
        ld.info.push_back(lat_dim_number("nu", 0, 3));
        ld.info.push_back(lat_dim_re_im());
        lat_data_alloc(ld);
        set_zero(ld);
        return ld;
    }

    inline LatData mk_two_point_momentum_no_projection_sum_table(const Coordinate &total_site, const std::vector<std::vector<int>> &Momentum_Targets)
    {
        LatData ld;
        ld.info.push_back(lat_dim_number("tsep", 0, total_site[3] - 1));
        ld.info.push_back(lat_dim_number("momentum_mode_snk", 0, Momentum_Targets.size() - 1));
        ld.info.push_back(lat_dim_number("momentum_mode_src", 0, Momentum_Targets.size() - 1));
        ld.info.push_back(lat_dim_number("mu", 0, 3));
        ld.info.push_back(lat_dim_number("nu", 0, 3));
        ld.info.push_back(lat_dim_re_im());
        lat_data_alloc(ld);
        set_zero(ld);
        return ld;
    }

    inline LatData contract_two_point_function_psel_correction(const LatData &ld0, const std::vector<int> &psel_num_list,
                                                               const Geometry &geo)
    {
        const Coordinate &total_site = geo.total_site();
        LatData ld = mk_two_point_sum_table(total_site);
        set_zero(ld);
        Vector<ComplexD> c_sum = lat_data_cget(ld);
        for (int tsep = 0; tsep < total_site[3]; ++tsep)
        {
            int pselsum = 0;
            for (int t_snk = 0; t_snk < total_site[3]; ++t_snk)
            {
                const int t_src = mod(t_snk - tsep, total_site[3]);
                pselsum += psel_num_list[t_snk] * psel_num_list[t_src];
                const ComplexD &c_t_snk_tsep = lat_data_complex_get_const(ld0, make_array(t_snk, tsep))[0];
                c_sum[tsep] += c_t_snk_tsep;
            }
            const double coef = -1.0 / (double)pselsum;
            //  * (double)total_site[0] * (double)total_site[1] * (double)total_site[2];

            c_sum[tsep] *= coef;
        }
        return ld;
    }

    inline void contract_two_point_function_psel_correction(
        LatData &meson_twop_data,
        LatData &baryon_twop_data,
        LatData &baryon_twop_data_no_projection,
        const std::vector<ComplexD> &meson_two_raw,
        const std::vector<ComplexD> &baryon_two_raw,
        const std::vector<SpinMatrix> &baryon_two_no_projection_raw,
        const std::vector<int> &psel_meson_num_list,
        const std::vector<int> &psel_baryon_num_list,
        const Coordinate &total_site)
    {
        set_zero(meson_twop_data);
        set_zero(baryon_twop_data);
        set_zero(baryon_twop_data_no_projection);

        Vector<ComplexD> meson_twop_sum = lat_data_cget(meson_twop_data);
        Vector<ComplexD> baryon_twop_sum = lat_data_cget(baryon_twop_data);
        for (int tsep = 0; tsep < total_site[3]; ++tsep)
        {

            Vector<ComplexD> baryon_twop_sum_no_projection = lat_data_complex_get(baryon_twop_data_no_projection, make_array(tsep));

            int psel_baryon_sum = 0;
            int psel_meson_sum = 0;

            for (int t_snk = 0; t_snk < total_site[3]; ++t_snk)
            {
                const int t_src = mod(t_snk - tsep, total_site[3]);
                psel_meson_sum += psel_meson_num_list[t_snk] * psel_meson_num_list[t_src];
                psel_baryon_sum += psel_baryon_num_list[t_snk] * psel_baryon_num_list[t_src];
            }
            const double coef_meson = -1.0 / (double)psel_meson_sum * (double)total_site[0] * (double)total_site[1] * (double)total_site[2];
            const double coef_baryon = -1.0 / (double)psel_baryon_sum * (double)total_site[0] * (double)total_site[1] * (double)total_site[2];

            meson_twop_sum[tsep] = coef_meson * meson_two_raw[tsep];
            baryon_twop_sum[tsep] = coef_baryon * baryon_two_raw[tsep];
            for (int mu = 0; mu < 4; mu++)
            {
                for (int nu = 0; nu < 4; nu++)
                {
                    baryon_twop_sum_no_projection[mu * 4 + nu] = coef_baryon * baryon_two_no_projection_raw[tsep](mu, nu);
                }
            }
        }
    }
}