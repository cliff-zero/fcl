#pragma once

#include "data-load.h"
#include "baryon-pion-block-acc.h"
#include "baryon-pion-coefficient.h"
#include "compute-utils.h"

namespace qlat
{
    inline void contract_baryon_pi_unshifted_psel_acc_x_typeI(Vector<Complex> &baryon_pi_data_temp, const Complex coef, const WilsonMatrix &wm_x_snk,
                                                              const WilsonMatrix &wm_x_y, const Baryon_Block_Psel_Psel &block)
    {
        const SpinMatrix &gamma5 = SpinMatrixConstants::get_gamma5();
        const WilsonMatrix wm_snk_x =
            gamma5 * (WilsonMatrix)matrix_adjoint(wm_x_snk) * gamma5;

        for (int gram = 0; gram < 5; gram++)
        {
            const WilsonMatrix wm_pp_y_tsrc_tsnk_x = block.block_pp[gram] * wm_snk_x;
            const WilsonMatrix wm_pp = wm_x_y * gamma5 * wm_pp_y_tsrc_tsnk_x;
            baryon_pi_data_temp[gram] += coef * matrix_trace(wm_pp, gamma5);
        }
        return;
    }

    inline void contract_baryon_pi_unshifted_psel_acc_x_typeII(
        Vector<Complex> &baryon_pi_data_temp, const Complex coef, const WilsonMatrix &wm_x_snk,
        const WilsonMatrix &wm_x_src, const Baryon_Block_Psel_Psel &block)
    {
        const SpinMatrix &gamma5 = SpinMatrixConstants::get_gamma5();
        const WilsonMatrix wm_snk_x =
            gamma5 * (WilsonMatrix)matrix_adjoint(wm_x_snk) * gamma5;

        for (int gram = 0; gram < 5; gram++)
        {
            const WilsonMatrix wm_sp = wm_x_src * block.block_pp[gram] * wm_snk_x;
            baryon_pi_data_temp[gram + 5] += coef * matrix_trace(wm_sp, gamma5);
        }
    }

    inline void contract_baryon_pi_unshifted_psel_typeI_core(LatData &baryon_pi_data_omp,
                                                             const std::string &job_tag, const int &traj,
                                                             const PselProp &prop_x_y, const Coordinate &xg_y,
                                                             const Coordinate &xg_x, const Coordinate &total_site,
                                                             const PselProp &prop_x_snk, const long &xg_x_idx,
                                                             const int &dt, const int &tsep, const int &t_baryon_src, const int &t_baryon_snk, const Baryon_Block_Psel_Psel &block1, const Complex coef)
    {
        const Coordinate xg_sep = smod(xg_x - xg_y, total_site);
        const WilsonMatrix &wm_x_y = prop_x_y.get_elem(xg_x_idx);
        const WilsonMatrix &wm_x_snk = prop_x_snk.get_elem(xg_x_idx);
        // std::vector<Vector<Complex>> vv(5);
        // for (int gram = 0; gram < (int)baryon_pi.size(); gram++)
        // {
        //     vv[gram] = baryon_pi[gram].get_elems(xg_x_idx);
        // }
        Vector<Complex> baryon_pi_data_temp = lat_data_complex_get(baryon_pi_data_omp, make_array(omp_get_thread_num(), dt));

        if (tsep >= 0)
        {
            if (t_baryon_src <= t_baryon_snk)
            {
                contract_baryon_pi_unshifted_psel_acc_x_typeI(baryon_pi_data_temp, coef, wm_x_snk, wm_x_y, block1);
                // baryon_pi_data_omp +=baryon_pi_data_temp;
            }
            else
            {
                const Complex coef1 = -1.0 * coef;
                contract_baryon_pi_unshifted_psel_acc_x_typeI(baryon_pi_data_temp, coef1, wm_x_snk, wm_x_y, block1);
            }
        }
        else
        {
            if (mod(t_baryon_src + tsep, total_site[3]) <= mod(t_baryon_snk - tsep, total_site[3]))
            {
                contract_baryon_pi_unshifted_psel_acc_x_typeI(baryon_pi_data_temp, coef, wm_x_snk, wm_x_y, block1);
                // baryon_pi_data_omp +=baryon_pi_data_temp;
            }
            else
            {
                const Complex coef1 = -1.0 * coef;
                contract_baryon_pi_unshifted_psel_acc_x_typeI(baryon_pi_data_temp, coef1, wm_x_snk, wm_x_y, block1);
            }
        }
    }

    inline void contract_baryon_pi_unshifted_psel_typeI(LatData &baryon_pi_data_omp,
                                                        const std::string &job_tag, const int &traj,
                                                        const std::vector<int> &psel_baryon_num_list,
                                                        const std::vector<int> &list_n_from_idx_baryon,
                                                        const PselProp &prop_x_y, const Coordinate &xg_y,
                                                        const long &xg_y_psel_idx, const int &tsep, const int &ti, const PointSelection &psel_baryon, const PointSelection &psel_meson, const FieldSelection &fsel,
                                                        const std::vector<std::vector<long>> baryon_pi_two_func_coeff, const int dtmin, const std::vector<std::vector<std::vector<Baryon_Block_Psel_Psel>>> &block,
                                                        const std::vector<std::vector<Baryon_Block_Psel_Psel>> &ppblock, const bool is_baryon_smear, const bool is_meson_smear)
    {
        TIMER_VERBOSE("contract_baryon_pi_unshifted_psel_typeI");
        const Geometry &geo = fsel.f_rank.geo();
        const Coordinate &total_site = geo.total_site();
        const int num_dtxy = get_baryon_pi_max_sep(total_site[3]);
        const int t_meson = xg_y[3];
        const long &psel_meson_num = psel_meson.size();
        const long &psel_baryon_num = psel_baryon.size();

        for (int dt = 0; dt < num_dtxy; ++dt)
        {
            const int dtxy = dt - num_dtxy / 2;
            const int total_time_length = abs(dtxy) + 2 * abs(tsep);
            
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
            const int dt_y = mod(t_meson - t_baryon_src, total_site[3]);
            const Complex coef = 1.0 / (double)baryon_pi_two_func_coeff[ti][dt];

            if (!(total_time_length == (tsep >= 0 ? abs(smod(t_baryon_snk - t_baryon_src, total_site[3])) : abs(smod(xt - t_meson, total_site[3])))))
            {
                displayln(ssprintf("thread_id: %d/%d, tsep: %d, t_baryon_snk: %d, xt: %d, t_meson: %d, t_baryon_src: %d,total_time_length: %d,,dtxy: %d,%d", omp_get_thread_num(), omp_get_num_threads(), tsep, t_baryon_snk, xt, t_meson, t_baryon_src, total_time_length, dtxy, (tsep >= 0 ? std::abs(smod(t_baryon_snk - t_baryon_src, total_site[3])) : abs(smod(xt - t_meson, total_site[3])))));
            }
            qassert(total_time_length == (tsep >= 0 ? abs(smod(t_baryon_snk - t_baryon_src, total_site[3])) : abs(smod(xt - t_meson, total_site[3]))));

            for (long nsnk = 0; nsnk < psel_baryon_num; ++nsnk)
            {
                const long xg_snk_psel_idx = nsnk;
                const Coordinate &xg_snk = psel_baryon[xg_snk_psel_idx];
                if (xg_snk[3] != t_baryon_snk)
                {
                    continue;
                }
                const PselProp &prop_x_snk = get_psel_prop(job_tag, traj, xg_snk, is_meson_smear, is_baryon_smear);
                int idx_snk = list_n_from_idx_baryon[xg_snk_psel_idx];
                const Baryon_Block_Psel_Psel &block1 = ppblock[dt][idx_snk];

                // idx_snk += 1;

                if (block1.xg_snk_psel_idx == -1)
                {
                    continue;
                }
                qassert(block1.is_build);
                if (!(block1.xg_snk_psel_idx == xg_snk_psel_idx))
                {
                    displayln(ssprintf("thread_id: %d/%d, tsep: %d, t_baryon_snk: %d, xt: %d, t_meson: %d, t_baryon_src: %d, xg_snk_psel_idx: %d, block1.xg_snk_psel_idx: %d, list_n_from_idx_baryon[xg_snk_psel_idx]: %d, list_n_from_idx_baryon[block1.xg_snk_psel_idx]: %d, idx_snk:%d, dt:%d", omp_get_thread_num(), omp_get_num_threads(), tsep, t_baryon_snk, xt, t_meson, t_baryon_src, xg_snk_psel_idx, block1.xg_snk_psel_idx, list_n_from_idx_baryon[xg_snk_psel_idx], list_n_from_idx_baryon[block1.xg_snk_psel_idx], idx_snk, dt));
                }
                qassert(block1.xg_snk_psel_idx == xg_snk_psel_idx);

#pragma omp parallel for
                for (long idx = 0; idx < psel_meson_num; ++idx)
                {
                    const long xg_x_idx = idx;
                    const Coordinate &xg_x = psel_meson[xg_x_idx];
                    // const long index = fsel.indices[xg_x_idx];
                    // const Coordinate xl = geo.coordinate_from_index(index);
                    // const Coordinate xg = geo.coordinate_g_from_l(xl);
                    // const Coordinate &xg_x = xg;
                    if (xg_x[3] != xt)
                    {
                        continue;
                    }
                    if (tsep == 0 && is_baryon_smear == is_meson_smear)
                    {
                        if (dtxy >= 0)
                        {
                            if (xg_snk_psel_idx == xg_x_idx)
                            {
                                continue;
                            }
                            contract_baryon_pi_unshifted_psel_typeI_core(baryon_pi_data_omp,
                                                                         job_tag, traj,
                                                                         prop_x_y, xg_y,
                                                                         xg_x, total_site,
                                                                         prop_x_snk, xg_x_idx,
                                                                         dt, tsep, t_baryon_src, t_baryon_snk, block1, coef);
                        }
                        else
                        {
                            // if (block1.xg_src_psel_idx == xg_x_idx)
                            // {
                            //     continue;
                            // }
                            Baryon_Block_Psel_Psel ppblock_cor, block2;
                            qassert(psel_baryon == psel_meson);
                            contract_baryon_pselpsel_y_block_coincide_corr(block, ppblock_cor, job_tag, traj, geo, psel_baryon, psel_baryon_num_list, list_n_from_idx_baryon, xg_snk_psel_idx, xg_x_idx, xg_y_psel_idx, xg_y, tsep, dtxy, dtmin, is_baryon_smear, is_meson_smear);

                            qassert(ppblock_cor.xg_snk_psel_idx == block1.xg_snk_psel_idx);
                            qassert(ppblock_cor.n_snk == block1.n_snk);
                            qassert(ppblock_cor.xg_src_psel_idx == xg_y_psel_idx);
                            qassert(block1.xg_src_psel_idx == xg_y_psel_idx);
                            qassert(ppblock_cor.n_src == block1.n_src);
                            qassert(ppblock_cor.is_build == true);
                            qassert(block1.is_build == true);

                            for (int gram = 0; gram < (int)block1.block_pp.size(); ++gram)
                            {
                                block2.block_pp[gram] = block1.block_pp[gram] - ppblock_cor.block_pp[gram];
                            }

                            block2.xg_snk_psel_idx = block1.xg_snk_psel_idx;
                            block2.n_snk = block1.n_snk;
                            block2.xg_src_psel_idx = block1.xg_src_psel_idx;
                            block2.n_src = block1.n_src;
                            block2.is_build = true;

                            contract_baryon_pi_unshifted_psel_typeI_core(baryon_pi_data_omp,
                                                                         job_tag, traj,
                                                                         prop_x_y, xg_y,
                                                                         xg_x, total_site,
                                                                         prop_x_snk, xg_x_idx,
                                                                         dt, tsep, t_baryon_src, t_baryon_snk, block2, coef);
                        }
                    }
                    else
                    {
                        contract_baryon_pi_unshifted_psel_typeI_core(baryon_pi_data_omp,
                                                                     job_tag, traj,
                                                                     prop_x_y, xg_y,
                                                                     xg_x, total_site,
                                                                     prop_x_snk, xg_x_idx,
                                                                     dt, tsep, t_baryon_src, t_baryon_snk, block1, coef);
                    }
                }
            }
        }
        return;
    }

    inline void contract_baryon_pi_unshifted_psel_typeII(LatData &baryon_pi_data_omp,
                                                         const std::string &job_tag, const int &traj,
                                                         const std::vector<int> &psel_baryon_num_list,
                                                         const std::vector<int> &list_n_from_idx_baryon,
                                                         const Coordinate &xg_y,
                                                         const long &xg_y_psel_idx, const int &tsep, const int &ti, const PointSelection &psel_baryon, const PointSelection &psel_meson, const FieldSelection &fsel,
                                                         const std::vector<std::vector<long>> baryon_pi_two_func_coeff,
                                                         const int dtmin,
                                                         const std::vector<std::vector<Baryon_Block_Psel_Psel>> &spblock, const bool is_baryon_smear, const bool is_meson_smear)
    {
        TIMER_VERBOSE("contract_baryon_pi_unshifted_psel_typeII");
        const Geometry &geo = fsel.f_rank.geo();
        const Coordinate &total_site = geo.total_site();
        const int num_dtxy = get_baryon_pi_max_sep(total_site[3]);
        const int t_meson = xg_y[3];
        const long &psel_meson_num = psel_meson.size();
        const long &psel_baryon_num = psel_baryon.size();

        // const int multiplicity = 8 * 8;
        // clear(baryon_pi);
        // baryon_pi.resize(5);

        const int x_sparsity = 1;

        // for (int gram = 0; gram < (int)baryon_pi.size(); gram++)
        // {
        //     baryon_pi[gram].init(fsel, multiplicity);
        //     set_zero(baryon_pi[gram]);
        // }
        for (int dt = 0; dt < num_dtxy; ++dt)
        {
            if ((int)spblock[dt].size() == 0)
            {
                continue;
            }
            const int dtxy = dt - num_dtxy / 2;
            const int total_time_length = abs(dtxy) + 2 * abs(tsep);
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

            qassert(total_time_length == (tsep >= 0 ? abs(smod(t_baryon_snk - t_baryon_src, total_site[3])) : abs(smod(xt - t_meson, total_site[3]))));

            const int dt_y = mod(t_meson - t_baryon_src, total_site[3]);
            const Complex coef = 1.0 / (double)baryon_pi_two_func_coeff[ti][dt] * (double)x_sparsity;
            const int num_src = psel_baryon_num_list[t_baryon_src];

            // int idx_snk = 0;

            for (long nsnk = 0; nsnk < psel_baryon_num; ++nsnk)
            {
                const long xg_snk_psel_idx = nsnk;
                const Coordinate &xg_snk = psel_baryon[xg_snk_psel_idx];
                int idx_snk = list_n_from_idx_baryon[xg_snk_psel_idx];
                if (xg_snk[3] != t_baryon_snk)
                {
                    continue;
                }
                const PselProp &prop_x_snk = get_psel_prop(job_tag, traj, xg_snk, is_meson_smear, is_baryon_smear);
                // int idx_src = 0;
                for (long nsrc = 0; nsrc < psel_baryon_num; ++nsrc)
                {
                    const long xg_src_psel_idx = nsrc;
                    const Coordinate &xg_src = psel_baryon[xg_src_psel_idx];
                    if (xg_src[3] != t_baryon_src)
                    {
                        continue;
                    }
                    const PselProp &prop_x_src = get_psel_prop(job_tag, traj, xg_src, is_meson_smear, is_baryon_smear);

                    int idx_src = list_n_from_idx_baryon[xg_src_psel_idx];
                    if (tsep == 0 && is_baryon_smear == is_meson_smear)
                    {
                        if (dtxy >= 0)
                        {
                            if (xg_src_psel_idx == xg_y_psel_idx)
                            {
                                continue;
                            }
                        }
                        else
                        {
                            if (xg_snk_psel_idx == xg_y_psel_idx)
                            {
                                continue;
                            }
                        }
                    }

                    const int idx_src_snk = idx_src + idx_snk * num_src;
                    const Baryon_Block_Psel_Psel &block1 = spblock[dt][idx_src_snk];
                    qassert(block1.is_build);
                    qassert(block1.xg_snk_psel_idx == xg_snk_psel_idx);
                    if (!(block1.xg_src_psel_idx == xg_src_psel_idx))
                    {
                        displayln(ssprintf("thread_id: %d/%d, tsep: %d, t_baryon_snk: %d, xt: %d, t_meson: %d, t_baryon_src: %d, xg_snk_psel_idx: %d, block1.xg_snk_psel_idx: %d, xg_src_psel_idx: %d, block1.xg_src_psel_idx: %d, idx_snk: %d, idx_snk: %d, idx_src: %d, idx_src: %d, dt: %d", omp_get_thread_num(), omp_get_num_threads(), tsep, t_baryon_snk, xt, t_meson, t_baryon_src, xg_snk_psel_idx, block1.xg_snk_psel_idx, xg_src_psel_idx, block1.xg_src_psel_idx,
                                           idx_snk, list_n_from_idx_baryon[block1.xg_snk_psel_idx],
                                           idx_src, list_n_from_idx_baryon[block1.xg_src_psel_idx], dt));
                    }
                    qassert(block1.xg_src_psel_idx == xg_src_psel_idx);

#pragma omp parallel for
                    for (long idx = 0; idx < psel_meson_num; ++idx)
                    {
                        if (mod(idx, x_sparsity) != 0)
                        {
                            continue;
                        }
                        const long xg_x_idx = idx;
                        // const long index = fsel.indices[xg_x_idx];
                        // const Coordinate xl = geo.coordinate_from_index(index);
                        // const Coordinate xg = geo.coordinate_g_from_l(xl);
                        const Coordinate &xg_x = psel_meson[xg_x_idx];
                        if (xg_x[3] != xt)
                        {
                            continue;
                        }
                        const Coordinate xg_sep = smod(xg_x - xg_y, total_site);

                        if (tsep == 0 && is_baryon_smear == is_meson_smear)
                        {
                            if (dtxy >= 0)
                            {
                                if (xg_snk_psel_idx == xg_x_idx)
                                {
                                    continue;
                                }
                            }
                            else
                            {
                                if (block1.xg_src_psel_idx == xg_x_idx)
                                {
                                    continue;
                                }
                            }
                        }
                        const WilsonMatrix &wm_x_snk = prop_x_snk.get_elem(xg_x_idx);
                        const WilsonMatrix &wm_x_src = prop_x_src.get_elem(xg_x_idx);
                        // std::vector<Vector<Complex>> vv(5);
                        // for (int gram = 0; gram < (int)baryon_pi.size(); gram++)
                        // {
                        //     vv[gram] = baryon_pi[gram].get_elems(xg_x_idx);
                        // }
                        Vector<Complex> baryon_pi_data_temp = lat_data_complex_get(baryon_pi_data_omp, make_array(omp_get_thread_num(), dt));

                        if (tsep >= 0)
                        {
                            if (t_baryon_src <= t_baryon_snk)
                            {
                                contract_baryon_pi_unshifted_psel_acc_x_typeII(baryon_pi_data_temp, coef, wm_x_snk, wm_x_src, block1);
                            }
                            else
                            {
                                const Complex coef1 = -1.0 * coef;
                                contract_baryon_pi_unshifted_psel_acc_x_typeII(baryon_pi_data_temp, coef1, wm_x_snk, wm_x_src, block1);
                            }
                        }
                        else
                        {
                            if (mod(t_baryon_src + tsep, total_site[3]) <= mod(t_baryon_snk - tsep, total_site[3]))
                            {
                                contract_baryon_pi_unshifted_psel_acc_x_typeII(baryon_pi_data_temp, coef, wm_x_snk, wm_x_src, block1);
                            }
                            else
                            {
                                const Complex coef1 = -1.0 * coef;
                                contract_baryon_pi_unshifted_psel_acc_x_typeII(baryon_pi_data_temp, coef1, wm_x_snk, wm_x_src, block1);
                            }
                        }
                    }
                    // idx_src += 1;
                }
                // idx_snk += 1;
            }
        }
        return;
    }

    inline void
    contract_baryon_pi_two_point_all_psel_acc_psel_typeI(
        LatData &baryon_pi_data_omp,
        const std::string &job_tag, const int &traj,
        const std::vector<int> &psel_baryon_num_list,
        const std::vector<int> &list_n_from_idx_baryon,
        const PselProp &prop_x_y, const Coordinate &xg_y,
        const long &xg_y_psel_idx, const int &tsep, const int &ti,
        const PointSelection &psel_baryon, const PointSelection &psel_meson,
        const FieldSelection &fsel,
        const std::vector<std::vector<long>> &baryon_pi_two_func_coeff,
        const int dtmin, const std::vector<std::vector<std::vector<Baryon_Block_Psel_Psel>>> &block,
        const std::vector<std::vector<Baryon_Block_Psel_Psel>> &ppblock, const bool is_baryon_smear, const bool is_meson_smear)
    {
        TIMER_VERBOSE("contract_baryon_pi_two_point_all_psel_acc_psel");
        const Geometry &geo = fsel.f_rank.geo();
        qassert(psel_meson[xg_y_psel_idx] == xg_y);
        contract_baryon_pi_unshifted_psel_typeI(baryon_pi_data_omp, job_tag, traj, psel_baryon_num_list, list_n_from_idx_baryon, prop_x_y, xg_y, xg_y_psel_idx,
                                                tsep, ti, psel_baryon, psel_meson, fsel, baryon_pi_two_func_coeff, dtmin, block, ppblock, is_baryon_smear, is_meson_smear);
    }

    inline void
    contract_baryon_pi_two_point_all_psel_acc_psel_typeII(
        LatData &baryon_pi_data_omp,
        const std::string &job_tag, const int &traj,
        const std::vector<int> &psel_baryon_num_list,
        const std::vector<int> &list_n_from_idx_baryon,
        const PselProp &prop_x_y, const Coordinate &xg_y,
        const long &xg_y_psel_idx, const int &tsep, const int &ti, const PointSelection &psel_baryon, const PointSelection &psel_meson,
        const FieldSelection &fsel,
        const std::vector<std::vector<long>> &baryon_pi_two_func_coeff,
        const int dtmin,
        const std::vector<std::vector<Baryon_Block_Psel_Psel>> &spblock, const bool is_baryon_smear, const bool is_meson_smear)
    {
        TIMER_VERBOSE("contract_baryon_pi_two_point_all_psel_acc_psel");
        const Geometry &geo = fsel.f_rank.geo();
        qassert(psel_meson[xg_y_psel_idx] == xg_y);
        contract_baryon_pi_unshifted_psel_typeII(baryon_pi_data_omp, job_tag, traj, psel_baryon_num_list, list_n_from_idx_baryon, xg_y, xg_y_psel_idx, tsep, ti, psel_baryon, psel_meson,
                                                 fsel, baryon_pi_two_func_coeff, dtmin, spblock, is_baryon_smear, is_meson_smear);
    }

    inline std::string get_baryon_pi_two_point_all_psel_func_psel_path(const std::string &job_tag,
                                                                       const int &traj)
    {
        return ssprintf("analysis/baryon_pi_two_point/%s/results=%d", job_tag.c_str(),
                        traj);
    }

    inline LatData mk_baryon_pi_two_point(const std::vector<int> &gram_convert_list, const Coordinate &total_site)
    {
        LatData ld;
        ld.info.push_back(lat_dim_number("txy", 0, get_baryon_pi_max_sep(total_site[3]) - 1));
        ld.info.push_back(lat_dim_number("gram", 0, gram_convert_list.size() - 1));
        // ld.info.push_back(lat_dim_number("mu", 0, 8 - 1));
        // ld.info.push_back(lat_dim_number("nu", 0, 8 - 1));
        ld.info.push_back(lat_dim_re_im());
        lat_data_alloc(ld);
        set_zero(ld);
        return ld;
    }

    inline LatData mk_baryon_pi_two_point_for_omp(const std::vector<int> &gram_convert_list, const Coordinate &total_site)
    {
        LatData ld;
        ld.info.push_back(lat_dim_number("omp", 0, omp_get_max_threads() - 1));
        ld.info.push_back(lat_dim_number("txy", 0, get_baryon_pi_max_sep(total_site[3]) - 1));
        ld.info.push_back(lat_dim_number("gram", 0, gram_convert_list.size() - 1));
        // ld.info.push_back(lat_dim_number("mu", 0, 8 - 1));
        // ld.info.push_back(lat_dim_number("nu", 0, 8 - 1));
        ld.info.push_back(lat_dim_re_im());
        lat_data_alloc(ld);
        set_zero(ld);
        return ld;
    }

    inline void compute_baryon_pi_two_point_all_psel_func_psel_type(const std::string &job_tag, const int &traj, const std::vector<int> &gram_convert_list, const std::vector<int> &tsep_list, const bool is_baryon_smear, const bool is_meson_smear)
    {
        check_sigterm();
        check_time_limit();
        Timer::autodisplay();
        const int gram_num = 10;
        const int tsep_num = (int)tsep_list.size();
        qassert(gram_convert_list.size() == gram_num);
        const std::string path = get_baryon_pi_two_point_all_psel_func_psel_path(job_tag, traj);
        std::vector<std::string> fn_gram_list(tsep_num);
        // std::vector<std::string> fn_gram_list_omp(tsep_num);
        for (int ti = 0; ti < tsep_num; ++ti)
        {
            const int &tsep = tsep_list[ti];
            fn_gram_list[ti] = path + ssprintf("/baryon_%s_pi_%s_two_point_dt%d.field", (is_baryon_smear ? "smear" : "point"), (is_meson_smear ? "smear" : "point"), tsep);
            // fn_gram_list_omp[ti] = path + ssprintf("/baryon_pi_two_point_omp_dt%d.field", tsep);
        }

        bool is_complete = true;
        for (int ti = 0; ti < tsep_num; ++ti)
        {
            if (fn_gram_list[ti] != "")
            {
                if (not is_d_field(fn_gram_list[ti]))
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
        TIMER_VERBOSE("compute_baryon_pi_two_point_all_psel_func_psel_type");
        const PointSelection &psel = get_point_selection(job_tag, traj);
        const long n_points = psel.size();
        const PointSelection &psel_smear = get_point_selection_smear(job_tag, traj);
        const long n_points_smear = psel_smear.size();

        const PointSelection &psel_meson = (!is_meson_smear ? psel : psel_smear);
        const long &psel_meson_num = (!is_meson_smear ? n_points : n_points_smear);
        const PointSelection &psel_baryon = (!is_baryon_smear ? psel : psel_smear);
        const long &psel_baryon_num = (!is_baryon_smear ? n_points : n_points_smear);

        const FieldSelection &fsel = get_field_selection(job_tag, traj);
        // std::vector<std::vector<FieldM<Complex, 8 * 8>>> baryon_pi_data_omp(tsep_num);
        std::vector<LatData> baryon_pi_data_omp(tsep_num);

        const Geometry &geo = fsel.f_rank.geo();
        const Coordinate total_site = geo.total_site();
        const int num_dtxy = get_baryon_pi_max_sep(total_site[3]);
        const int dtmax = num_dtxy / 2 + 2 * (tsep_list[tsep_num - 1] > 0 ? tsep_list[tsep_num - 1] : 0);
        const int dtmin = 2 * (tsep_list[0] > 0 ? tsep_list[0] : 0); // block 所需时间片长度

        for (int ti = 0; ti < tsep_num; ++ti)
        {
            baryon_pi_data_omp[ti] = mk_baryon_pi_two_point_for_omp(gram_convert_list, total_site);
        }

        std::vector<std::vector<long>> baryon_pi_two_func_coeff;
        std::vector<int> psel_baryon_num_list;
        std::vector<int> psel_meson_num_list;
        std::vector<int> list_n_from_idx_baryon(psel_meson_num);
        std::vector<int> list_n_from_idx_meson(psel_meson_num);

        std::vector<std::vector<std::vector<Baryon_Block_Psel_Psel>>> block;

        baryon_pion_two_point_coeffience(
            job_tag, traj, baryon_pi_two_func_coeff,
            psel_baryon_num_list, psel_meson_num_list,
            list_n_from_idx_baryon, list_n_from_idx_meson,
            geo, psel_baryon, psel_meson, tsep_list, is_baryon_smear, is_meson_smear);

        displayln_info("baryon_pion_two_point_coeffience_finish;");

        contract_baryon_pselpsel_block(block, psel_baryon_num_list, list_n_from_idx_baryon, job_tag, traj, geo, psel_baryon, dtmax, dtmin, is_baryon_smear);

        displayln_info("contract_baryon_pselpsel_block_finish;");

        const int y_idx_sparsity = 1;
        long iter = 0;
        // #pragma omp parallel for
        for (int t_meson = 0; t_meson < total_site[3]; ++t_meson)
        {
            // load_prop_psrc_baryon_pi(job_tag, traj, tsep_list.back(), num_dtxy, t_meson, psel, geo);
            for (long ny = 0; ny < psel_meson_num; ++ny)
            {
                const long xg_y_psel_idx = ny;
                const Coordinate &xg_y = psel_meson[xg_y_psel_idx];
                if (t_meson != xg_y[3])
                {
                    continue;
                }
                if (mod((int)xg_y_psel_idx, y_idx_sparsity) != 0)
                {
                    continue;
                }
                // if (get_point_src_info(job_tag, traj, xg_y, 0).size() == 0)
                // {
                //     continue;
                // }
                Timer::autodisplay();
                TIMER_VERBOSE("compute_baryon_pi_two_point_all_psel_func_noama_type-iter");

                iter += 1;
                // const SelProp &prop_x_y = get_prop_psrc(job_tag, traj, xg_y, 0, 0);
                // const ShiftShufflePlan ssp = make_shift_shuffle_plan(fsel, -xg_y);
                const PselProp &prop_x_y = get_psel_prop(job_tag, traj, xg_y, is_meson_smear, is_meson_smear);

                displayln_info(fname + ssprintf(":n=%ld iter=%ld", ny, iter));

                for (int ti = 0; ti < tsep_num; ++ti)
                {
                    const int tsep = tsep_list[ti];
                    std::vector<std::vector<Baryon_Block_Psel_Psel>> ppblock(num_dtxy);
                    std::vector<std::vector<Baryon_Block_Psel_Psel>> spblock(num_dtxy);
                    contract_baryon_pselpsel_y_block(block, ppblock, psel_baryon_num_list, list_n_from_idx_baryon, job_tag, traj, geo, psel_baryon, psel_meson, xg_y_psel_idx, xg_y, dtmin, tsep, is_baryon_smear, is_meson_smear);

                    // displayln_info(fname + ssprintf(":n=%ld iter=%ld", ny, iter) + "y_block");

                    contract_baryon_pi_two_point_all_psel_acc_psel_typeI(baryon_pi_data_omp[ti], job_tag, traj, psel_baryon_num_list, list_n_from_idx_baryon,
                                                                         prop_x_y, xg_y, xg_y_psel_idx, tsep, ti,
                                                                         psel_baryon, psel_meson, fsel, baryon_pi_two_func_coeff, dtmin, block, ppblock, is_baryon_smear, is_meson_smear);
                    // displayln_info(fname + ssprintf(":n=%ld iter=%ld", ny, iter) + "all_psel_acc_psel_typeI ");

                    contract_baryon_sequential_meson_y_block(spblock, psel_baryon_num_list, list_n_from_idx_baryon, job_tag, traj, geo, psel_baryon, psel_meson,
                                                             xg_y_psel_idx, xg_y, dtmin, tsep, is_baryon_smear, is_meson_smear);

                    // displayln_info(fname + ssprintf(":n=%ld iter=%ld", ny, iter) + "sp_y_block ");

                    contract_baryon_pi_two_point_all_psel_acc_psel_typeII(baryon_pi_data_omp[ti], job_tag, traj, psel_baryon_num_list, list_n_from_idx_baryon, prop_x_y, xg_y,
                                                                          xg_y_psel_idx, tsep, ti, psel_baryon, psel_meson, fsel,
                                                                          baryon_pi_two_func_coeff, dtmin, spblock, is_baryon_smear, is_meson_smear);

                    // displayln_info(fname + ssprintf(":n=%ld iter=%ld", ny, iter) + "all_psel_acc_psel_typeII");
                }
            }
        }

        for (int ti = 0; ti < tsep_num; ++ti)
        {
            LatData baryon_pi_data = mk_baryon_pi_two_point(gram_convert_list, total_site);

            // lat_data_save_info(fn_gram_list_omp[ti], baryon_pi_data_omp[ti]);
            for (int txy = 0; txy < get_baryon_pi_max_sep(total_site[3]); txy++)
            {
                Vector<Complex> temp_vector = lat_data_complex_get(baryon_pi_data, make_array(txy));
                for (int thread_id = 0; thread_id < omp_get_max_threads(); thread_id++)
                {
                    for (int gram = 0; gram < gram_convert_list.size(); gram++)
                    {
                        temp_vector[gram_convert_list[gram] - 1] += lat_data_complex_get(baryon_pi_data_omp[ti], make_array(thread_id, txy, gram))[0] * (double)y_idx_sparsity;
                    }
                }
            }
            lat_data_save_info(fn_gram_list[ti], baryon_pi_data);
        }
    }

    inline void compute_baryon_pi_two_point_all_psel_func_psel(const std::string &job_tag,
                                                               const int &traj, const bool is_baryon_smear, const bool is_meson_smear)
    {
        check_sigterm();
        Timer::autodisplay();
        TIMER_VERBOSE("compute_baryon_pi_two_point_all_psel_func_psel");
        const std::string path = get_baryon_pi_two_point_all_psel_func_psel_path(job_tag, traj);
        const std::string path_checkpoint = path + "/checkpoint.txt";
        if (does_file_exist_sync_node(path_checkpoint))
        {
            return;
        }
        if (does_file_exist_sync_node(path + ssprintf("/baryon_%s_pi_%s_two_point_dt%d.field", (is_baryon_smear ? "smear" : "point"), (is_meson_smear ? "smear" : "point"), 0)))
        {
            return;
        }
        if (not(check_prop_psrc(job_tag, traj, 0)))
        {
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
        qmkdir_info("analysis/baryon_pi_two_point");
        qmkdir_info(ssprintf("analysis/baryon_pi_two_point/%s", job_tag.c_str()));
        qmkdir_info(path);
        std::vector<int> gram_convert_list;
        gram_convert_list.push_back(3);
        gram_convert_list.push_back(4);
        gram_convert_list.push_back(8);
        gram_convert_list.push_back(9);
        gram_convert_list.push_back(10);
        gram_convert_list.push_back(1);
        gram_convert_list.push_back(2);
        gram_convert_list.push_back(5);
        gram_convert_list.push_back(6);
        gram_convert_list.push_back(7);

        std::vector<int> tsep_list;
        // tsep_list.push_back(-1);
        tsep_list.push_back(0);
        // tsep_list.push_back(1);
        
        compute_baryon_pi_two_point_all_psel_func_psel_type(job_tag, traj, gram_convert_list, tsep_list, is_baryon_smear, is_meson_smear);

        release_lock();
    }
}
