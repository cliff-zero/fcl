#pragma once

#include "data-load.h"
#include "contract_proton.h"
#include "compute-utils.h"

namespace qlat
{
#ifndef PROTONPPBLOCK
#define PROTONPPBLOCK
    struct ProtonPPBlock
    {
        bool is_build;
        array<WilsonMatrix, 5> block_pp;
        long xg_snk_psel_idx;
        long xg_src_psel_idx;
        int n_snk;
        int n_src;

        ProtonPPBlock()
        {
            is_build = false;
            xg_snk_psel_idx = -1;
            xg_src_psel_idx = -1;
            n_snk = -1;
            n_src = -1;
            for (int i = 0; i < 5; i++)
            {
                set_zero(block_pp[i]);
            }
        }
    };
    struct ProtonSPBlock
    {
        bool is_build;
        array<WilsonMatrix, 8 * 5> block_sp;
        long xg_snk_psel_idx;
        long xg_src_psel_idx;
        int n_snk;
        int n_src;

        ProtonSPBlock()
        {
            is_build = false;
            xg_snk_psel_idx = -1;
            xg_src_psel_idx = -1;
            n_snk = -1;
            n_src = -1;
            for (int i = 0; i < 5; ++i)
            {
                for (int j = 0; j < 8; ++j)
                {
                    set_zero(block_sp[j + 5 * i]);
                }
            }
        }
    };
#endif
    inline void load_prop_smear_fourp(const std::string &job_tag, const int traj, const int tsep,
                                      const int num_dtxy, const int yt, const PointSelection &psel,
                                      const Geometry &geo)
    {
        TIMER_VERBOSE("load_prop_smear_fourp");
        const long n_points = psel.size();
        const Coordinate total_site = geo.total_site();
        const int num_t = 2 * tsep + num_dtxy;
        // const int num_t = total_site[3];
        const int t_begin = mod(yt - num_t / 2, total_site[3]);
        for (int ti = 0; ti < num_t; ++ti)
        {
            const int tsrc = mod(t_begin + ti, total_site[3]);
            for (long nsrc = 0; nsrc < n_points; ++nsrc)
            {
                const long xg_src_psel_idx = nsrc;
                const Coordinate &xg_src = psel[xg_src_psel_idx];
                if (xg_src[3] != tsrc)
                {
                    continue;
                }
                const SelProp &prop = get_prop_smear(job_tag, traj, xg_src, 0, 0, 0);
                qassert(prop.initialized);
            }
        }
        for (long nsrc = 0; nsrc < n_points; ++nsrc)
        {
            const long xg_src_psel_idx = nsrc;
            const Coordinate xg_src = psel[xg_src_psel_idx];
            if (xg_src[3] != t_begin)
            {
                continue;
            }
            const SelProp &prop = get_prop_smear(job_tag, traj, xg_src, 0, 0, 0);
            qassert(prop.initialized);
            break;
        }
    }

    inline void contract_proton_pp_block_smear(std::vector<std::vector<std::vector<ProtonPPBlock>>> &block,
                                               std::vector<int> &psel_num_list_smear,
                                               std::vector<std::vector<double>> &psel_dt_num_list,
                                               const std::string &job_tag, const int traj,
                                               const Geometry &geo, const PointSelection &psel,
                                               const PointSelection &psel_smear,
                                               const int dtmax, const int dtmin)
    {
        TIMER_VERBOSE("contract_proton_pp_block_smear");
        const int num_dt = dtmax - dtmin + 1;
        qassert(num_dt > 0);
        const long n_points = psel.size();
        const long n_points_smear = psel_smear.size();
        const Coordinate total_site = geo.total_site();

        block.resize(n_points_smear);
        for (long n = 0; n < n_points_smear; ++n)
        {
            block[n].resize(num_dt);
        }

        psel_num_list_smear.resize(total_site[3]);
        set_zero(psel_num_list_smear);
        psel_dt_num_list.resize(num_dt);
        for (int dt = 0; dt < num_dt; ++dt)
        {
            const int del_t = dtmin + dt;
            psel_dt_num_list[dt].resize(del_t + 1);
            set_zero(psel_dt_num_list[dt]);
        }
        std::vector<int> list_n_from_idx(n_points);
        for (long n = 0; n < n_points_smear; ++n)
        {
            const long xg_psel_idx = n;
            const Coordinate xg = psel_smear[xg_psel_idx];
            const int t = xg[3];
            list_n_from_idx[xg_psel_idx] = psel_num_list_smear[t];
            psel_num_list_smear[t] += 1;
        }
        for (long nsnk = 0; nsnk < n_points_smear; ++nsnk)
        {
            const long xg_snk_psel_idx = nsnk;
            const Coordinate &xg_snk = psel_smear[xg_snk_psel_idx];
            const int tsnk = xg_snk[3];
            for (int dt = 0; dt < num_dt; ++dt)
            {
                const int del_t = dtmin + dt;
                const int tsrc = mod(tsnk - del_t, total_site[3]);
                block[xg_snk_psel_idx][dt].resize(psel_num_list_smear[tsrc]);
            }
        }

        for (long nsrc = 0; nsrc < n_points_smear; ++nsrc)
        {
            const long xg_src_psel_idx = nsrc;
            const Coordinate &xg_src = psel_smear[xg_src_psel_idx];
            const int tsrc = xg_src[3];
            for (long nsnk = 0; nsnk < n_points_smear; ++nsnk)
            {
                const long xg_snk_psel_idx = nsnk;
                const Coordinate &xg_snk = psel_smear[xg_snk_psel_idx];
                const int tsnk = xg_snk[3];
                const int del_t = mod(tsnk - tsrc, total_site[3]);
                const int dt = del_t - dtmin;
                if ((dt < 0) || (dt >= num_dt))
                {
                    continue;
                }
                for (long ny = 0; ny < n_points; ++ny)
                {
                    const long xg_y_psel_idx = ny;
                    const Coordinate &xg_y = psel[xg_y_psel_idx];
                    const int ty = xg_y[3];
                    const int dt_y_src = mod(ty - tsrc, total_site[3]);
                    if ((dt_y_src < 0) || (dt_y_src > del_t))
                    {
                        continue;
                    }
                    psel_dt_num_list[dt][dt_y_src] += 1.0;
                }
            }
        }

        for (long nsrc = 0; nsrc < n_points_smear; ++nsrc)
        {
            const long xg_src_psel_idx = nsrc;
            const Coordinate &xg_src = psel_smear[xg_src_psel_idx];
            const int tsrc = xg_src[3];
            const PselProp &ssprop_src = get_psel_prop_smear(job_tag, traj, xg_src, 0, 0, 1);

#pragma omp parallel for
            for (long nsnk = 0; nsnk < n_points_smear; ++nsnk)
            {
                const long xg_snk_psel_idx = nsnk;
                const Coordinate &xg_snk = psel_smear[xg_snk_psel_idx];
                const int tsnk = xg_snk[3];
                const int del_t = mod(tsnk - tsrc, total_site[3]);
                const int dt = del_t - dtmin;
                if ((dt < 0) || (dt >= num_dt))
                {
                    continue;
                }
                ProtonPPBlock ppblock;
                const WilsonMatrix &wm = ssprop_src.get_elem(xg_snk_psel_idx);
                proton_fourp_block_pp(ppblock.block_pp, wm, wm);
                ppblock.xg_snk_psel_idx = xg_snk_psel_idx;
                ppblock.xg_src_psel_idx = xg_src_psel_idx;
                ppblock.n_snk = list_n_from_idx[xg_snk_psel_idx];
                ppblock.n_src = list_n_from_idx[xg_src_psel_idx];
                ppblock.is_build = true;

                const int n_t = list_n_from_idx[xg_src_psel_idx];
                block[xg_snk_psel_idx][dt][n_t] = ppblock;
            }
        }
        return;
    }

    inline void contract_proton_pp_y_block_smear(const std::vector<std::vector<std::vector<ProtonPPBlock>>> &block,
                                                 std::vector<std::vector<ProtonPPBlock>> &ppblock,
                                                 const std::vector<int> &psel_num_list_smear,
                                                 const std::string &job_tag, const int traj,
                                                 const Geometry &geo, const PointSelection &psel_smear,
                                                 const long &xg_y_psel_idx, const Coordinate &xg_y,
                                                 const int dtmin, const int tsep_src, const int tsep_snk)
    {
        TIMER_VERBOSE("contract_proton_pp_y_block_smear");
        const Coordinate &total_site = geo.total_site();
        const int num_dtxy = get_proton_fourp_max_sep(total_site[3]);
        qassert(num_dtxy == (int)ppblock.size());
        const int yt = xg_y[3];
        const long n_points = psel_smear.size();
        for (int dt = 0; dt < num_dtxy; dt++)
        {
            const int dtxy = dt - num_dtxy / 2;
            int xt = mod(yt + dtxy, total_site[3]);
            int tsrc, tsnk;
            if (dtxy >= 0)
            {
                tsrc = mod(yt - tsep_src, total_site[3]);
                tsnk = mod(xt + tsep_snk, total_site[3]);
            }
            else
            {
                tsrc = mod(xt - tsep_src, total_site[3]);
                tsnk = mod(yt + tsep_snk, total_site[3]);
            }
            const int num_psel = psel_num_list_smear[tsrc];
            const int num_psel_snk = psel_num_list_smear[tsnk];
            ppblock[dt].resize(num_psel_snk);

            for (long nsnk = 0; nsnk < n_points; ++nsnk)
            {
                const long xg_snk_psel_idx = nsnk;
                const Coordinate &xg_snk = psel_smear[xg_snk_psel_idx];
                if (xg_snk[3] != tsnk)
                {
                    continue;
                }
                const int dt_idx = mod(tsnk - tsrc, total_site[3]) - dtmin;
                for (int n_t = 0; n_t < num_psel; ++n_t)
                {
                    const ProtonPPBlock &block1 = block[xg_snk_psel_idx][dt_idx][n_t];
                    qassert(block1.is_build);
                    qassert(block1.xg_snk_psel_idx == xg_snk_psel_idx);
                    const int i_snk = block1.n_snk;
                    const long xg_src_psel_idx = block1.xg_src_psel_idx;
                    const Coordinate &xg_src = psel_smear[xg_src_psel_idx];
                    const PselProp &psprop = get_psel_prop_smear(job_tag, traj, xg_src, 0, 0, 0);
                    const WilsonMatrix &wm = psprop.get_elem(xg_y_psel_idx);
                    for (int gram = 0; gram < (int)block1.block_pp.size(); ++gram)
                    {
                        ppblock[dt][i_snk].block_pp[gram] += (wm * block1.block_pp[gram]);
                    }
                    ppblock[dt][i_snk].xg_snk_psel_idx = xg_snk_psel_idx;
                    ppblock[dt][i_snk].n_snk = i_snk;
                    ppblock[dt][i_snk].is_build = true;
                }
            }
        }
    }

    inline void get_all_wilson_matrix_smear(std::vector<WilsonMatrix> &wm_snk_src_all,
                                            std::vector<WilsonMatrix> &wm_y_src_all,
                                            std::vector<WilsonMatrix> &wm_snk_y_all,
                                            const std::string &job_tag, const int traj,
                                            const long &xg_y_psel_idx,
                                            const int tsrc, const int tsnk, const PointSelection &psel_smear,
                                            const int num_src, const int num_snk,
                                            std::vector<long> &list_idx_from_isrc,
                                            std::vector<long> &list_idx_from_isnk)
    {
        if (num_src == 0)
        {
            return;
        }
        if (num_snk == 0)
        {
            return;
        }

        const SpinMatrix &gamma5 = SpinMatrixConstants::get_gamma5();
        const long n_points = psel_smear.size();
        int i_snk = 0;
        int i_src = 0;
        wm_snk_src_all.resize(num_src * num_snk);
        wm_snk_y_all.resize(num_snk);
        wm_y_src_all.resize(num_src);
        list_idx_from_isnk.resize(num_snk);
        list_idx_from_isrc.resize(num_src);
        for (long n = 0; n < n_points; ++n)
        {
            const long xg_snk_psel_idx = n;
            const Coordinate &xg_snk = psel_smear[xg_snk_psel_idx];
            if (xg_snk[3] == tsnk)
            {
                const PselProp &prop_snk = get_psel_prop_smear(job_tag, traj, xg_snk, 0, 0, 0);
                const WilsonMatrix &wm_y_snk = prop_snk.get_elem(xg_y_psel_idx);
                wm_snk_y_all[i_snk] =
                    gamma5 * (WilsonMatrix)matrix_adjoint(wm_y_snk) * gamma5;
                list_idx_from_isnk[i_snk] = xg_snk_psel_idx;
                i_snk += 1;
                continue;
            }
            const long xg_src_psel_idx = n;
            const Coordinate &xg_src = psel_smear[xg_src_psel_idx];
            if (xg_src[3] == tsrc)
            {
                const PselProp &prop_src = get_psel_prop_smear(job_tag, traj, xg_src, 0, 0, 0);
                const PselProp &prop_src_smear = get_psel_prop_smear(job_tag, traj, xg_src, 0, 0, 1);
                const WilsonMatrix &wm_y_src = prop_src.get_elem(xg_y_psel_idx);
                wm_y_src_all[i_src] = wm_y_src;
                list_idx_from_isrc[i_src] = xg_src_psel_idx;
                int i_snk1 = 0;
                for (long nsnk = 0; nsnk < n_points; ++nsnk)
                {
                    const long xg_snk_psel_idx1 = nsnk;
                    const Coordinate &xg_snk1 = psel_smear[xg_snk_psel_idx1];
                    if (xg_snk1[3] == tsnk)
                    {
                        const WilsonMatrix &wm_snk_src = prop_src_smear.get_elem(xg_snk_psel_idx1);
                        const int i_src_snk = i_src + num_src * i_snk1;
                        wm_snk_src_all[i_src_snk] = wm_snk_src;
                        i_snk1 += 1;
                    }
                }
                i_src += 1;
            }
        }
        return;
    }

    inline void contract_proton_sp_y_block_smear(std::vector<std::vector<ProtonSPBlock>> &spblock,
                                                 const std::vector<int> &psel_num_list_smear,
                                                 const std::string &job_tag, const int traj,
                                                 const Geometry &geo, const PointSelection &psel_smear,
                                                 const long &xg_y_psel_idx, const Coordinate &xg_y,
                                                 const int dtmin, const int tsep_src, const int tsep_snk)
    {
        TIMER_VERBOSE("contract_proton_sp_y_block_smear");
        const Coordinate &total_site = geo.total_site();
        const int num_dtxy = get_proton_fourp_max_sep(total_site[3]);
        qassert(num_dtxy == (int)spblock.size());
        const int yt = xg_y[3];
        for (int dt = 0; dt < num_dtxy; dt++)
        {
            const int dtxy = dt - num_dtxy / 2;
            int xt = mod(yt + dtxy, total_site[3]);
            int tsrc, tsnk;
            if (dtxy >= 0)
            {
                tsrc = mod(yt - tsep_src, total_site[3]);
                tsnk = mod(xt + tsep_snk, total_site[3]);
            }
            else
            {
                tsrc = mod(xt - tsep_src, total_site[3]);
                tsnk = mod(yt + tsep_snk, total_site[3]);
            }
            const int num_psel_src = psel_num_list_smear[tsrc];
            const int num_psel_snk = psel_num_list_smear[tsnk];
            const int num_psel = num_psel_src * num_psel_snk;
            spblock[dt].resize(num_psel);

            std::vector<WilsonMatrix> wm_snk_src_all;
            std::vector<WilsonMatrix> wm_y_src_all;
            std::vector<WilsonMatrix> wm_snk_y_all;
            std::vector<long> list_idx_from_isrc;
            std::vector<long> list_idx_from_isnk;
            get_all_wilson_matrix_smear(wm_snk_src_all, wm_y_src_all, wm_snk_y_all, job_tag, traj,
                                        xg_y_psel_idx, tsrc, tsnk, psel_smear, num_psel_src,
                                        num_psel_snk, list_idx_from_isrc, list_idx_from_isnk);

#pragma omp parallel for
            for (int i_src_snk = 0; i_src_snk < num_psel; ++i_src_snk)
            {
                const int i_src = mod(i_src_snk, num_psel_src);
                const int i_snk = (i_src_snk - i_src) / num_psel_src;
                const long &xg_src_psel_idx = list_idx_from_isrc[i_src];
                const long &xg_snk_psel_idx = list_idx_from_isnk[i_snk];
                const Coordinate &xg_src = psel_smear[xg_src_psel_idx];
                const Coordinate &xg_snk = psel_smear[xg_snk_psel_idx];
                const WilsonMatrix &wm_snk_src = wm_snk_src_all[i_src_snk];
                const WilsonMatrix &wm_y_src = wm_y_src_all[i_src];
                const WilsonMatrix &wm_snk_y = wm_snk_y_all[i_snk];
                ProtonSPBlock block1;
                proton_fourp_block_sp(block1.block_sp, wm_snk_src, wm_y_src, wm_snk_y);
                block1.xg_snk_psel_idx = xg_snk_psel_idx;
                block1.xg_src_psel_idx = xg_src_psel_idx;
                block1.n_snk = i_snk;
                block1.n_src = i_src;
                block1.is_build = true;
                spblock[dt][i_src_snk] = block1;
            }
        }
    }

    inline void contract_proton_fourp_unshifted_smear_acc_x_typeI(
        std::vector<Vector<Complex>> vv, const Complex coef, const WilsonMatrix &wm_x_snk,
        const WilsonMatrix &wm_x_y, const ProtonPPBlock &block)
    {
        const array<SpinMatrix, 8> &va_ms = get_va_matrices();
        const SpinMatrix &gamma5 = SpinMatrixConstants::get_gamma5();
        const WilsonMatrix wm_snk_x =
            gamma5 * (WilsonMatrix)matrix_adjoint(wm_x_snk) * gamma5;

        for (int gram = 0; gram < 5; gram++)
        {
            const WilsonMatrix wm_pp_y_tsrc_tsnk_x = block.block_pp[gram] * wm_snk_x;
            for (int mu = 0; mu < 8; mu++)
            {
                const WilsonMatrix wm_pp = wm_x_y * va_ms[mu] * wm_pp_y_tsrc_tsnk_x;
                for (int nu = 0; nu < 8; nu++)
                {
                    const int mu_nu = 8 * nu + mu;
                    vv[gram][mu_nu] += coef * matrix_trace(wm_pp, va_ms[nu]);
                }
            }
        }
        return;
    }

    inline void contract_proton_fourp_unshifted_smear_acc_x_typeII(
        std::vector<Vector<Complex>> vv, const Complex coef, const WilsonMatrix &wm_x_snk,
        const WilsonMatrix &wm_x_src, const ProtonSPBlock &block)
    {
        const SpinMatrix &gamma5 = SpinMatrixConstants::get_gamma5();
        const WilsonMatrix wm_snk_x =
            gamma5 * (WilsonMatrix)matrix_adjoint(wm_x_snk) * gamma5;
        const array<SpinMatrix, 8> &va_ms = get_va_matrices();

        for (int gram = 0; gram < 5; gram++)
        {
            for (int mu = 0; mu < 8; mu++)
            {
                const WilsonMatrix wm_sp = wm_x_src * block.block_sp[mu + gram * 8] * wm_snk_x;
                for (int nu = 0; nu < 8; nu++)
                {
                    const int mu_nu = 8 * nu + mu;
                    vv[gram][mu_nu] += coef * matrix_trace(wm_sp, va_ms[nu]);
                }
            }
        }
    }

    inline void contract_proton_fourp_unshifted_smear_typeI(std::vector<SelectedField<Complex>> &proton_fourp,
                                                            const std::string &job_tag, const int traj,
                                                            const SelProp &prop_x_y, const Coordinate &xg_y,
                                                            const long xg_y_psel_idx, const int tsep_src, const int tsep_snk,
                                                            const PointSelection &psel_smear, const FieldSelection &fsel,
                                                            const std::vector<std::vector<double>> psel_dt_num_list, const int dtmin,
                                                            const std::vector<std::vector<ProtonPPBlock>> &ppblock)
    {
        TIMER_VERBOSE("contract_proton_fourp_unshifted_smear_typeI");
        const Geometry &geo = fsel.f_rank.geo();
        const Coordinate &total_site = geo.total_site();
        const int num_dtxy = get_proton_fourp_max_sep(total_site[3]);
        const int num_dxxy = 33;
        const int yt = xg_y[3];
        const long n_points = psel_smear.size();
        const int multiplicity = 8 * 8;
        clear(proton_fourp);
        proton_fourp.resize(5);

        for (int gram = 0; gram < (int)proton_fourp.size(); gram++)
        {
            proton_fourp[gram].init(fsel, multiplicity);
            set_zero(proton_fourp[gram]);
        }
        for (int dt = 0; dt < num_dtxy; ++dt)
        {
            const int dtxy = dt - num_dtxy / 2;
            const int i_dt = abs(dtxy) + tsep_src + tsep_snk - dtmin;
            int xt = mod(yt + dtxy, total_site[3]);
            int tsrc, tsnk;
            if (dtxy >= 0)
            {
                tsrc = mod(yt - tsep_src, total_site[3]);
                tsnk = mod(xt + tsep_snk, total_site[3]);
            }
            else
            {
                tsrc = mod(xt - tsep_src, total_site[3]);
                tsnk = mod(yt + tsep_snk, total_site[3]);
            }
            const int dt_y = mod(yt - tsrc, total_site[3]);
            qassert(psel_dt_num_list[i_dt][dt_y] != 0.0);
            const Complex coef = 1.0 / psel_dt_num_list[i_dt][dt_y];
            // const Complex coef = 1.0;

            int idx_snk = 0;
            for (long nsnk = 0; nsnk < n_points; ++nsnk)
            {
                const long xg_snk_psel_idx = nsnk;
                const Coordinate &xg_snk = psel_smear[xg_snk_psel_idx];
                if (xg_snk[3] != tsnk)
                {
                    continue;
                }
                const SelProp &prop_x_snk = get_prop_smear(job_tag, traj, xg_snk, 0, 0, 0);
                const ProtonPPBlock &block1 = ppblock[dt][idx_snk];
                idx_snk += 1;
                if (block1.xg_snk_psel_idx == -1)
                {
                    continue;
                }
                qassert(block1.is_build);
                qassert(block1.xg_snk_psel_idx == xg_snk_psel_idx);
#pragma omp parallel for
                for (long idx = 0; idx < fsel.n_elems; ++idx)
                {
                    const long xg_x_idx = idx;
                    const long index = fsel.indices[xg_x_idx];
                    const Coordinate xl = geo.coordinate_from_index(index);
                    const Coordinate xg = geo.coordinate_g_from_l(xl);
                    const Coordinate &xg_x = xg;
                    if (xg_x[3] != xt)
                    {
                        continue;
                    }
                    const Coordinate xg_sep = smod(xg_x - xg_y, total_site);
                    if ((abs(xg_sep[0]) > num_dxxy / 2) || (abs(xg_sep[1]) > num_dxxy / 2) ||
                        (abs(xg_sep[2]) > num_dxxy / 2))
                    {
                        continue;
                    }
                    const WilsonMatrix &wm_x_y = prop_x_y.get_elem(idx);
                    const WilsonMatrix &wm_x_snk = prop_x_snk.get_elem(idx);
                    std::vector<Vector<Complex>> vv(5);
                    for (int gram = 0; gram < (int)proton_fourp.size(); gram++)
                    {
                        vv[gram] = proton_fourp[gram].get_elems(idx);
                    }
                    if (tsrc <= tsnk)
                    {
                        contract_proton_fourp_unshifted_smear_acc_x_typeI(vv, coef, wm_x_snk, wm_x_y, block1);
                    }
                    else
                    {
                        const Complex coef1 = -1.0 * coef;
                        contract_proton_fourp_unshifted_smear_acc_x_typeI(vv, coef1, wm_x_snk, wm_x_y, block1);
                    }
                }
            }
        }
        return;
    }

    inline void contract_proton_fourp_unshifted_smear_typeII(std::vector<SelectedField<Complex>> &proton_fourp,
                                                             const std::string &job_tag, const int traj,
                                                             const Coordinate &xg_y,
                                                             const long xg_y_psel_idx, const int tsep_src, const int tsep_snk,
                                                             const PointSelection &psel_smear, const FieldSelection &fsel,
                                                             const std::vector<int> psel_num_list_smear,
                                                             const std::vector<std::vector<double>> psel_dt_num_list,
                                                             const int dtmin,
                                                             const std::vector<std::vector<ProtonSPBlock>> &spblock)
    {
        TIMER_VERBOSE("contract_proton_fourp_unshifted_smear_typeII");
        const Geometry &geo = fsel.f_rank.geo();
        const Coordinate &total_site = geo.total_site();
        const int num_dtxy = get_proton_fourp_max_sep(total_site[3]);
        const int num_dxxy = 33;
        const int yt = xg_y[3];
        const long n_points = psel_smear.size();
        const int multiplicity = 8 * 8;
        clear(proton_fourp);
        proton_fourp.resize(5);

        const int x_sparsity = 4;

        for (int gram = 0; gram < (int)proton_fourp.size(); gram++)
        {
            proton_fourp[gram].init(fsel, multiplicity);
            set_zero(proton_fourp[gram]);
        }
        for (int dt = 0; dt < num_dtxy; ++dt)
        {
            if ((int)spblock[dt].size() == 0)
            {
                continue;
            }
            const int dtxy = dt - num_dtxy / 2;
            const int i_dt = abs(dtxy) + tsep_src + tsep_snk - dtmin;
            int xt = mod(yt + dtxy, total_site[3]);
            int tsrc, tsnk;
            if (dtxy >= 0)
            {
                tsrc = mod(yt - tsep_src, total_site[3]);
                tsnk = mod(xt + tsep_snk, total_site[3]);
            }
            else
            {
                tsrc = mod(xt - tsep_src, total_site[3]);
                tsnk = mod(yt + tsep_snk, total_site[3]);
            }
            const int dt_y = mod(yt - tsrc, total_site[3]);
            qassert(psel_dt_num_list[i_dt][dt_y] != 0.0);
            const Complex coef = 1.0 / psel_dt_num_list[i_dt][dt_y] * (double)x_sparsity;
            // const Complex coef = 1.0;
            const int num_src = psel_num_list_smear[tsrc];

            int idx_snk = 0;
            for (long nsnk = 0; nsnk < n_points; ++nsnk)
            {
                const long xg_snk_psel_idx = nsnk;
                const Coordinate &xg_snk = psel_smear[xg_snk_psel_idx];
                if (xg_snk[3] != tsnk)
                {
                    continue;
                }
                const SelProp &prop_x_snk = get_prop_smear(job_tag, traj, xg_snk, 0, 0, 0);
                int idx_src = 0;
                for (long nsrc = 0; nsrc < n_points; ++nsrc)
                {
                    const long xg_src_psel_idx = nsrc;
                    const Coordinate &xg_src = psel_smear[xg_src_psel_idx];
                    if (xg_src[3] != tsrc)
                    {
                        continue;
                    }
                    const SelProp &prop_x_src = get_prop_smear(job_tag, traj, xg_src, 0, 0, 0);
                    const int idx_src_snk = idx_src + idx_snk * num_src;
                    const ProtonSPBlock &block1 = spblock[dt][idx_src_snk];
                    qassert(block1.is_build);
                    qassert(block1.xg_snk_psel_idx == xg_snk_psel_idx);
                    qassert(block1.xg_src_psel_idx == xg_src_psel_idx);

#pragma omp parallel for
                    for (long idx = 0; idx < fsel.n_elems; ++idx)
                    {
                        if (mod(idx, x_sparsity) != 0)
                        {
                            continue;
                        }
                        const long xg_x_idx = idx;
                        const long index = fsel.indices[xg_x_idx];
                        const Coordinate xl = geo.coordinate_from_index(index);
                        const Coordinate xg = geo.coordinate_g_from_l(xl);
                        const Coordinate &xg_x = xg;
                        if (xg_x[3] != xt)
                        {
                            continue;
                        }
                        const Coordinate xg_sep = smod(xg_x - xg_y, total_site);
                        if ((abs(xg_sep[0]) > num_dxxy / 2) || (abs(xg_sep[1]) > num_dxxy / 2) ||
                            (abs(xg_sep[2]) > num_dxxy / 2))
                        {
                            continue;
                        }
                        const WilsonMatrix &wm_x_snk = prop_x_snk.get_elem(idx);
                        const WilsonMatrix &wm_x_src = prop_x_src.get_elem(idx);
                        std::vector<Vector<Complex>> vv(5);
                        for (int gram = 0; gram < (int)proton_fourp.size(); gram++)
                        {
                            vv[gram] = proton_fourp[gram].get_elems(idx);
                        }
                        if (tsrc <= tsnk)
                        {
                            contract_proton_fourp_unshifted_smear_acc_x_typeII(vv, coef, wm_x_snk, wm_x_src, block1);
                        }
                        else
                        {
                            const Complex coef1 = -1.0 * coef;
                            contract_proton_fourp_unshifted_smear_acc_x_typeII(vv, coef1, wm_x_snk, wm_x_src, block1);
                        }
                    }

                    idx_src += 1;
                }
                idx_snk += 1;
            }
        }
        return;
    }

    inline void
    contract_proton_four_point_acc_smear_typeI(
        std::vector<FieldM<Complex, 8 * 8>> &fourp_data,
        const std::string &job_tag, const int traj,
        const SelProp &prop_x_y, const Coordinate &xg_y,
        const long xg_y_psel_idx, const int tsep_src, const int tsep_snk, const PointSelection &psel_smear,
        const FieldSelection &fsel, const ShiftShufflePlan &ssp,
        const std::vector<int> &psel_num_list_smear,
        const std::vector<std::vector<double>> &psel_dt_num_list, const int dtmin,
        const std::vector<std::vector<ProtonPPBlock>> &ppblock)
    // xg_y = psel[xg_y_psel_idx] is the point src location for prop3_x_y
    // ssp = make_shift_shuffle_plan(fsel, -xg_y);
    {
        TIMER_VERBOSE("contract_proton_four_point_acc_psel");
        const Geometry &geo = fsel.f_rank.geo();
        qassert(fourp_data.size() == 10);
        qassert(geo == prop_x_y.geo());
        qassert(fsel.n_elems == prop_x_y.n_elems);
        qassert(is_initialized(prop_x_y));
        qassert(ssp.shift == -xg_y);
        qassert(ssp.is_reflect == false);
        std::vector<SelectedField<Complex>> sfs;
        contract_proton_fourp_unshifted_smear_typeI(sfs, job_tag, traj, prop_x_y, xg_y, xg_y_psel_idx,
                                                    tsep_src, tsep_snk, psel_smear, fsel, psel_dt_num_list, dtmin, ppblock);
        qassert(sfs.size() == 5);
        qassert(fsel.prob == ssp.fsel.prob);
        const Complex coef = 1.0 / fsel.prob;
        for (int gram = 0; gram < 5; gram++)
        {
            SelectedField<Complex> &s_data = sfs[gram];
            acc_field(fourp_data[gram], coef, s_data, ssp);
        }
    }

    inline void
    contract_proton_four_point_acc_smear_typeII(
        std::vector<FieldM<Complex, 8 * 8>> &fourp_data,
        const std::string &job_tag, const int traj,
        const SelProp &prop_x_y, const Coordinate &xg_y,
        const long xg_y_psel_idx, const int tsep_src, const int tsep_snk, const PointSelection &psel_smear,
        const FieldSelection &fsel, const ShiftShufflePlan &ssp,
        const std::vector<int> &psel_num_list_smear,
        const std::vector<std::vector<double>> &psel_dt_num_list, const int dtmin,
        const std::vector<std::vector<ProtonSPBlock>> &spblock)
    // xg_y = psel[xg_y_psel_idx] is the point src location for prop3_x_y
    // ssp = make_shift_shuffle_plan(fsel, -xg_y);
    {
        TIMER_VERBOSE("contract_proton_four_point_acc_psel");
        const Geometry &geo = fsel.f_rank.geo();
        qassert(fourp_data.size() == 10);
        qassert(geo == prop_x_y.geo());
        qassert(fsel.n_elems == prop_x_y.n_elems);
        qassert(is_initialized(prop_x_y));
        qassert(ssp.shift == -xg_y);
        qassert(ssp.is_reflect == false);
        std::vector<SelectedField<Complex>> sfs;
        contract_proton_fourp_unshifted_smear_typeII(sfs, job_tag, traj, xg_y, xg_y_psel_idx, tsep_src, tsep_snk, psel_smear,
                                                     fsel, psel_num_list_smear, psel_dt_num_list, dtmin, spblock);
        qassert(sfs.size() == 5);
        qassert(fsel.prob == ssp.fsel.prob);
        const Complex coef = 1.0 / fsel.prob;
        for (int gram = 0; gram < 5; gram++)
        {
            SelectedField<Complex> &s_data = sfs[gram];
            acc_field(fourp_data[gram + 5], coef, s_data, ssp);
        }
    }

    inline std::string get_proton_four_point_func_smear_path(const std::string &job_tag,
                                                             const int traj)
    {
        return ssprintf("analysis/proton_fourp_smear/%s/results=%d", job_tag.c_str(),
                        traj);
    }

    inline void compute_proton_four_point_func_smear_type(const std::string &job_tag, const int traj,
                                                          const std::vector<int> &gram_convert_list,
                                                          const std::vector<int> &tsep_list_src,
                                                          const std::vector<int> &tsep_list_snk)
    {
        check_sigterm();
        check_time_limit();
        Timer::autodisplay();
        const int gram_num = 10;
        const int tsep_num = (int)tsep_list_src.size();
        const int tsep_num_snk = (int)tsep_list_snk.size();
        qassert(tsep_num == tsep_num_snk);
        qassert(gram_convert_list.size() == gram_num);
        const std::string path = get_proton_four_point_func_smear_path(job_tag, traj);
        std::vector<std::string> fn_gram_list(gram_num * tsep_num);
        for (int gram = 0; gram < gram_num; ++gram)
        {
            const int &gram_name = gram_convert_list[gram];
            for (int ti = 0; ti < tsep_num; ++ti)
            {
                const int &tsep_src = tsep_list_src[ti];
                const int &tsep_snk = tsep_list_snk[ti];
                fn_gram_list[gram + gram_num * ti] = path + ssprintf("/fourp-data-dt%d-%d-%d.field", tsep_snk, tsep_src, gram_name);
            }
        }
        bool is_complete = true;
        for (int gram = 0; gram < gram_num; ++gram)
        {
            for (int ti = 0; ti < tsep_num; ++ti)
            {
                if (fn_gram_list[gram + gram_num * ti] != "")
                {
                    if (not is_d_field(fn_gram_list[gram + gram_num * ti]))
                    {
                        is_complete = false;
                        break;
                    }
                }
            }
        }
        if (is_complete)
        {
            return;
        }

        TIMER_VERBOSE("compute_proton_four_point_func_smear_type");
        const PointSelection &psel = get_point_selection(job_tag, traj);
        const long n_points = psel.size();
        const PointSelection &psel_smear = get_point_selection_smear(job_tag, traj);
        const long n_points_smear = psel_smear.size();
        const FieldSelection &fsel = get_field_selection(job_tag, traj);
        std::vector<std::vector<FieldM<Complex, 8 * 8>>> fourp_data(tsep_num);
        const Geometry &geo = fsel.f_rank.geo();
        const Coordinate total_site = geo.total_site();
        const int num_dtxy = get_proton_fourp_max_sep(total_site[3]);
        const int dtmax = num_dtxy / 2 + 2 * tsep_list_src[tsep_num - 1];
        const int dtmin = 2 * tsep_list_src[0];
        for (int t = 0; t < tsep_num; ++t)
        {
            fourp_data[t].resize(gram_num);
        }

        std::vector<int> psel_num_list_smear;
        std::vector<std::vector<double>> psel_dt_num_list;
        std::vector<std::vector<std::vector<ProtonPPBlock>>> block;
        contract_proton_pp_block_smear(block, psel_num_list_smear, psel_dt_num_list, job_tag, traj, geo, psel, psel_smear, dtmax, dtmin);

        const int y_idx_sparsity = 1;
        long iter = 0;
        for (int yt = 0; yt < total_site[3]; ++yt)
        {
            load_prop_smear_fourp(job_tag, traj, tsep_list_src.back(), num_dtxy, yt, psel_smear, geo);
            for (long ny = 0; ny < n_points; ++ny)
            {
                const long xg_y_psel_idx = ny;
                const Coordinate &xg_y = psel[xg_y_psel_idx];
                if (yt != xg_y[3])
                {
                    continue;
                }
                Timer::autodisplay();
                TIMER_VERBOSE("compute_proton_four_point_func_noama_type-iter");
                iter += 1;
                const SelProp &prop_x_y = get_prop_psrc(job_tag, traj, xg_y, 0, 0);
                const ShiftShufflePlan ssp = make_shift_shuffle_plan(fsel, -xg_y);
                displayln_info(fname + ssprintf(":n=%ld iter=%ld", ny,
                                                iter));

                for (int ti = 0; ti < tsep_num; ++ti)
                {
                    const int tsep_src = tsep_list_src[ti];
                    const int tsep_snk = tsep_list_snk[ti];
                    std::vector<std::vector<ProtonPPBlock>> ppblock(num_dtxy);
                    std::vector<std::vector<ProtonSPBlock>> spblock(num_dtxy);
                    contract_proton_pp_y_block_smear(block, ppblock, psel_num_list_smear, job_tag, traj, geo, psel_smear,
                                                     xg_y_psel_idx, xg_y, dtmin, tsep_src, tsep_snk);
                    contract_proton_four_point_acc_smear_typeI(fourp_data[ti], job_tag, traj, prop_x_y, xg_y,
                                                               xg_y_psel_idx, tsep_src, tsep_snk, psel_smear, fsel, ssp, psel_num_list_smear,
                                                               psel_dt_num_list, dtmin, ppblock);
                    const int type = mod((int)xg_y_psel_idx, y_idx_sparsity);
                    if (type == 0)
                    {
                        contract_proton_sp_y_block_smear(spblock, psel_num_list_smear, job_tag, traj, geo, psel_smear,
                                                         xg_y_psel_idx, xg_y, dtmin, tsep_src, tsep_snk);
                        contract_proton_four_point_acc_smear_typeII(fourp_data[ti], job_tag, traj, prop_x_y, xg_y,
                                                                    xg_y_psel_idx, tsep_src, tsep_snk, psel_smear, fsel, ssp, psel_num_list_smear,
                                                                    psel_dt_num_list, dtmin, spblock);
                    }
                }
            }
        }

        for (int ti = 0; ti < tsep_num; ++ti)
        {
            for (int gram = 0; gram < 5; gram++)
            {
                const double coef_wall = (double)total_site[0] * (double)total_site[1] * (double)total_site[2];
                const double coef = coef_wall * coef_wall;
                fourp_data[ti][gram] *= coef;
                const double coef1 = coef * (double)y_idx_sparsity;
                fourp_data[ti][gram + 5] *= coef1;
            }
        }
        FieldM<Complex, 8 * 8> avg;
        for (int gram = 0; gram < gram_num; gram++)
        {
            for (int ti = 0; ti < tsep_num; ++ti)
            {
                qassert(is_initialized(fourp_data[ti][gram]));
                qassert(fn_gram_list[gram + gram_num * ti] != "");
                avg = fourp_data[ti][gram];
                write_field_float_from_double(avg, fn_gram_list[gram + gram_num * ti]);
            }
        }
    }

    inline void compute_proton_four_point_func_smear(const std::string &job_tag,
                                                     const int traj)
    {
        check_sigterm();
        Timer::autodisplay();
        TIMER_VERBOSE("compute_proton_four_point_func_smear");
        const std::string path = get_proton_four_point_func_smear_path(job_tag, traj);
        const std::string path_checkpoint = path + "/checkpoint.txt";
        if (does_file_exist_sync_node(path_checkpoint))
        {
            return;
        }
        if (does_file_exist_sync_node(path + "/fourp-data-dt3-2-7.field"))
        {
            return;
        }
        if (not(check_prop_smear(job_tag, traj, 0)))
        {
            return;
        }
        check_sigterm();
        check_time_limit();
        if (not obtain_lock(
                ssprintf("lock-proton-fourp-smear-%s-dt6-%d", job_tag.c_str(), traj)))
        {
            return;
        }
        setup(job_tag, traj);
        qmkdir_info("analysis/proton_fourp_smear");
        qmkdir_info(ssprintf("analysis/proton_fourp_smear/%s", job_tag.c_str()));
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
        std::vector<int> tsep_list_src, tsep_list_snk;
        //tsep_list_src.push_back(3); // small to large
        tsep_list_src.push_back(2); // small to large
        tsep_list_src.push_back(3); // small to large
        //tsep_list_src.push_back(4); // small to large
        //tsep_list_src.push_back(4); // small to large
	//tsep_list_snk.push_back(4);
	//tsep_list_snk.push_back(4);
	tsep_list_snk.push_back(3);
	tsep_list_snk.push_back(2);
	//tsep_list_snk.push_back(3);
        compute_proton_four_point_func_smear_type(job_tag, traj, gram_convert_list, tsep_list_src, tsep_list_snk);
        release_lock();
    }
}
