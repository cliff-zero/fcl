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
#endif

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

    inline void load_prop_psrc_fourp(const std::string &job_tag, const int traj, const int tsep,
                                     const int num_dtxy, const int yt, const PointSelection &psel,
                                     const Geometry &geo)
    {
        TIMER_VERBOSE("load_prop_psrc_fourp");
        const long n_points = psel.size();
        const Coordinate total_site = geo.total_site();
        const int num_t = 2 * tsep + num_dtxy;
        //const int num_t = total_site[3];
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
                const SelProp &prop = get_prop_psrc(job_tag, traj, xg_src, 0, 0);
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
            const SelProp &prop = get_prop_psrc(job_tag, traj, xg_src, 0, 0);
            qassert(prop.initialized);
            break;
        }
    }

    inline void contract_proton_pp_block(std::vector<std::vector<std::vector<ProtonPPBlock>>> &block,
                                         std::vector<int> &psel_num_list,
                                         std::vector<std::vector<long>> &psel_dt_num_list,
                                         const std::string &job_tag, const int traj,
                                         const Geometry &geo, const PointSelection &psel,
                                         const int dtmax, const int dtmin)
    {
        TIMER_VERBOSE("contract_proton_pp_block");
        const int num_dt = dtmax - dtmin + 1;
        qassert(num_dt > 0);
        const long n_points = psel.size();
        const Coordinate total_site = geo.total_site();
        block.resize(n_points);
        for (long n = 0; n < n_points; ++n)
        {
            block[n].resize(num_dt);
        }

        psel_num_list.resize(total_site[3]);
        set_zero(psel_num_list);
        psel_dt_num_list.resize(num_dt);
        for (int dt = 0; dt < num_dt; ++dt)
        {
            const int del_t = dtmin + dt;
            psel_dt_num_list[dt].resize(del_t + 1);
            set_zero(psel_dt_num_list[dt]);
        }
        std::vector<int> list_n_from_idx(n_points);
        for (long n = 0; n < n_points; ++n)
        {
            const long xg_psel_idx = n;
            const Coordinate xg = psel[xg_psel_idx];
            const int t = xg[3];
            list_n_from_idx[xg_psel_idx] = psel_num_list[t];
            psel_num_list[t] += 1;
        }
        for (long nsnk = 0; nsnk < n_points; ++nsnk)
        {
            const long xg_snk_psel_idx = nsnk;
            const Coordinate &xg_snk = psel[xg_snk_psel_idx];
            const int tsnk = xg_snk[3];
            for (int dt = 0; dt < num_dt; ++dt)
            {
                const int del_t = dtmin + dt;
                const int tsrc = mod(tsnk - del_t, total_site[3]);
                block[xg_snk_psel_idx][dt].resize(psel_num_list[tsrc]);
            }
        }
        for (int tsrc = 0; tsrc < total_site[3]; ++tsrc)
        {
            for (int dt = 0; dt < num_dt; ++dt)
            {
                const int del_t = dtmin + dt;
                const int tsnk = mod(tsrc + del_t, total_site[3]);
                const int num_src_snk = psel_num_list[tsrc] * psel_num_list[tsnk];
                for (int dt_y = 0; dt_y <= del_t; ++dt_y)
                {
                    const int yt = mod(tsrc + dt_y, total_site[3]);
                    const int num_src_y_snk = num_src_snk * psel_num_list[yt];
                    psel_dt_num_list[dt][dt_y] += num_src_y_snk;
                }
            }
        }

        for (long nsrc = 0; nsrc < n_points; ++nsrc)
        {
            const long xg_src_psel_idx = nsrc;
            const Coordinate &xg_src = psel[xg_src_psel_idx];
            const int tsrc = xg_src[3];
            const PselProp &prop_src = get_psel_prop_psrc(job_tag, traj, xg_src, 0, 0);

#pragma omp parallel for
            for (long nsnk = 0; nsnk < n_points; ++nsnk)
            {
                const long xg_snk_psel_idx = nsnk;
                const Coordinate &xg_snk = psel[xg_snk_psel_idx];
                const int tsnk = xg_snk[3];
                const int del_t = mod(tsnk - tsrc, total_site[3]);
                const int dt = del_t - dtmin;
                if ((dt < 0) || (dt >= num_dt))
                {
                    continue;
                }
                ProtonPPBlock ppblock;
                const WilsonMatrix &wm = prop_src.get_elem(xg_snk_psel_idx);
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

    inline void contract_proton_fourp_unshifted_psel_acc_x_typeI(
        std::vector<Vector<Complex>> vv, const Complex coef, const WilsonMatrix &wm_x_src,
        const WilsonMatrix &wm_x_snk, const ProtonPPBlock &block)
    {
        const array<SpinMatrix, 8> &va_ms = get_va_matrices();
        const SpinMatrix &gamma5 = SpinMatrixConstants::get_gamma5();
        const WilsonMatrix wm_snk_x =
            gamma5 * (WilsonMatrix)matrix_adjoint(wm_x_snk) * gamma5;

        for (int gram = 0; gram < 5; gram++)
        {
            const WilsonMatrix wm_pp_tsrc_tsnk_x = block.block_pp[gram] * wm_snk_x;
            const WilsonMatrix wm_pp = wm_x_src * wm_pp_tsrc_tsnk_x;
            for (int mu = 0; mu < 8; mu++)
            {
                vv[gram][mu] += coef * matrix_trace(wm_pp, va_ms[mu]);
            }
        }
        return;
    }

    inline void contract_proton_fourp_unshifted_psel_typeI(std::vector<SelectedField<Complex>> &proton_fourp,
                                                           const std::string &job_tag, const int traj,
                                                           const SelProp &prop_x_y, const Coordinate &xg_y,
                                                           const long xg_y_psel_idx, const int tsep,
                                                           const PointSelection &psel, const FieldSelection &fsel,
                                                           const std::vector<int> &psel_num_list,
                                                           const std::vector<std::vector<long>> psel_dt_num_list, const int dtmin,
                                                           const std::vector<std::vector<std::vector<ProtonPPBlock>>> &block)
    {
        TIMER_VERBOSE("contract_proton_fourp_unshifted_psel_typeI");
        const Geometry &geo = fsel.f_rank.geo();
        const Coordinate &total_site = geo.total_site();
        const int num_dtxy = 13;
        const int num_dxxy = total_site[0];
        const int yt = xg_y[3];
        const long n_points = psel.size();
        const int multiplicity = 8;
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
            const int i_dt = abs(dtxy) + 2 * tsep - dtmin;
            int xt = mod(yt + dtxy, total_site[3]);
            int tsrc, tsnk;
            if (dtxy >= 0)
            {
                tsrc = mod(yt - tsep, total_site[3]);
                tsnk = mod(xt + tsep, total_site[3]);
            }
            else
            {
                tsrc = mod(xt - tsep, total_site[3]);
                tsnk = mod(yt + tsep, total_site[3]);
            }
            const int dt_y = mod(yt - tsrc, total_site[3]);
            //const Complex coef = 1.0 / (double)psel_dt_num_list[i_dt][dt_y];
            const Complex coef = 1.0;

            int idx_snk = 0;
            for (long nsnk = 0; nsnk < n_points; ++nsnk)
            {
                const long xg_snk_psel_idx = nsnk;
                const Coordinate &xg_snk = psel[xg_snk_psel_idx];
                if (xg_snk[3] != tsnk)
                {
                    continue;
                }
                const SelProp &prop_x_snk = get_prop_psrc(job_tag, traj, xg_snk, 0, 0);
                idx_snk += 1;
                const int num_psel_src = psel_num_list[tsrc];
                for (int nt = 0; nt < num_psel_src; ++nt)
                {
                    const int i_dtss = mod(tsnk - tsrc, total_site[3]) - dtmin;
                    const ProtonPPBlock &block1 = block[xg_snk_psel_idx][i_dtss][nt];
                    qassert(block1.is_build);
                    qassert(block1.xg_snk_psel_idx == xg_snk_psel_idx);
                    const int i_snk = block1.n_snk;
                    const long xg_src_psel_idx = block1.xg_src_psel_idx;
                    const Coordinate &xg_src = psel[xg_src_psel_idx];
                    if (block1.xg_snk_psel_idx == -1)
                    {
                        continue;
                    }
                    const SelProp &prop_x_src = get_prop_psrc(job_tag, traj, xg_src, 0, 0);
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
                        const WilsonMatrix &wm_x_src = prop_x_src.get_elem(idx);
                        const WilsonMatrix &wm_x_snk = prop_x_snk.get_elem(idx);
                        std::vector<Vector<Complex>> vv(5);
                        for (int gram = 0; gram < (int)proton_fourp.size(); gram++)
                        {
                            vv[gram] = proton_fourp[gram].get_elems(idx);
                        }
                        if (tsrc <= tsnk)
                        {
                            contract_proton_fourp_unshifted_psel_acc_x_typeI(vv, coef, wm_x_src, wm_x_snk, block1);
                        }
                        else
                        {
                            const Complex coef1 = -1.0 * coef;
                            contract_proton_fourp_unshifted_psel_acc_x_typeI(vv, coef1, wm_x_src, wm_x_snk, block1);
                        }
                    }
                }
            }
        }
        return;
    }

    inline void
    contract_proton_four_point_acc_psel_typeI(
        std::vector<FieldM<Complex, 8>> &fourp_data,
        const std::string &job_tag, const int traj,
        const SelProp &prop_x_y, const Coordinate &xg_y,
        const long xg_y_psel_idx, const int tsep, const PointSelection &psel,
        const FieldSelection &fsel, const ShiftShufflePlan &ssp,
        const std::vector<int> &psel_num_list,
        const std::vector<std::vector<long>> &psel_dt_num_list, const int dtmin,
        const std::vector<std::vector<std::vector<ProtonPPBlock>>> &block)
    // xg_y = psel[xg_y_psel_idx] is the point src location for prop3_x_y
    // ssp = make_shift_shuffle_plan(fsel, -xg_y);
    {
        TIMER_VERBOSE("contract_proton_four_point_acc_psel");
        const Geometry &geo = fsel.f_rank.geo();
        qassert(fourp_data.size() == 5);
        qassert(geo == prop_x_y.geo());
        qassert(fsel.n_elems == prop_x_y.n_elems);
        qassert(is_initialized(prop_x_y));
        qassert(psel[xg_y_psel_idx] == xg_y);
        qassert(ssp.shift == -xg_y);
        qassert(ssp.is_reflect == false);
        std::vector<SelectedField<Complex>> sfs;
        contract_proton_fourp_unshifted_psel_typeI(sfs, job_tag, traj, prop_x_y, xg_y, xg_y_psel_idx,
                                                   tsep, psel, fsel, psel_num_list, psel_dt_num_list, dtmin, block);
        qassert(sfs.size() == 5);
        qassert(fsel.prob == ssp.fsel.prob);
        const Complex coef = 1.0 / fsel.prob;
        for (int gram = 0; gram < 5; gram++)
        {
            SelectedField<Complex> &s_data = sfs[gram];
            acc_field(fourp_data[gram], coef, s_data, ssp);
        }
    }

    inline std::string get_proton_four_point_func_psel_path(const std::string &job_tag,
                                                            const int traj)
    {
        return ssprintf("analysis/proton_threep_from_fourp/%s/results=%d", job_tag.c_str(),
                        traj);
    }

    inline void compute_proton_four_point_func_psel_type(const std::string &job_tag, const int traj,
                                                         const std::vector<int> &gram_convert_list,
                                                         const std::vector<int> &tsep_list)
    {
        check_sigterm();
        check_time_limit();
        Timer::autodisplay();
        const int gram_num = 5;
        const int tsep_num = (int)tsep_list.size();
        qassert(gram_convert_list.size() == gram_num);
        const std::string path = get_proton_four_point_func_psel_path(job_tag, traj);
        std::vector<std::string> fn_gram_list(gram_num * tsep_num);
        for (int gram = 0; gram < gram_num; ++gram)
        {
            const int &gram_name = gram_convert_list[gram];
            for (int ti = 0; ti < tsep_num; ++ti)
            {
                const int &tsep = tsep_list[ti];
                fn_gram_list[gram + gram_num * ti] = path + ssprintf("/fourp-data-dt%d-%d.field", tsep, gram_name);
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
        TIMER_VERBOSE("compute_proton_four_point_func_psel_type");
        const PointSelection &psel = get_point_selection(job_tag, traj);
        const long n_points = psel.size();
        const FieldSelection &fsel = get_field_selection(job_tag, traj);
        std::vector<std::vector<FieldM<Complex, 8>>> fourp_data(tsep_num);
        const Geometry &geo = fsel.f_rank.geo();
        const Coordinate total_site = geo.total_site();
        const int num_dtxy = 13;
        const int dtmax = num_dtxy / 2 + 2 * tsep_list[tsep_num - 1];
        const int dtmin = 2 * tsep_list[0];
        for (int t = 0; t < tsep_num; ++t)
        {
            fourp_data[t].resize(gram_num);
        }

        std::vector<int> psel_num_list;
        std::vector<std::vector<long>> psel_dt_num_list;
        std::vector<std::vector<std::vector<ProtonPPBlock>>> block;
        contract_proton_pp_block(block, psel_num_list, psel_dt_num_list, job_tag, traj, geo, psel, dtmax, dtmin);

        long iter = 0;
        for (int yt = 0; yt < total_site[3]; ++yt)
        {
            load_prop_psrc_fourp(job_tag, traj, tsep_list.back(), num_dtxy, yt, psel, geo);
            for (long ny = 0; ny < n_points; ++ny)
            {
                const long xg_y_psel_idx = ny;
                const Coordinate &xg_y = psel[xg_y_psel_idx];
                if (yt != xg_y[3])
                {
                    continue;
                }
                if (get_point_src_info(job_tag, traj, xg_y, 0).size() == 0)
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
                    const int tsep = tsep_list[ti];
                    contract_proton_four_point_acc_psel_typeI(fourp_data[ti], job_tag, traj, prop_x_y, xg_y,
                                                              xg_y_psel_idx, tsep, psel, fsel, ssp, psel_num_list,
                                                              psel_dt_num_list, dtmin, block);
                }
            }
        }

        for (int ti = 0; ti < tsep_num; ++ti)
        {
            for (int gram = 0; gram < 5; gram++)
            {
                const double coef_wall = (double)total_site[0] * (double)total_site[1] * (double)total_site[2];
                const double coef = 1.0; // coef_wall * coef_wall;
                fourp_data[ti][gram] *= coef;
            }
        }
        FieldM<Complex, 8> avg;
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

    inline void compute_proton_four_point_func_psel(const std::string &job_tag,
                                                    const int traj)
    {
        check_sigterm();
        Timer::autodisplay();
        TIMER_VERBOSE("compute_proton_four_point_func_psel");
        const std::string path = get_proton_four_point_func_psel_path(job_tag, traj);
        const std::string path_checkpoint = path + "/checkpoint.txt";
        if (does_file_exist_sync_node(path_checkpoint))
        {
            return;
        }
        if (does_file_exist_sync_node(path + "/fourp-data-dt2-5.field"))
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
                ssprintf("lock-proton-threep-from-fourp-%s-%d", job_tag.c_str(), traj)))
        {
            return;
        }
        setup(job_tag, traj);
        qmkdir_info("analysis/proton_threep_from_fourp");
        qmkdir_info(ssprintf("analysis/proton_threep_from_fourp/%s", job_tag.c_str()));
        qmkdir_info(path);
        std::vector<int> gram_convert_list;
        gram_convert_list.push_back(1);
        gram_convert_list.push_back(2);
        gram_convert_list.push_back(3);
        gram_convert_list.push_back(4);
        gram_convert_list.push_back(5);
        std::vector<int> tsep_list;
        tsep_list.push_back(2);
        compute_proton_four_point_func_psel_type(job_tag, traj, gram_convert_list, tsep_list);
        release_lock();
    }
}
