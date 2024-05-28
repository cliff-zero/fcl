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
    inline void load_prop_smear_threep(const std::string &job_tag, const int traj,
                                       const int dtmax, const int tsrc, const PointSelection &psel,
                                       const Geometry &geo)
    {
        TIMER_VERBOSE("load_prop_smear_threep");
        const long n_points = psel.size();
        const Coordinate total_site = geo.total_site();
        const int num_t = dtmax + 1;

        for (int ti = 0; ti < num_t; ++ti)
        {
            const int tslice = mod(tsrc + ti, total_site[3]);
            for (long nsrc = 0; nsrc < n_points; ++nsrc)
            {
                const long xg_src_psel_idx = nsrc;
                const Coordinate &xg_src = psel[xg_src_psel_idx];
                if (xg_src[3] != tslice)
                {
                    continue;
                }
                const SelProp &prop = get_prop_smear(job_tag, traj, xg_src, 0, 0, 0);
                qassert(prop.initialized);
            }
        }
    }

    inline void contract_proton_pp_block(std::vector<std::vector<std::vector<ProtonPPBlock>>> &block,
                                         std::vector<int> &psel_num_list,
                                         std::vector<double> &psel_dt_num_list,
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
        set_zero(psel_dt_num_list);

        std::vector<int> list_n_from_idx(n_points);
        for (long n = 0; n < n_points; ++n)
        {
            const long xg_psel_idx = n;
            const Coordinate xg = psel[xg_psel_idx];
            const int t = xg[3];
            list_n_from_idx[xg_psel_idx] = psel_num_list[t];
            psel_num_list[t] += 1;
        }
        for (long nsrc = 0; nsrc < n_points; ++nsrc)
        {
            const long xg_src_psel_idx = nsrc;
            const Coordinate &xg_src = psel[xg_src_psel_idx];
            const int tsrc = xg_src[3];
            for (int dt = 0; dt < num_dt; ++dt)
            {
                const int del_t = dtmin + dt;
                const int tsnk = mod(tsrc + del_t, total_site[3]);
                block[xg_src_psel_idx][dt].resize(psel_num_list[tsnk]);
            }
        }
        for (long nsrc = 0; nsrc < n_points; ++nsrc)
        {
            const long xg_src_psel_idx = nsrc;
            const Coordinate &xg_src = psel[xg_src_psel_idx];
            const int tsrc = xg_src[3];
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
                psel_dt_num_list[dt] += 1.0;
            }
        }

        for (long nsrc = 0; nsrc < n_points; ++nsrc)
        {
            const long xg_src_psel_idx = nsrc;
            const Coordinate &xg_src = psel[xg_src_psel_idx];
            const int tsrc = xg_src[3];
            const PselProp &prop_src = get_psel_prop_smear(job_tag, traj, xg_src, 0, 0, 1);

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

                const int n_t = list_n_from_idx[xg_snk_psel_idx];
                block[xg_src_psel_idx][dt][n_t] = ppblock;
            }
        }
        return;
    }

    inline void contract_proton_threep_unshifted_smear_acc_x(
        std::vector<Vector<Complex>> vv, const Complex coef, const WilsonMatrix &wm_x_snk,
        const WilsonMatrix &wm_x_src, const ProtonPPBlock &block)
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

    inline void contract_proton_threep_unshifted_smear(std::vector<SelectedField<Complex>> &proton_threep,
                                                       const std::string &job_tag, const int traj,
                                                       const SelProp &prop_x_src, const Coordinate &xg_src,
                                                       const long xg_src_psel_idx, const PointSelection &psel,
                                                       const FieldSelection &fsel,
                                                       double &psel_dt_num_list, const int tsnk,
                                                       const std::vector<ProtonPPBlock> &ppblock)
    {
        TIMER_VERBOSE("contract_proton_threep_unshifted_smear");
        const Geometry &geo = fsel.f_rank.geo();
        const Coordinate &total_site = geo.total_site();
        const int tsrc = xg_src[3];
        const int del_t = mod(tsnk - tsrc, total_site[3]);
        const long n_points = psel.size();
        const int multiplicity = 8;
        clear(proton_threep);
        proton_threep.resize(5);

        for (int gram = 0; gram < (int)proton_threep.size(); gram++)
        {
            proton_threep[gram].init(fsel, multiplicity);
            set_zero(proton_threep[gram]);
        }

        const Complex coef = 1.0 / psel_dt_num_list;

        int idx_snk = 0;
        for (long nsnk = 0; nsnk < n_points; ++nsnk)
        {
            const long xg_snk_psel_idx = nsnk;
            const Coordinate &xg_snk = psel[xg_snk_psel_idx];
            if (xg_snk[3] != tsnk)
            {
                continue;
            }
            const SelProp &prop_x_snk = get_prop_smear(job_tag, traj, xg_snk, 0, 0, 0);
            const ProtonPPBlock &block1 = ppblock[idx_snk];
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
                const int dt_src_x = mod(xg_x[3] - tsrc, total_site[3]);
                if ((dt_src_x <= 0) || (dt_src_x >= del_t))
                {
                    continue;
                }
                const WilsonMatrix &wm_x_src = prop_x_src.get_elem(idx);
                const WilsonMatrix &wm_x_snk = prop_x_snk.get_elem(idx);
                std::vector<Vector<Complex>> vv(5);
                for (int gram = 0; gram < (int)proton_threep.size(); gram++)
                {
                    vv[gram] = proton_threep[gram].get_elems(idx);
                }
                if (tsrc <= tsnk)
                {
                    contract_proton_threep_unshifted_smear_acc_x(vv, coef, wm_x_snk, wm_x_src, block1);
                }
                else
                {
                    const Complex coef1 = -1.0 * coef;
                    contract_proton_threep_unshifted_smear_acc_x(vv, coef1, wm_x_snk, wm_x_src, block1);
                }
            }
        }

        return;
    }

    inline void
    contract_proton_three_point_acc_smear(
        std::vector<FieldM<Complex, 8>> &threep_data,
        const std::string &job_tag, const int traj,
        const SelProp &prop_x_src, const Coordinate &xg_src,
        const long xg_src_psel_idx, const PointSelection &psel,
        const FieldSelection &fsel, const ShiftShufflePlan &ssp,
        const std::vector<int> &psel_num_list,
        double &psel_dt_num_list, const int tsnk,
        const std::vector<ProtonPPBlock> &ppblock)
    // xg_y = psel[xg_y_psel_idx] is the point src location for prop3_x_y
    // ssp = make_shift_shuffle_plan(fsel, -xg_y);
    {
        TIMER_VERBOSE("contract_proton_three_point_acc_smear");
        const Geometry &geo = fsel.f_rank.geo();
        qassert(threep_data.size() == 5);
        qassert(geo == prop_x_src.geo());
        qassert(fsel.n_elems == prop_x_src.n_elems);
        qassert(is_initialized(prop_x_src));
        qassert(psel[xg_src_psel_idx] == xg_src);
        qassert(ssp.shift == -xg_src);
        qassert(ssp.is_reflect == false);
        std::vector<SelectedField<Complex>> sfs;
        contract_proton_threep_unshifted_smear(sfs, job_tag, traj, prop_x_src, xg_src, xg_src_psel_idx,
                                               psel, fsel, psel_dt_num_list, tsnk, ppblock);
        qassert(sfs.size() == 5);
        qassert(fsel.prob == ssp.fsel.prob);
        const Complex coef = 1.0 / fsel.prob;
        for (int gram = 0; gram < 5; gram++)
        {
            SelectedField<Complex> &s_data = sfs[gram];
            acc_field(threep_data[gram], coef, s_data, ssp);
        }
    }

    inline std::string get_proton_three_point_func_smear_path(const std::string &job_tag,
                                                              const int traj)
    {
        return ssprintf("analysis/proton_threep_smear/%s/results=%d", job_tag.c_str(),
                        traj);
    }

    inline void compute_proton_three_point_func_smear_type(const std::string &job_tag, const int traj,
                                                           const std::vector<int> &gram_convert_list)
    {
        check_sigterm();
        check_time_limit();
        Timer::autodisplay();

        const int dtmax = 10;
        const int dtmin = 2;
        qassert(dtmin <= dtmax);
        const int num_t = dtmax - dtmin + 1;

        const int gram_num = 5;
        qassert(gram_convert_list.size() == gram_num);
        const std::string path = get_proton_three_point_func_smear_path(job_tag, traj);
        std::vector<std::vector<std::string>> fn_gram_list(num_t);
        for (int dt = 0; dt < num_t; ++dt)
        {
            const int tsep = dt + dtmin;
            fn_gram_list[dt].resize(5);
            for (int gram = 0; gram < gram_num; ++gram)
            {
                const int &gram_name = gram_convert_list[gram];

                fn_gram_list[dt][gram] = path + ssprintf("/threep-data-dt%d-%d.field", tsep, gram_name);
            }
        }
        bool is_complete = true;
        for (int dt = 0; dt < num_t; ++dt)
        {
            for (int gram = 0; gram < gram_num; ++gram)
            {
                if (fn_gram_list[dt][gram] != "")
                {
                    if (not is_d_field(fn_gram_list[dt][gram]))
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

        TIMER_VERBOSE("compute_proton_three_point_func_smear_type");
        const PointSelection &psel = get_point_selection_smear(job_tag, traj);
        const long n_points = psel.size();
        const FieldSelection &fsel = get_field_selection(job_tag, traj);
        std::vector<std::vector<FieldM<Complex, 8>>> threep_data(num_t);
        for (int dt = 0; dt < num_t; ++dt)
        {
            threep_data[dt].resize(gram_num);
        }
        const Geometry &geo = fsel.f_rank.geo();
        const Coordinate total_site = geo.total_site();

        std::vector<int> psel_num_list;
        std::vector<double> psel_dt_num_list;
        std::vector<std::vector<std::vector<ProtonPPBlock>>> block;
        contract_proton_pp_block(block, psel_num_list, psel_dt_num_list, job_tag, traj, geo, psel, dtmax, dtmin);
        displayln_info(fname + ssprintf("test_num: %.1f", psel_dt_num_list[0]));
        long iter = 0;
        for (int tsrc = 0; tsrc < total_site[3]; ++tsrc)
        {
            load_prop_smear_threep(job_tag, traj, dtmax, tsrc, psel, geo);
            for (long nsrc = 0; nsrc < n_points; ++nsrc)
            {
                const long xg_src_psel_idx = nsrc;
                const Coordinate &xg_src = psel[xg_src_psel_idx];
                if (tsrc != xg_src[3])
                {
                    continue;
                }
                Timer::autodisplay();
                TIMER_VERBOSE("compute_proton_three_point_func_noama_type-iter");
                iter += 1;
                const SelProp &prop_x_src = get_prop_smear(job_tag, traj, xg_src, 0, 0, 0);
                const ShiftShufflePlan ssp = make_shift_shuffle_plan(fsel, -xg_src);
                displayln_info(fname + ssprintf(":n=%ld iter=%ld", nsrc,
                                                iter));

                for (int tsep = dtmin; tsep <= dtmax; ++tsep)
                {
                    const int tsnk = mod(tsrc + tsep, total_site[3]);
                    const int dt = tsep - dtmin;
                    const std::vector<ProtonPPBlock> &block_isrc_dt = block[xg_src_psel_idx][dt];
                    contract_proton_three_point_acc_smear(threep_data[dt], job_tag, traj, prop_x_src, xg_src,
                                                          xg_src_psel_idx, psel, fsel, ssp, psel_num_list,
                                                          psel_dt_num_list[dt], tsnk, block_isrc_dt);
                }
            }
        }

        for (int dt = 0; dt < num_t; ++dt)
        {
            for (int gram = 0; gram < gram_num; ++gram)
            {
                const double coef_wall = (double)total_site[0] * (double)total_site[1] * (double)total_site[2];
                const double coef = coef_wall * coef_wall;
                threep_data[dt][gram] *= coef;
            }
        }
        FieldM<Complex, 8> avg;
        for (int gram = 0; gram < gram_num; gram++)
        {
            for (int dt = 0; dt < num_t; ++dt)
            {
                qassert(is_initialized(threep_data[dt][gram]));
                qassert(fn_gram_list[dt][gram] != "");
                avg = threep_data[dt][gram];
                write_field_float_from_double(avg, fn_gram_list[dt][gram]);
            }
        }
    }

    inline void compute_proton_three_point_func_smear(const std::string &job_tag,
                                                      const int traj)
    {
        check_sigterm();
        Timer::autodisplay();
        TIMER_VERBOSE("compute_proton_three_point_func_smear");
        const std::string path = get_proton_three_point_func_smear_path(job_tag, traj);
        const std::string path_checkpoint = path + "/checkpoint.txt";
        if (does_file_exist_sync_node(path_checkpoint))
        {
            return;
        }
        if (does_file_exist_sync_node(path + "/threep-data-dt10-5.field"))
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
                ssprintf("lock-proton-threep-smear-%s-%d", job_tag.c_str(), traj)))
        {
            return;
        }
        setup(job_tag, traj);
        qmkdir_info("analysis/proton_threep_smear");
        qmkdir_info(ssprintf("analysis/proton_threep_smear/%s", job_tag.c_str()));
        qmkdir_info(path);
        std::vector<int> gram_convert_list;
        gram_convert_list.push_back(1);
        gram_convert_list.push_back(2);
        gram_convert_list.push_back(3);
        gram_convert_list.push_back(4);
        gram_convert_list.push_back(5);

        compute_proton_three_point_func_smear_type(job_tag, traj, gram_convert_list);
        release_lock();
    }
}
