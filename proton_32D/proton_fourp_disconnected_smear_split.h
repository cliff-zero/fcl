#pragma once

#include "data-load.h"
#include "contract_proton.h"
#include "compute-utils.h"

namespace qlat
{
#ifndef PROTONTWOPBLOCK
#define PROTONTWOPBLOCK
    struct ProtonTwopBlock
    {
        bool is_build;
        Complex c_twop;
        long xg_snk_psel_idx;
        long xg_src_psel_idx;
        int n_snk;
        int n_src;

        ProtonTwopBlock()
        {
            is_build = false;
            xg_snk_psel_idx = -1;
            xg_src_psel_idx = -1;
            n_snk = -1;
            n_src = -1;
            c_twop = 0.0;
        }
    };
#endif

    inline void contract_proton_twop_block(std::vector<std::vector<std::vector<ProtonTwopBlock>>> &block,
                                           std::vector<int> &psel_num_list,
                                           std::vector<std::vector<double>> &psel_dt_num_list,
                                           const std::string &job_tag, const int traj,
                                           const Geometry &geo, const PointSelection &psel, const PointSelection &psel_smear,
                                           const int dtmax, const int dtmin)
    {
        TIMER_VERBOSE("contract_proton_pp_block");
        const int num_dt = dtmax - dtmin + 1;
        qassert(num_dt > 0);
        const long n_points = psel_smear.size();
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
            const Coordinate xg = psel_smear[xg_psel_idx];
            const int t = xg[3];
            list_n_from_idx[xg_psel_idx] = psel_num_list[t];
            psel_num_list[t] += 1;
        }
        for (long nsnk = 0; nsnk < n_points; ++nsnk)
        {
            const long xg_snk_psel_idx = nsnk;
            const Coordinate &xg_snk = psel_smear[xg_snk_psel_idx];
            const int tsnk = xg_snk[3];
            for (int dt = 0; dt < num_dt; ++dt)
            {
                const int del_t = dtmin + dt;
                const int tsrc = mod(tsnk - del_t, total_site[3]);
                block[xg_snk_psel_idx][dt].resize(psel_num_list[tsrc]);
            }
        }

        for (long nsrc = 0; nsrc < n_points; ++nsrc)
        {
            const long xg_src_psel_idx = nsrc;
            const Coordinate &xg_src = psel_smear[xg_src_psel_idx];
            const int tsrc = xg_src[3];
            for (long nsnk = 0; nsnk < n_points; ++nsnk)
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
                    const int dt_snk_y = mod(tsnk - ty, total_site[3]);
                    if ((dt_y_src < 0) || (dt_y_src > del_t))
                    {
                        continue;
                    }
                    psel_dt_num_list[dt][dt_y_src] += 1;
                }
            }
        }

        for (long nsrc = 0; nsrc < n_points; ++nsrc)
        {
            const long xg_src_psel_idx = nsrc;
            const Coordinate &xg_src = psel_smear[xg_src_psel_idx];
            const int tsrc = xg_src[3];
            const PselProp &prop_src = get_psel_prop_smear(job_tag, traj, xg_src, 0, 0, 1);

#pragma omp parallel for
            for (long nsnk = 0; nsnk < n_points; ++nsnk)
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
                ProtonTwopBlock tpblock;
                const WilsonMatrix &wm = prop_src.get_elem(xg_snk_psel_idx);
                tpblock.c_twop -= proton_twop_block(wm, wm, wm, 0);
                tpblock.c_twop += proton_twop_block(wm, wm, wm, 1);
                tpblock.xg_snk_psel_idx = xg_snk_psel_idx;
                tpblock.xg_src_psel_idx = xg_src_psel_idx;
                tpblock.n_snk = list_n_from_idx[xg_snk_psel_idx];
                tpblock.n_src = list_n_from_idx[xg_src_psel_idx];
                tpblock.is_build = true;

                const int n_t = list_n_from_idx[xg_src_psel_idx];
                block[xg_snk_psel_idx][dt][n_t] = tpblock;
            }
        }
        return;
    }

    inline void contract_proton_fourp_disc_unshifted_psel_acc_x_coef(Vector<Complex> vv, const Complex coef,
                                                                     const WilsonMatrix &wm_x_y,
                                                                     const Complex &c_twop)
    {
        const array<SpinMatrix, 8> &va_ms = get_va_matrices();
        const SpinMatrix &gamma5 = SpinMatrixConstants::get_gamma5();
        const WilsonMatrix wm_y_x =
            gamma5 * (WilsonMatrix)matrix_adjoint(wm_x_y) * gamma5;

        const WilsonMatrix wm_pp_y_tsrc_tsnk_x = c_twop * wm_y_x;
        for (int mu = 0; mu < 8; mu++)
        {
            const WilsonMatrix wm_pp = wm_x_y * va_ms[mu] * wm_pp_y_tsrc_tsnk_x;
            for (int nu = 0; nu < 8; nu++)
            {
                const int mu_nu = 8 * nu + mu;
                vv[mu_nu] += coef * matrix_trace(wm_pp, va_ms[nu]);
            }
        }

        return;
    }

    inline void contract_proton_fourp_disc_unshifted_psel_coef(
        SelectedField<Complex> &proton_fourp,
        const std::string &job_tag, const int traj,
        const SelProp &prop_x_y,
        const Coordinate &xg_y,
        const long xg_y_psel_idx, const int tsep_src, const int tsep_snk,
        const PointSelection &psel, const FieldSelection &fsel,
        const std::vector<int> psel_num_list,
        const std::vector<std::vector<double>> psel_dt_num_list,
        const int dtmin,
        const std::vector<std::vector<std::vector<ProtonTwopBlock>>> &block)
    {
        TIMER_VERBOSE("contract_proton_fourp_unshifted_psel_typeII");
        const Geometry &geo = fsel.f_rank.geo();
        const Coordinate &total_site = geo.total_site();
        const int num_dtxy = get_proton_fourp_max_sep(total_site[3]);
        const int num_dxxy = 25;
        const int yt = xg_y[3];
        const long n_points = psel.size();
        const int multiplicity = 8 * 8;

        const int x_sparsity = 1;

        proton_fourp.init(fsel, multiplicity);
        set_zero(proton_fourp);

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
            qassert(psel_dt_num_list[i_dt][dt_y] != 0);
            const Complex coef_psel = 1.0 / psel_dt_num_list[i_dt][dt_y] * (double)x_sparsity;
            Complex c_twop = 0.0;
            for (long nsnk = 0; nsnk < n_points; ++nsnk)
            {
                const long xg_snk_psel_idx = nsnk;
                const Coordinate &xg_snk = psel[xg_snk_psel_idx];
                if (xg_snk[3] != tsnk)
                {
                    continue;
                }
                int idx_src = 0;
                for (long nsrc = 0; nsrc < n_points; ++nsrc)
                {
                    const long xg_src_psel_idx = nsrc;
                    const Coordinate &xg_src = psel[xg_src_psel_idx];
                    if (xg_src[3] != tsrc)
                    {
                        continue;
                    }
                    const Complex coef = coef_psel;
                    const ProtonTwopBlock &tpblock = block[xg_snk_psel_idx][i_dt][idx_src];
                    qassert(tpblock.is_build);
                    qassert(tpblock.xg_snk_psel_idx == xg_snk_psel_idx);
                    qassert(tpblock.xg_src_psel_idx == xg_src_psel_idx);
                    c_twop += tpblock.c_twop * coef;

                    idx_src += 1;
                }
            }
            const Complex coef = 1.0;
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
                const WilsonMatrix &wm_x_y = prop_x_y.get_elem(idx);
                Vector<Complex> vv;

                vv = proton_fourp.get_elems(idx);

                if (tsrc <= tsnk)
                {
                    contract_proton_fourp_disc_unshifted_psel_acc_x_coef(vv, coef, wm_x_y, c_twop);
                }
                else
                {
                    const Complex coef1 = -1.0 * coef;
                    contract_proton_fourp_disc_unshifted_psel_acc_x_coef(vv, coef1, wm_x_y, c_twop);
                }
            }
        }

        return;
    }

    inline void
    contract_proton_four_point_disc_acc_psel_coef(
        FieldM<Complex, 8 * 8> &fourp_data,
        const std::string &job_tag, const int traj,
        const SelProp &prop_x_y, const Coordinate &xg_y,
        const long xg_y_psel_idx, const int tsep_src, const int tsep_snk, const PointSelection &psel,
        const FieldSelection &fsel, const ShiftShufflePlan &ssp,
        const std::vector<int> &psel_num_list,
        const std::vector<std::vector<double>> &psel_dt_num_list, const int dtmin,
        const std::vector<std::vector<std::vector<ProtonTwopBlock>>> &block)
    // xg_y = psel[xg_y_psel_idx] is the point src location for prop3_x_y
    // ssp = make_shift_shuffle_plan(fsel, -xg_y);
    {
        TIMER_VERBOSE("contract_proton_four_point_acc_psel");
        const Geometry &geo = fsel.f_rank.geo();
        qassert(geo == prop_x_y.geo());
        qassert(fsel.n_elems == prop_x_y.n_elems);
        qassert(is_initialized(prop_x_y));
        qassert(ssp.shift == -xg_y);
        qassert(ssp.is_reflect == false);
        SelectedField<Complex> sfs;
        contract_proton_fourp_disc_unshifted_psel_coef(sfs, job_tag, traj, prop_x_y, xg_y, xg_y_psel_idx,
                                                       tsep_src, tsep_snk, psel, fsel, psel_num_list, psel_dt_num_list,
                                                       dtmin, block);
        qassert(fsel.prob == ssp.fsel.prob);
        const Complex coef = 1.0 / fsel.prob;

        SelectedField<Complex> &s_data = sfs;
        acc_field(fourp_data, coef, s_data, ssp);
    }

    inline std::string get_proton_four_point_disc_func_noama_path_coef(const std::string &job_tag,
                                                                       const int traj)
    {
        return ssprintf("analysis/proton_fourp_disconnected_smear_split/%s/results=%d", job_tag.c_str(),
                        traj);
    }

    inline void compute_proton_four_point_disc_func_type_coef(const std::string &job_tag, const int traj,
                                                              const std::vector<int> &tsep_list_src,
                                                              const std::vector<int> &tsep_list_snk)
    {
        check_sigterm();
        check_time_limit();
        Timer::autodisplay();
        const int tsep_num = (int)tsep_list_src.size();
        const int tsep_num_snk = (int)tsep_list_snk.size();
        qassert(tsep_num == tsep_num_snk);
        TIMER_VERBOSE("compute_proton_four_point_disc_func_type_noama");
        const PointSelection &psel = get_point_selection(job_tag, traj);
        const PointSelection &psel_smear = get_point_selection_smear(job_tag, traj);
        const long n_points = psel.size();
        const FieldSelection &fsel = get_field_selection(job_tag, traj);
        const Geometry &geo = fsel.f_rank.geo();
        const Coordinate total_site = geo.total_site();
        const std::string path = get_proton_four_point_disc_func_noama_path_coef(job_tag, traj);
        std::vector<std::string> fn_gram_list(tsep_num * total_site[3]);

        for (int ti = 0; ti < tsep_num; ++ti)
        {
            for (int ty = 0; ty < total_site[3]; ++ty)
            {
                const int &tsep_src = tsep_list_src[ti];
                const int &tsep_snk = tsep_list_snk[ti];
                fn_gram_list[ti + ty * tsep_num] = path + ssprintf("/fourp-data-dt%d-%d-disconnected-t%02d.field", tsep_snk, tsep_src, ty);
            }
        }

        bool is_complete = true;

        for (int ti = 0; ti < tsep_num; ++ti)
        {
            for (int ty = 0; ty < total_site[3]; ++ty)
            {
                if (fn_gram_list[ti + ty * tsep_num] != "")
                {
                    if (not is_d_field(fn_gram_list[ti + ty * tsep_num]))
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
        const int num_dtxy = get_proton_fourp_max_sep(total_site[3]);
        const int dtmax = num_dtxy / 2 + 2 * tsep_list_src[tsep_num - 1];
        const int dtmin = 2 * tsep_list_src[0];

        LatData coef_fourp;
        const std::string path_coef_table = ssprintf("coef-list/fourp/coef-fourp-80slice.lat");
        if (does_file_exist_sync_node(path_coef_table))
        {
            coef_fourp = lat_data_load_info(path_coef_table);
        }
        else
        {
            qassert(false);
        }

        std::vector<int> psel_num_list;
        std::vector<std::vector<double>> psel_dt_num_list;
        std::vector<std::vector<std::vector<ProtonTwopBlock>>> block;
        contract_proton_twop_block(block, psel_num_list, psel_dt_num_list, job_tag, traj, geo, psel, psel_smear, dtmax, dtmin);
        for (int ty = 0; ty < total_site[3]; ++ty)
        {
            std::vector<FieldM<Complex, 8 * 8>> fourp_data(tsep_num);

            long iter = 0;
            for (long n = 0; n < n_points; ++n)
            {
                const long xg_y_psel_idx = n;
                const Coordinate &xg_y = psel[xg_y_psel_idx];
                if (ty != xg_y[3])
                {
                    continue;
                }
                if (get_point_src_info(job_tag, traj, xg_y, 0).size() == 0)
                {
                    continue;
                }
                Timer::autodisplay();
                TIMER_VERBOSE("compute_proton_four_point_disc_func_type_noama-iter");
                iter += 1;
                const ShiftShufflePlan ssp = make_shift_shuffle_plan(fsel, -xg_y);
                const SelProp &prop_x_y = get_prop_psrc(job_tag, traj, xg_y, 0, 0);

                displayln_info(fname + ssprintf(":n=%ld iter=%ld", n,
                                                iter));
                for (int ti = 0; ti < tsep_num; ++ti)
                {
                    const int tsep_src = tsep_list_src[ti];
                    const int tsep_snk = tsep_list_snk[ti];
                    contract_proton_four_point_disc_acc_psel_coef(fourp_data[ti], job_tag, traj, prop_x_y, xg_y, xg_y_psel_idx,
                                                                  tsep_src, tsep_snk, psel_smear, fsel, ssp, psel_num_list, psel_dt_num_list,
                                                                  dtmin, block);
                }
            }

            for (int ti = 0; ti < tsep_num; ++ti)
            {
                const double coef = (double)total_site[0] * (double)total_site[1] * (double)total_site[2];
                fourp_data[ti] *= (coef * coef);
            }

            FieldM<Complex, 8 * 8> avg;

            for (int ti = 0; ti < tsep_num; ++ti)
            {
                qassert(is_initialized(fourp_data[ti]));
                qassert(fn_gram_list[ti] != "");
                avg = fourp_data[ti];
                write_field_float_from_double(avg, fn_gram_list[ti + ty * tsep_num]);
            }
        }
    }

    inline void compute_proton_four_point_disconnected_func_coef(const std::string &job_tag,
                                                                 const int traj)
    {
        check_sigterm();
        Timer::autodisplay();
        TIMER_VERBOSE("compute_proton_four_point_disconnected_func_coef");
        const std::string path = get_proton_four_point_disc_func_noama_path_coef(job_tag, traj);
        const std::string path_checkpoint = path + "/checkpoint.txt";
        if (does_file_exist_sync_node(path_checkpoint))
        {
            return;
        }
        if (does_file_exist_sync_node(path + "/fourp-data-dt1-1-disconnected-t63.field"))
        {
            return;
        }
        if (not(check_prop_psrc(job_tag, traj, 0) and check_prop_smear(job_tag, traj, 0)))
        {
            return;
        }
        check_sigterm();
        check_time_limit();
        if (not obtain_lock(
                ssprintf("lock-proton-fourp-disc-smear-split-%s-%d", job_tag.c_str(), traj)))
        {
            return;
        }
        setup(job_tag, traj);
        qmkdir_info("analysis/proton_fourp_disconnected_smear_split");
        qmkdir_info(ssprintf("analysis/proton_fourp_disconnected_smear_split/%s", job_tag.c_str()));
        qmkdir_info(path);
        std::vector<int> tsep_list_src, tsep_list_snk;
        tsep_list_src.push_back(1); // small to large
        //tsep_list_src.push_back(3); // small to large
        tsep_list_snk.push_back(1);
        //tsep_list_snk.push_back(2);
        compute_proton_four_point_disc_func_type_coef(job_tag, traj, tsep_list_src, tsep_list_snk);
        release_lock();
    }
}
