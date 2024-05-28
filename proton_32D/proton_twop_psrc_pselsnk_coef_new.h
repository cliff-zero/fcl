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

    inline LatData mk_proton_two_point_sum_table_new(const Coordinate &total_site)
    {
        LatData ld;
        ld.info.push_back(lat_dim_number("tsep", 0, total_site[3] - 1));
        ld.info.push_back(lat_dim_re_im());
        lat_data_alloc(ld);
        set_zero(ld);
        return ld;
    }

    inline void contract_proton_twop_block_new(std::vector<std::vector<std::vector<ProtonTwopBlock>>> &block,
                                               std::vector<int> &psel_num_list,
                                               std::vector<std::vector<double>> &psel_dt_num_list,
                                               const std::string &job_tag, const int traj,
                                               const Geometry &geo, const PointSelection &psel,
                                               const int dtmax, const int dtmin,
                                               const LatData &coef_table)
    {
        TIMER_VERBOSE("contract_proton_twop_block_new");
        const int num_dt = dtmax - dtmin + 1;
        qassert(num_dt > 0);
        const long n_points = psel.size();
        const Coordinate total_site = geo.total_site();
        const double dx_max = sqrt(32.0 * 32.0 + 12.0 * 12.0 * 3.0);

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
                const Coordinate xg_sep1 = smod(xg_snk - xg_src, total_site);
                const int dist1 = ceil(sqrt(sum(xg_sep1 * xg_sep1)) / dx_max * 80.0);
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
                    const Coordinate xg_sep2 = smod(xg_y - xg_src, total_site);
                    const Coordinate xg_sep3 = smod(xg_snk - xg_y, total_site);
                    const int dist2 = ceil(sqrt(sum(xg_sep2 * xg_sep2)) / dx_max * 80.0);
                    const int dist3 = ceil(sqrt(sum(xg_sep3 * xg_sep3)) / dx_max * 80.0);
                    const double coef_sphere = real(lat_data_complex_get_const(coef_table,
                                                                               make_array(dist1, dist2, dist3))[0]);
                    psel_dt_num_list[dt][dt_y_src] += coef_sphere;
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

    inline LatData proton_twop_psnk_sum(
        const std::string &job_tag, const int traj,
        const Coordinate &xg_y,
        const long xg_y_psel_idx,
        const PointSelection &psel, const FieldSelection &fsel,
        const std::vector<int> &psel_num_list,
        const std::vector<std::vector<double>> &psel_dt_num_list,
        const int tsep_snk, const int tsep_src, const int dtmin,
        const std::vector<std::vector<std::vector<ProtonTwopBlock>>> &block,
        const LatData &coef_table)
    {
        TIMER_VERBOSE("contract_proton_fourp_unshifted_psel_typeII");
        const Geometry &geo = fsel.f_rank.geo();
        const Coordinate &total_site = geo.total_site();
        const int num_dtxy = get_proton_fourp_max_sep(total_site[3]);
        const int yt = xg_y[3];
        const long n_points = psel.size();
        const double dx_max = sqrt(32.0 * 32.0 + 12.0 * 12.0 * 3.0);
        std::vector<Complex> c_twop(total_site[3]);
        set_zero(c_twop);

        for (int dt = 0; dt < num_dtxy; ++dt)
        {
            const int dtxy = dt - num_dtxy / 2;
            const int i_dt = abs(dtxy) + tsep_snk + tsep_src - dtmin;
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

            const int t_src_snk = mod(tsnk - tsrc, total_site[3]);
            qassert(t_src_snk - dtmin == i_dt);
            const int dt_y = mod(yt - tsrc, total_site[3]);

            qassert(psel_dt_num_list[i_dt][dt_y] != 0);
            Complex coef_psel = 1.0 / psel_dt_num_list[i_dt][dt_y];

            if (tsrc > tsnk)
            {
                coef_psel *= -1.0;
            }

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
                    const Coordinate xg_sep1 = smod(xg_snk - xg_src, total_site);
                    const Coordinate xg_sep2 = smod(xg_y - xg_src, total_site);
                    const Coordinate xg_sep3 = smod(xg_snk - xg_y, total_site);
                    const int dist1 = ceil(sqrt(sum(xg_sep1 * xg_sep1)) / dx_max * 80.0);
                    const int dist2 = ceil(sqrt(sum(xg_sep2 * xg_sep2)) / dx_max * 80.0);
                    const int dist3 = ceil(sqrt(sum(xg_sep3 * xg_sep3)) / dx_max * 80.0);
                    const double coef_sphere = real(lat_data_complex_get_const(coef_table,
                                                                               make_array(dist1, dist2, dist3))[0]);
                    const Complex coef = coef_psel * coef_sphere;
                    const ProtonTwopBlock &tpblock = block[xg_snk_psel_idx][i_dt][idx_src];
                    qassert(tpblock.is_build);
                    qassert(tpblock.xg_snk_psel_idx == xg_snk_psel_idx);
                    qassert(tpblock.xg_src_psel_idx == xg_src_psel_idx);
                    c_twop[t_src_snk] += tpblock.c_twop * coef;

                    idx_src += 1;
                }
            }
        }
        LatData ld = mk_proton_two_point_sum_table_new(total_site);
        set_zero(ld);
        Vector<Complex> c_sum = lat_data_cget(ld);
        for (int tss = 0; tss < total_site[3]; ++tss)
        {
            c_sum[tss] = c_twop[tss];
        }
        return ld;
    }

    inline std::string get_proton_twop_psrc_pselsnk_path_coef_new(const std::string &job_tag,
                                                                  const int traj)
    {
        return ssprintf("analysis/proton_twop_pselsnk_coef_new/%s/results=%d", job_tag.c_str(),
                        traj);
    }

    inline void compute_proton_two_point_psnk_func_coef_new(const std::string &job_tag, const int traj, const int tsep_snk, const int tsep_src)
    {
        check_sigterm();
        check_time_limit();
        Timer::autodisplay();
        const std::string path = get_proton_twop_psrc_pselsnk_path_coef_new(job_tag, traj);
        const std::string path_two_point_psrc_sum =
            path + ssprintf("/two-point-psrc-pselsnk-sum-dt%d-%d.lat", tsep_snk, tsep_src);
        LatData ld_two_point_psrc_func_sum;

        TIMER_VERBOSE("compute_proton_two_point_psnk_func_coef_new");
        const PointSelection &psel = get_point_selection(job_tag, traj);
        const long n_points = psel.size();
        const FieldSelection &fsel = get_field_selection(job_tag, traj);
        const Geometry &geo = fsel.f_rank.geo();
        const Coordinate total_site = geo.total_site();
        const int num_dtxy = get_proton_fourp_max_sep(total_site[3]);
        const int dtmax = num_dtxy / 2 + tsep_snk + tsep_src;
        const int dtmin = tsep_snk + tsep_src;

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
        contract_proton_twop_block_new(block, psel_num_list, psel_dt_num_list, job_tag, traj, geo, psel, dtmax, dtmin, coef_fourp);

        long iter = 0;
        for (long n = 0; n < n_points; ++n)
        {
            const long xg_y_psel_idx = n;
            const Coordinate &xg_y = psel[xg_y_psel_idx];
            if (get_point_src_info(job_tag, traj, xg_y, 0).size() == 0)
            {
                continue;
            }
            Timer::autodisplay();
            TIMER_VERBOSE("compute_proton_two_point_psnk_func_coef_new-iter");
            iter += 1;

            displayln_info(fname + ssprintf(":n=%ld iter=%ld", n,
                                            iter));

            LatData ld = proton_twop_psnk_sum(job_tag, traj, xg_y,
                                              xg_y_psel_idx, psel, fsel, psel_num_list,
                                              psel_dt_num_list, tsep_snk, tsep_src, dtmin, block, coef_fourp);
            ld_two_point_psrc_func_sum += ld;
        }

        const double coef = (double)total_site[0] * (double)total_site[1] * (double)total_site[2];
        ld_two_point_psrc_func_sum *= coef;

        lat_data_save_info(path_two_point_psrc_sum, ld_two_point_psrc_func_sum);
    }

    inline void compute_proton_two_point_psrc_pselsnk_func_new(const std::string &job_tag,
                                                               const int traj)
    {
        check_sigterm();
        Timer::autodisplay();
        TIMER_VERBOSE("compute_proton_two_point_psrc_pselsnk_func_new");
        const int tsep_snk = 1;
        const int tsep_src = 1;
        const std::string path = get_proton_twop_psrc_pselsnk_path_coef_new(job_tag, traj);
        const std::string path_checkpoint = path + "/checkpoint.txt";
        if (does_file_exist_sync_node(path_checkpoint))
        {
            return;
        }
        if (does_file_exist_sync_node(path + ssprintf("/two-point-psrc-pselsnk-sum-dt%d-%d.lat", tsep_snk, tsep_src)))
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
                ssprintf("lock-two-point-pselsnk-func-coef-new-%s-%d", job_tag.c_str(), traj)))
        {
            return;
        }
        setup(job_tag, traj);
        qmkdir_info("analysis/proton_twop_pselsnk_coef_new");
        qmkdir_info(ssprintf("analysis/proton_twop_pselsnk_coef_new/%s", job_tag.c_str()));
        qmkdir_info(path);
        compute_proton_two_point_psnk_func_coef_new(job_tag, traj, tsep_snk, tsep_src);
        release_lock();
    }
}
