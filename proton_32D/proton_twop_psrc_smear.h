#pragma once

#include "data-load.h"
#include "contract_proton.h"
#include "compute-utils.h"

namespace qlat
{ //

    inline std::string get_proton_two_point_psrc_path(const std::string &job_tag,
                                                      const int traj)
    {
        return ssprintf("analysis/proton_twop_smear_fsel/%s/results=%d", job_tag.c_str(),
                        traj);
    }

    inline void compute_proton_two_point_psrc(const std::string &job_tag,
                                              const int traj)
    {
        check_sigterm();
        Timer::autodisplay();
        TIMER_VERBOSE("compute_proton_two_point_psrc");
        const PointSelection &psel = get_point_selection_smear(job_tag, traj);
        const long n_points = psel.size();
        const FieldSelection &fsel = get_field_selection(job_tag, traj);
        LatData ld_two_point_psrc_func;
        const std::string path = get_proton_two_point_psrc_path(job_tag, traj);
        const std::string path_two_point_psrc =
            path + ssprintf("/two-point-psrc.lat");

        long iter = 0;
        for (long n = 0; n < n_points; n++)
        {
            const long xg_src_psel_idx = n;
            const Coordinate &xg_src = psel[xg_src_psel_idx];
            const int &tslice = xg_src[3];
            iter += 1;
            const SelProp &pprop = get_prop_smear(job_tag, traj, xg_src, 0, 0, 0);
            LatData ld = contract_proton_two_point_psrc_function(pprop, tslice, fsel);
            ld_two_point_psrc_func += ld;
        }
        const long n_iter = iter;
        const double coef = 1.0 / (double)n_iter;
        ld_two_point_psrc_func *= coef;
        displayln_info(ssprintf("n_iter=%ld", n_iter));
        lat_data_save_info(path_two_point_psrc, ld_two_point_psrc_func);
    }

    inline void compute_proton_two_point_psrc_func(const std::string &job_tag, const int traj)
    {
        check_sigterm();
        Timer::autodisplay();
        TIMER_VERBOSE("compute_proton_two_point_psrc_func");
        const std::string path = get_proton_two_point_psrc_path(job_tag, traj);
        const std::string path_checkpoint = path + "/checkpoint.txt";
        if (does_file_exist_sync_node(path_checkpoint))
        {
            return;
        }
        if (does_file_exist_sync_node(
                path + "/two-point-psrc.lat"))
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
                ssprintf("lock-two-point-smear-fsel-func-%s-%d", job_tag.c_str(), traj)))
        {
            return;
        }
        setup(job_tag, traj);
        qmkdir_info("analysis/proton_twop_smear_fsel");
        qmkdir_info(ssprintf("analysis/proton_twop_smear_fsel/%s", job_tag.c_str()));
        qmkdir_info(path);
        compute_proton_two_point_psrc(job_tag, traj);
        release_lock();
    }

} // namespace qlat
