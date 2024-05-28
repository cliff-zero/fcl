#pragma once

#include "data-load.h"
#include "contract_proton.h"
#include "compute-utils.h"

namespace qlat
{
    inline std::string get_proton_two_point_va_func_noama_path(const std::string &job_tag,
                                                               const int traj)
    {
        return ssprintf("analysis/proton_twop_va/%s/results=%d", job_tag.c_str(),
                        traj);
    }

    inline void compute_proton_two_point_va_func_type_noama(const std::string &job_tag, const int traj)
    {
        check_sigterm();
        check_time_limit();
        Timer::autodisplay();
        const std::string path = get_proton_two_point_va_func_noama_path(job_tag, traj);
        std::string fn_gram;

        fn_gram = path + ssprintf("/twop-data-va.field");

        bool is_complete = true;

        if (fn_gram != "")
        {
            if (not is_d_field(fn_gram))
            {
                is_complete = false;
            }
        }

        if (is_complete)
        {
            return;
        }
        TIMER_VERBOSE("compute_proton_two_point_va_func_type_noama");
        const PointSelection &psel = get_point_selection(job_tag, traj);
        const long n_points = psel.size();
        const FieldSelection &fsel = get_field_selection(job_tag, traj);

        FieldM<Complex, 8 * 8> twop_data;

        long iter = 0;
        for (long n = 0; n < n_points; ++n)
        {
            const long xg_y_psel_idx = n;
            const Coordinate &xg_y = psel[xg_y_psel_idx];
            Timer::autodisplay();
            TIMER_VERBOSE("compute_proton_two_point_va_func_type_noama-iter");
            iter += 1;
            const ShiftShufflePlan ssp = make_shift_shuffle_plan(fsel, -xg_y);

            displayln_info(fname + ssprintf(":n=%ld iter=%ld", n,
                                            iter));

            contract_va_two_point_acc_noama(twop_data, job_tag, traj, xg_y, xg_y_psel_idx,
                                            psel, fsel, ssp);
        }
        const Geometry &geo = fsel.f_rank.geo();
        const Coordinate total_site = geo.total_site();

        const double coef = (double)total_site[0] * (double)total_site[1] * (double)total_site[2] / (double)iter;
        twop_data *= coef;

        FieldM<Complex, 8 * 8> avg;

        qassert(is_initialized(twop_data));
        qassert(fn_gram != "");
        avg = twop_data;
        write_field_float_from_double(avg, fn_gram);
    }

    inline void compute_proton_two_point_va_func_noama(const std::string &job_tag,
                                                       const int traj)
    {
        check_sigterm();
        Timer::autodisplay();
        TIMER_VERBOSE("compute_proton_two_point_va_func_noama");
        const std::string path = get_proton_two_point_va_func_noama_path(job_tag, traj);
        const std::string path_checkpoint = path + "/checkpoint.txt";
        if (does_file_exist_sync_node(path_checkpoint))
        {
            return;
        }
        if (does_file_exist_sync_node(path + "/twop-data-va.field"))
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
                ssprintf("lock-proton-twop-va-noama-%s-%d", job_tag.c_str(), traj)))
        {
            return;
        }
        setup(job_tag, traj);
        qmkdir_info("analysis/proton_twop_va");
        qmkdir_info(ssprintf("analysis/proton_twop_va/%s", job_tag.c_str()));
        qmkdir_info(path);
        compute_proton_two_point_va_func_type_noama(job_tag, traj);
        release_lock();
    }
}
