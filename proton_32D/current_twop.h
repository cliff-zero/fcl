#pragma once

#include "data-load.h"
#include "contract_proton.h"
#include "compute-utils.h"

namespace qlat
{
    inline std::string get_va_two_point_func_path(const std::string &job_tag,
                                                  const int traj)
    {
        return ssprintf("analysis/va_twop/%s/results=%d", job_tag.c_str(),
                        traj);
    }

    inline void compute_va_two_point_func_type(const std::string &job_tag, const int traj)
    {
        check_sigterm();
        check_time_limit();
        Timer::autodisplay();
        const std::string path = get_va_two_point_func_path(job_tag, traj);

        const std::string path_va_two_point = path + ssprintf("/va-twop-data.field");

        bool is_complete = true;

        if (not is_d_field(path_va_two_point))
        {
            is_complete = false;
        }

        if (is_complete)
        {
            return;
        }
        TIMER_VERBOSE("compute_va_two_point_func_type");
        const PointSelection &psel = get_point_selection(job_tag, traj);
        const long n_points = psel.size();
        const FieldSelection &fsel = get_field_selection(job_tag, traj);
        FieldM<Complex, 8 * 8> data_va_twop;
        set_zero(data_va_twop);
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
            TIMER_VERBOSE("compute_va_two_point_func_type-iter");
            iter += 1;
            const ShiftShufflePlan ssp = make_shift_shuffle_plan(fsel, -xg_y);

            displayln_info(fname + ssprintf(":n=%ld iter=%ld", n,
                                            iter));
            contract_va_two_point_acc(data_va_twop, job_tag, traj, xg_y, xg_y_psel_idx, psel,
                                      fsel, ssp);
        }
        const long n_iter = iter;

        const double coef = 1.0 / (double)n_iter;
        data_va_twop *= coef;

        qassert(is_initialized(data_va_twop));

        write_field_float_from_double(data_va_twop, path_va_two_point);
    }

    inline void compute_va_two_point_func(const std::string &job_tag,
                                          const int traj)
    {
        check_sigterm();
        Timer::autodisplay();
        TIMER_VERBOSE("compute_va_two_point_func");
        const std::string path = get_va_two_point_func_path(job_tag, traj);
        const std::string path_checkpoint = path + "/checkpoint.txt";
        if (does_file_exist_sync_node(path_checkpoint))
        {
            return;
        }
        if (does_file_exist_sync_node(path + "/va-twop-data.field"))
        {
            return;
        }
        if (not(check_wall_src_props(job_tag, traj, 0) and
                check_prop_psrc(job_tag, traj, 0)))
        {
            return;
        }
        check_sigterm();
        check_time_limit();
        if (not obtain_lock(
                ssprintf("lock-va-twop-%s-%d", job_tag.c_str(), traj)))
        {
            return;
        }
        setup(job_tag, traj);
        qmkdir_info("analysis/va_twop");
        qmkdir_info(ssprintf("analysis/va_twop/%s", job_tag.c_str()));
        qmkdir_info(path);
        compute_va_two_point_func_type(job_tag, traj);
        release_lock();
    }
}