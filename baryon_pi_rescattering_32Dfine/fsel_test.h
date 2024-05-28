#pragma once

#include "data-load.h"
// #include "baryon-pion-block-acc.h"
// #include "baryon-pion-coefficient.h"
#include "compute-utils.h"

namespace qlat
{
    inline void fsel_prop_test(const std::string &job_tag, const int &traj)
    {
        check_sigterm();
        check_time_limit();
        Timer::autodisplay();
       
        TIMER_VERBOSE("compute_baryon_pi_rescattering_fsel_arbitrarty_func_psel_type");
        const PointsSelection &psel = get_point_selection(job_tag, traj);
        const long n_points = psel.size();
        const PointsSelection &psel_smear = get_point_selection_smear(job_tag, traj);
        const long n_points_smear = psel_smear.size();

        const FieldSelection &fsel = get_field_selection(job_tag, traj);

        const Geometry &geo = fsel.get_geo();
        // const Coordinate total_site = geo.total_site();

    for(int i =0 ;i< n_points; i++)
        {
        const Coordinate &xg = psel[i];
        const SelProp &prop = get_prop_psrc(job_tag, traj, xg, 0, 0);

        // const Coordinate &xg_smear = psel_smear[i];
        // const SelProp &prop_smear = get_prop_smear(job_tag, traj, xg_smear, 0, 0, 0);
        // const SelProp &prop_smear = get_prop_smear(job_tag, traj, xg, 0, 0, 0);
        }
    }

    inline void compute_fsel_prop_test(const std::string &job_tag, const int &traj)
    {
        check_sigterm();
        Timer::autodisplay();
        TIMER_VERBOSE("compute_baryon_pi_rescattering_fsel_arbitrarty_func_psel");

        displayln_info(get_prop_smear_path(job_tag, traj, 0));
        displayln_info(get_psel_prop_smear_path(job_tag, traj, 0));
        if (not(check_prop_smear(job_tag, traj, 0)))
        {
            displayln_info("err_check_prop_smear;");
            displayln_info(get_prop_smear_path(job_tag, traj, 0));
            displayln_info(get_psel_prop_smear_path(job_tag, traj, 0));
            return;
        }
        check_sigterm();
        check_time_limit();
        if (not obtain_lock(
                ssprintf("lock-baryon-pi-%s-%d", job_tag.c_str(), traj)))
        {
            return;
        }
        setup(job_tag, traj);

        fsel_prop_test(job_tag, traj);

        release_lock();
    }
}
