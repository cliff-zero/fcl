// #include "baryon-pion-two-point-function.h"
// #include "baryon-pion-two-point-function-fsel.h"
// #include "baryon-pion-two-point-function-disconnect.h"
// #include "baryon-pion-two-point-function-disconnect-fsel.h"
// #include "proton_twop_smear.h"
// #include <qlat/qlat.h>  
// #include "baryon-pion-two-point-function.h"

// #include "baryon-pion-two-point-function-arbitrary.h"
// #include "baryon-pion-two-point-function-disconnect-arbitrary.h"

// #include "baryon-pion-two-point-function-arbitrary-fsel.h"
// #include "baryon-pion-two-point-function-arbitrary-fsel-MATE-ULTRA.h"
// #include "baryon-pion-two-point-function-arbitrary-fsel-MATE-ULTRA-MAGIC.h"
// #include "baryon-pion-two-point-function-arbitrary-fsel-MATE-ULTRA-MAGIC-PRO.h"
// #include "baryon-pion-two-point-function-arbitrary-fsel-disconnect-MATE-ULTRA-MAGIC-PRO.h"
// #include "fsel_test.h"
// #include "proton-pion-two-point.h"
#include "two-point-pro.h"

namespace qlat
{ //

  inline void compute_traj(const std::string &job_tag, const int traj)
  {
    setup(job_tag);
    TIMER_VERBOSE("compute_traj");
    displayln_info(ssprintf("Computing %s %d", job_tag.c_str(), traj));
    // compute_proton_two_point_smear_func(job_tag, traj);

    // compute_baryon_pi_rescattering_all_psel_arbitrarty_func_psel(job_tag, traj, true, true);
    // compute_baryon_pi_rescattering_disconnect_all_psel_arbitrarty_func_psel(job_tag, traj, true, true);

    // compute_baryon_pi_rescattering_fsel_arbitrarty_func_psel(job_tag, traj, false, false);
    // compute_baryon_pi_rescattering_fsel_arbitrarty_func_psel(job_tag, traj, true, false);
    // compute_baryon_pi_rescattering_fsel_arbitrarty_func_psel(job_tag, traj, true, true);
    // compute_baryon_pi_rescattering_fsel_disconnect_arbitrarty_func_psel(job_tag, traj, true, true);
    // compute_fsel_prop_test(job_tag, traj);

    // compute_two_point_correlation_function_all_psel_func_psel(job_tag, traj, true, true);
    compute_two_point_correlation_function_pro_all_psel_func_psel(job_tag, traj, true, true);

    // compute_baryon_pi_two_point_all_psel_func_psel(job_tag, traj, true, true);
    // compute_baryon_pi_two_point_all_psel_func_psel(job_tag, traj, true, true);
    // compute_baryon_pi_two_point_disconnect_all_psel_func_psel(job_tag, traj, true, true);

    // compute_baryon_pi_two_point_all_psel_func_psel(job_tag, traj, true, false);
    // compute_baryon_pi_two_point_disconnect_all_psel_func_psel(job_tag, traj, true, false);

    // compute_baryon_pi_two_point_all_fsel_psel_func_psel(job_tag, traj);
    // compute_baryon_pi_two_point_disconnect_fsel_psel_func_psel(job_tag, traj, true, false);
    
    clear_all_data_cache();
  }

  inline void compute(const std::string &job_tag)
  {
    TIMER_VERBOSE("compute");
    const std::vector<int> trajs = get_data_trajs(job_tag);
    for (int i = 0; i < (int)trajs.size(); ++i)
    {
      const int traj = trajs[i];
      compute_traj(job_tag, traj);
      if (get_total_time() > 1.0)
      {
        Timer::display();
      }
      Timer::reset();
    }
  }

  inline void test()
  {
    TIMER_VERBOSE("test");
    // compute_traj("32Dfine", 520);
    // compute_traj("32Dfine", 530);
    // compute_traj("32Dfine", 540);
    // compute_traj("32Dfine", 550);
    // compute_traj("32Dfine", 560);
    // compute_traj("32Dfine", 570);
    // compute_traj("32Dfine", 580);
    // compute_traj("32Dfine", 590);
    // compute_traj("32Dfine", 740);
    // compute_traj("32Dfine", 780);
    // compute_traj("32Dfine", 800);
    // compute_traj("32Dfine", 880);
    compute_traj("32Dfine", 1410);
    compute_traj("32Dfine", 680);
    compute_traj("32Dfine", 740);
    Timer::display();
    sync_node();
    // ssleep(3.0);
    exit(0);
  }

} // namespace qlat

int main(int argc, char *argv[])
{
  using namespace qlat;
  std::vector<Coordinate> size_node_list;
  size_node_list.push_back(Coordinate(1, 1, 1, 1));
  size_node_list.push_back(Coordinate(1, 1, 1, 2));
  size_node_list.push_back(Coordinate(1, 1, 1, 4));
  size_node_list.push_back(Coordinate(1, 1, 1, 8));
  size_node_list.push_back(Coordinate(1, 1, 1, 16));
  size_node_list.push_back(Coordinate(1, 1, 2, 16));
  size_node_list.push_back(Coordinate(1, 2, 2, 16));
  size_node_list.push_back(Coordinate(2, 2, 2, 16));
  size_node_list.push_back(Coordinate(2, 2, 4, 16));
  size_node_list.push_back(Coordinate(2, 4, 4, 16));
  size_node_list.push_back(Coordinate(4, 4, 4, 16));
  
  begin(&argc, &argv, size_node_list);
  display_geometry_node();
  setup_log_idx();
  setup();
  qmkdir_info(ssprintf("analysis"));

  // test();

  std::vector<std::string> job_tags;
  // SADJUST ME
  job_tags.push_back("32Dfine");
  
  for (int k = 0; k < (int)job_tags.size(); ++k)
  {
    const std::string &job_tag = job_tags[k];
    compute(job_tag);
  }
  end();
  return 0;
}
