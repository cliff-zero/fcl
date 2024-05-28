#include "proton_fourp_disconnected_smear.h"
#include "proton_twop_psrc_pselsnk_pz.h"
#include "proton_fourp_psrc_split.h"
#include "proton_threep_smear.h"
#include "proton_fourp_smear_9to11.h"
#include "current_twop_sp.h"
#include "contract_proton.h"

namespace qlat
{ //

  inline void compute_traj(const std::string &job_tag, const int traj)
  {
    setup(job_tag);
    TIMER_VERBOSE("compute_traj");
    displayln_info(ssprintf("Computing %s %d", job_tag.c_str(), traj));
    compute_proton_two_point_psrc_pselsnk_func(job_tag, traj);
    //compute_pion_two_point_psrc_func(job_tag, traj);
    //compute_proton_four_point_func_smear(job_tag, traj);
    //  compute_proton_four_point_split(job_tag, traj);
    // compute_proton_three_point_func_smear(job_tag, traj);
    //  compute_proton_four_point_func_psel(job_tag, traj);
    //compute_proton_four_point_disconnected_func_coef(job_tag, traj);
    //compute_proton_two_point_psrc_pselsnk_func(job_tag, traj);
    //  compute_proton_two_point_psrc_pselsnk_func_new(job_tag, traj);
    //compute_proton_two_point_va_func_noama(job_tag, traj);
    //compute_pion_2pt_smear_main(job_tag, traj);
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
    compute_traj("24D", 1010);
    compute_traj("24D", 1900);
    compute_traj("24D", 1010);
    compute_traj("24D", 1900);
    Timer::display();
    sync_node();
    ssleep(3.0);
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
  //
  // test();
  //
  std::vector<std::string> job_tags;
  // SADJUST ME
  job_tags.push_back("32Dfine");
  //
  for (int k = 0; k < (int)job_tags.size(); ++k)
  {
    const std::string &job_tag = job_tags[k];
    compute(job_tag);
  }
  end();
  return 0;
}
