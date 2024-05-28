#pragma once

#include <qlat/qlat.h>
#include "th-path.h"

namespace qlat
{ //

  inline std::string get_point_src_info_path(const std::string &job_tag, const int traj)
  {
    return prop_path_default + ssprintf("point-src-info/%s/traj=%d.txt", job_tag.c_str(), traj);
  }

  inline std::string get_point_selection_path(const std::string &job_tag, const int traj)
  {
    return prop_path_default + ssprintf("point-selection/%s/traj=%d.txt", job_tag.c_str(), traj);
  }

  inline std::string get_point_selection_smear_path(const std::string &job_tag, const int traj)
  {
    return prop_path_default + ssprintf("point-selection-smear/%s/traj-%d.txt", job_tag.c_str(), traj);
  }

  inline std::string get_point_src_info_simulated_path(const std::string &job_tag, const int traj)
  {
    return prop_path_default + ssprintf("point-src-info-simulated/%s/traj=%d.txt", job_tag.c_str(), traj);
  }

  inline std::string get_point_distribution_path(const std::string &job_tag)
  {
    return prop_path_default + ssprintf("point-distribution/%s/point-distribution.txt", job_tag.c_str());
  }

  inline std::string get_field_selection_path(const std::string &job_tag, const int traj)
  {
    return prop_path_default + ssprintf("field-selection/%s/traj=%d.field", job_tag.c_str(), traj);
  }

  inline std::string get_gauge_transform_path(const std::string &job_tag, const int traj)
  {
    return prop_path_default + ssprintf("gauge-transform/%s/traj=%d.field", job_tag.c_str(), traj);
  }

  inline std::string get_prop_psrc_light_path(const std::string &job_tag, const int traj)
  {
    return prop_path_default + ssprintf("prop-psrc-light/%s/traj=%d", job_tag.c_str(), traj);
  }

  inline std::string get_prop_smear_light_path(const std::string &job_tag, const int traj)
  {
    return prop_path_default + ssprintf("prop-smear-light/%s/traj-%d", job_tag.c_str(), traj);
  }

  inline std::string get_psel_prop_psrc_light_path(const std::string &job_tag, const int traj)
  {
    return prop_path_default + ssprintf("psel-prop-psrc-light/%s/traj=%d", job_tag.c_str(), traj);
  }

  inline std::string get_psel_prop_smear_light_path(const std::string &job_tag, const int traj)
  {
    return prop_path_default + ssprintf("psel-prop-smear-light/%s/traj-%d", job_tag.c_str(), traj);
  }

  inline std::string get_prop_psrc_strange_path(const std::string &job_tag, const int traj)
  {
    return prop_path_default + ssprintf("prop-psrc-strange/%s/traj=%d", job_tag.c_str(), traj);
  }

  inline std::string get_prop_smear_strange_path(const std::string &job_tag, const int traj)
  {
    return prop_path_default + ssprintf("prop-smear-strange/%s/traj=%d", job_tag.c_str(), traj);
  }

  inline std::string get_psel_prop_psrc_strange_path(const std::string &job_tag, const int traj)
  {
    return prop_path_default + ssprintf("psel-prop-psrc-strange/%s/traj=%d", job_tag.c_str(), traj);
  }

  inline std::string get_psel_prop_smear_strange_path(const std::string &job_tag, const int traj)
  {
    return prop_path_default + ssprintf("psel-prop-smear-strange/%s/traj=%d", job_tag.c_str(), traj);
  }

  inline std::string get_prop_psrc_exact_path(const std::string &job_tag, const int traj)
  {
    return prop_path_default + ssprintf("prop-psrc-exact/%s/traj=%d", job_tag.c_str(), traj);
  }

  inline std::string get_psel_prop_psrc_exact_path(const std::string &job_tag, const int traj)
  {
    return prop_path_default + ssprintf("psel-prop-psrc-exact/%s/traj=%d", job_tag.c_str(), traj);
  }

  inline std::string get_psrc_tag(const Coordinate &xg, const int type, const int accuracy)
  {
    return ssprintf("xg=(%d,%d,%d,%d) ; type=%d ; accuracy=%d", xg[0], xg[1],
                    xg[2], xg[3], type, accuracy);
  }

  inline std::string get_smear_tag(const Coordinate &xg, const int type, const int accuracy, const int smeared)
  {
    if (smeared == 0)
    {
      return ssprintf("smear ; xg=(%d,%d,%d,%d) ; type=%d ; accuracy=%d", xg[0], xg[1],
                      xg[2], xg[3], type, accuracy);
    }
    else
    {
      return ssprintf("smear ; xg=(%d,%d,%d,%d) ; type=%d ; accuracy=%d ; smear-snk", xg[0], xg[1],
                      xg[2], xg[3], type, accuracy);
    }
  }

  inline std::string get_prop_wsrc_light_path(const std::string &job_tag, const int traj)
  {
    return prop_path_default + ssprintf("prop-wsrc-light/%s/traj=%d", job_tag.c_str(), traj);
  }

  inline std::string get_psel_prop_wsrc_light_path(const std::string &job_tag, const int traj)
  {
    return prop_path_default + ssprintf("psel-prop-wsrc-light/%s/traj=%d", job_tag.c_str(), traj);
  }

  inline std::string get_prop_wsrc_strange_path(const std::string &job_tag, const int traj)
  {
    return prop_path_default + ssprintf("prop-wsrc-strange/%s/traj=%d", job_tag.c_str(), traj);
  }

  inline std::string get_psel_prop_wsrc_strange_path(const std::string &job_tag, const int traj)
  {
    return prop_path_default + ssprintf("psel-prop-wsrc-strange/%s/traj=%d", job_tag.c_str(), traj);
  }

  inline std::string get_wsrc_tag(const int tslice, const int type, const int accuracy)
  {
    return ssprintf("tslice=%d ; type=%d ; accuracy=%d", tslice, type, accuracy);
  }

  inline std::string get_prop_psrc_path(const std::string &job_tag, const int traj, const int type)
  {
    if (type == 0)
    {
      return get_prop_psrc_light_path(job_tag, traj);
    }
    else if (type == 1)
    {
      return get_prop_psrc_strange_path(job_tag, traj);
    }
    else
    {
      qassert(false);
      return "";
    }
  }

  inline std::string get_prop_smear_path(const std::string &job_tag, const int traj, const int type)
  {
    if (type == 0)
    {
      return get_prop_smear_light_path(job_tag, traj);
    }
    else if (type == 1)
    {
      return get_prop_smear_strange_path(job_tag, traj);
    }
    else
    {
      qassert(false);
      return "";
    }
  }

  inline std::string get_psel_prop_psrc_path(const std::string &job_tag, const int traj, const int type)
  {
    if (type == 0)
    {
      return get_psel_prop_psrc_light_path(job_tag, traj);
    }
    else if (type == 1)
    {
      return get_psel_prop_psrc_strange_path(job_tag, traj);
    }
    else
    {
      qassert(false);
      return "";
    }
  }

  inline std::string get_psel_prop_smear_path(const std::string &job_tag, const int traj, const int type)
  {
    if (type == 0)
    {
      return get_psel_prop_smear_light_path(job_tag, traj);
    }
    else if (type == 1)
    {
      return get_psel_prop_smear_strange_path(job_tag, traj);
    }
    else
    {
      qassert(false);
      return "";
    }
  }

  inline std::string get_prop_wsrc_path(const std::string &job_tag, const int traj, const int type)
  {
    if (type == 0)
    {
      return get_prop_wsrc_light_path(job_tag, traj);
    }
    else if (type == 1)
    {
      return get_prop_wsrc_strange_path(job_tag, traj);
    }
    else
    {
      qassert(false);
      return "";
    }
  }

  inline std::string get_psel_prop_wsrc_path(const std::string &job_tag, const int traj, const int type)
  {
    if (type == 0)
    {
      return get_psel_prop_wsrc_light_path(job_tag, traj);
    }
    else if (type == 1)
    {
      return get_psel_prop_wsrc_strange_path(job_tag, traj);
    }
    else
    {
      qassert(false);
      return "";
    }
  }

} // namespace qlat
