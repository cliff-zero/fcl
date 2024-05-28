## 一些时间变量定义
```c++
    const int num_dtxy = get_proton_baryon_pi_max_sep(total_site[3]);
    const int dtmax = num_dtxy / 2 + 2 * tsep_list[tsep_num - 1];
    const int dtmin = 2 * tsep_list[0];
    const int num_dt = dtmax - dtmin + 1;
    const int num_dxxy = num_dtxy + 8;
    psel_num_list //每个时间片上有多少point-selection的点
    psel_dt_num_list // 每种时间长度的block，每个y位置，有多少point-selection点位置。

    inline LatData mk_baryon_pi_four_point(const std::vector<int> &gram_convert_list,const Coordinate &total_site)
    {
        LatData ld;
        ld.info.push_back(lat_dim_re_im("txy", 0, get_baryon_pi_max_sep(total_site[3])));
        ld.info.push_back(lat_dim_number("gram", 0, gram_convert_list.size() - 1));
        lat_data_alloc(ld);
        set_zero(ld);
        return ld;
    }
```

## block性质
 snk端是全空间1024个点都遍历，然后找出对应时间间隔src点

 ## 传播子读取
 ```c++
    const PselProp &prop_src = get_psel_prop_smear(job_tag, traj, xg_src, 0, 0, 1);
    // 函数里面先读src指标再读snk指标
 ```

 ## 奇怪的顺序
 `contract_proton_pp_y_block`对于y的遍历是在函数外面的
 有可能是为了尽量少读取大传播子p-f