# rsync -vaz --delete --exclude "transcation" --exclude "build/*" --exclude "analysis/*"  --exclude "lock*" --exclude "analysis_old/*" --exclude "core.*" --exclude "out/*" --exclude "err/*" --exclude "results/*" ~/Yusheng_Gao/Contraction/baryon_two_point_24D ~/Yusheng_Gao/Contraction/baryon_contraction_Npi_2pt_24D

rsync -vaz \
 --exclude ".git/*" \
 --exclude "analysis/*" \
 --exclude "analysis_old/*" \
 --exclude "build/*" \
 --exclude "build_debug/*" \
 --exclude "out/*" \
 --exclude "err/*" \
 --exclude "lock*" \
 --exclude "results/*" \
 --exclude "compute-utils.h" \
 --exclude "configs.h*" \
 --exclude "data-load-base.h" \
 --exclude "data-load.h" \
 --exclude "data-paths.h" \
 --exclude "main.cpp" \
 --exclude "psrc-distribution.h" \
 --exclude "psrc-sample.h" \
 --exclude "run*" \
 --exclude "setenv.sh" \
 --exclude "th-path.h" \
 --exclude "to_24D.sh" \
 ~/Yusheng_Gao/Contraction/baryon_contraction_Npi_2pt_32Dfine/ ~/Yusheng_Gao/Contraction/baryon_contraction_Npi_2pt_24D/
