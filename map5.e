rm fort.*

ln -s milank.dat fort.23
ln -s input$1.def fort.9
ln -s output$2.def fort.8
ln -s input$1.head fort.30
ln -s input$1.grid fort.31
ln -s input$1.diff fort.32
ln -s input$1.time fort.33
ln -s input$2.defdep fort.39
ln -s input$1.dep fort.40
ln -s input$1.h2o fort.41
ln -s input$1.coord fort.53
ln -s output$2.h2o fort.42
ln -s input$1.temp fort.43
ln -s output$2.temp fort.44
ln -s output$2.dep fort.45
ln -s out$2.time fort.34
ln -s bc$2.time fort.35
ln -s temp$2.time fort.36
ln -s tslice$2.time fort.37
ln -s tslice$2.xyz fort.38
ln -s outline.data fort.10
ln -s points.GMT fort.46
ln -s rad$1.data fort.4
ln -s out1$2.data fort.7
ln -s out2$2.data fort.12
ln -s out3$2.data fort.11
ln -s out4$2.data fort.21
ln -s linem.data fort.20
ln -s output$1.data fort.26
ln -s mat$2.data fort.13
ln -s vol$2.data fort.18
ln -s volt$2.data fort.17
ln -s mass$2.data fort.19
ln -s temp$2.data fort.29
ln -s input$1.bc fort.69
ln -s output$2.bc fort.70
ln -s color.map fort.71
ln -s tlist fort.72
ln -s slist fort.75
ln -s slist$2.new fort.77
ln -s tlist.new fort.73
ln -s prof.dat fort.74
ln -s payne-$2.data fort.89
ln -s afract.data fort.90
ln -s dep$2.time fort.88
ln -s bedcond$2.data fort.79
ln -s depl$2.stf fort.80
ln -s depw$2.stf fort.81
ln -s depr$2.stf fort.82
ln -s depd$2.stf fort.83
ln -s deplot$2.data fort.92
ln -s eismint.$2GR fort.93
ln -s eismint.$2HF fort.94
ln -s eismint.$2VF fort.95
ln -s eismint.$2TV fort.96
ln -s times.eismint fort.97
ln -s grip.data fort.98
ln -s console.$2 fort.99
ln -s ume_$2_p_tk fort.50
ln -s ume_$2_p_tp fort.51
ln -s ume_$2_p_uq fort.52
ln -s flist fort.54
ln -s flist$2.new fort.55
ln -s ume_$2_p_vmag fort.65
ln -s ume_$2_p_uvq fort.66
#ln -s ume_$2_p_af fort.54
#ln -s ume_$2_t_vo fort.55
ln -s ume_$2_t_ar fort.56
ln -s ume_$2_t_fr fort.57
ln -s ume_$2_t_tk fort.58
ln -s ume_$2_t_tp fort.59
ln -s ume_$2_d_tp fort.60
ln -s ume_$2_d_uv fort.61
ln -s ume_$2_d_vv fort.62
ln -s ume_$2_d_wv fort.63
ln -s ume_$2_d_af fort.64
time ./map5.x
