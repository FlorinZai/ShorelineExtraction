
d3d_qp('openfile','E:\Delft3D\urd\21\trih-21.dat')
d3d_qp('selectfield','cumulative total transport')
d3d_qp('station','(241,1)..(247,1)  ')
d3d_qp('exportdata','E:\postprocess\extract\21\cumulative total transport.mat')
d3d_qp('closefile')
d3d_qp('openfile','E:\Delft3D\urd\21\trim-21.dat')
d3d_qp('selectfield','water level (when dry: bed level)')
d3d_qp('editt',[ 1:3:1099 ])
d3d_qp('exportdata','E:\postprocess\extract\21\water level (when dry_ bed level).mat')
d3d_qp('selectfield','water level')
d3d_qp('exportdata','E:\postprocess\extract\21\water level.mat')
d3d_qp('selectfield','water depth')
d3d_qp('exportdata','E:\postprocess\extract\21\water depth.mat')
d3d_qp('selectfield','depth averaged velocity')
d3d_qp('component','magnitude')
d3d_qp('exportdata','E:\postprocess\extract\21\depth averaged velocity.mat')
d3d_qp('selectfield','total transport')
d3d_qp('component','magnitude')
d3d_qp('exportdata','E:\postprocess\extract\21\total transport.mat')
d3d_qp('selectfield','bed level in water level points')
d3d_qp('exportdata','E:\postprocess\extract\21\bed level in water level points.mat')
d3d_qp('closefile')

