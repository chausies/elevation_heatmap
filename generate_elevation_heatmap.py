# Get 'GEBCO_2020.nc' from
# https://www.bodc.ac.uk/data/open_download/gebco/gebco_2020/zip/
# and make sure it's in same directory as this script. 

import matplotlib.pylab as P # Used for matplotlib and numpy
import netCDF4 as nc # Used to read data in .nc (matlab) file


df = 6 # downsampling factor. 6 requires lots of memory. Consider downsampling by 20 or more for speed/stability.
dpi = 1000 # DPI of final image. Higher means higher resolution.

ds = nc.Dataset('GEBCO_2020.nc')
lon, lat, z = ( ds.variables[v] for v in ['lon', 'lat', 'elevation'] )
ind = 2446652862
ea, eb = ind//len(lon), ind%len(lon) # index of everest peak
ao, bo = ea%df, eb%df # offset to start at to include everest
z = P.array(z[ao::df, bo::df])
lon = P.array(lon[ao::df])
lat = P.array(lat[bo::df])
# lon[lon<0], lon[lon>=0] = lon[lon<0]+180, lon[lon>=0]-180
A = lon<0
B = lon>=0
l2 = P.hstack([lon[B]-180, lon[A]+180])
z2 = P.hstack([z[:, B], z[:, A]])

def l_func(p):
  return 0.870700 - 0.131979*p**2 - 0.013791*p**4 + 0.003971*p**10 - 0.001529*p**12
def d_func(p):
  return p*(1.007226 + 0.015085*p**2 - 0.044475*p**6 + 0.028874*p**8 - 0.005916*p**10)

# Natural Earth Projection
def transform(lon, lat):
  rad_lon = lon*P.pi/180
  rad_lat = lat*P.pi/180
  x = l_func(rad_lat)*rad_lon
  y = d_func(rad_lat)
  return (x, y)

fig = P.figure(figsize=(9,8), dpi=dpi)
xx, yy = P.meshgrid(lon, lat)
tx, ty = transform(xx, yy)
ax = fig.add_subplot(211)
cm = ax.pcolormesh(tx, ty, z, shading='gouraud')
ax.set_aspect(1)
P.axis("off")
fig.colorbar(cm, ax=ax)
xx, yy = P.meshgrid(l2, lat)
tx, ty = transform(xx, yy)
ax = fig.add_subplot(212)
cm = ax.pcolormesh(tx, ty, z2, shading='gouraud')
ax.set_aspect(1)
P.axis("off")
fig.colorbar(cm, ax=ax)
P.tight_layout()
print("saving...")
P.savefig(f'elevation_viridis_heatmap_{dpi}dpi.png', dpi=dpi)
# P.show()
