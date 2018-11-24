#!/usr/bin/env python


def make_xml(flist,nlist,times,size,attrs):
    """Make the xml file for the given time series"""

    topology = """<Topology TopologyType="3DRECTMesh" NumberOfElements="{:d} {:d} {:d}">
    </Topology>""".format(*[i+1 for i in size])

    geometry="""<Geometry GeometryType="VXVYVZ">
            <DataItem Name="xm3" Dimensions="{:d}" NumberType="Float" Precision="8" Format="HDF">
                {name}:/Data/xm1/
            </DataItem>
            <DataItem Name="xm2" Dimensions="{:d}" NumberType="Float" Precision="8" Format="HDF">
                {name}:/Data/xm2/
            </DataItem>
            <DataItem Name="xm1" Dimensions="{:d}" NumberType="Float" Precision="8" Format="HDF">
                {name}:/Data/xm3/
            </DataItem>
          </Geometry>""".format(*[i+1 for i in size[::-1]], name=flist[0])

    grids = ''.join([grid_block(f,n,t,size,attrs) for f,n,t in zip(flist,nlist,times)])

    output="""<?xml version="1.0" ?>
    <!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>
    <Xdmf>
      <Domain>
      {}
      {}
      <Grid Name="Sim" GridType="Collection" CollectionType="Temporal" >
      {}
      </Grid>
      </Domain>
    </Xdmf>""".format(topology, geometry, grids)
    return output


def attribute_block(fname,attr,size,prec=8):
    """Construct an attribute block for the given attribute"""
    output = """<Attribute Name="{a}" AttributeType="Scalar" Center="Cell">
            <DataItem Dimensions="{:d} {:d} {:d}" NumberType="Float" Precision="8" Format="HDF">
              {name}:/Data/{a}/
            </DataItem>
          </Attribute>""".format(*size,name=fname,a=attr)
    return output
def grid_block(fname,n,time,size,fields):
    """Construct a grid block"""

    attrs= ''.join([attribute_block(fname,f,size) for f in fields])

    output="""
    <Grid Name="T{:d}" GridType="Uniform">
          <Topology Reference="/Xdmf/Domain/Topology[1]"/>
          <Geometry Reference="/Xdmf/Domain/Geometry[1]"/>
          <Time Type="Single" Value="{:.5f}" />
          {}
        </Grid>""".format(n,time,attrs,name=fname)
    return output


def parse_range(arg,dtype=int):
    """Parse argument of type (0,10), [0,10]"""

    try:
        # Remove brackets and recombine
        x = ''.join(arg[1:-1])

        # Get rid of the comma
        x = x.split(',')

        # Convert to given data type
        return dtype(x[0]), dtype(x[-1])
    except ValueError:
        # We were given the default arg
        return dtype(arg[0]),dtype(arg[1])




if __name__ == "__main__":
    import argparse
    import h5py
    parser = argparse.ArgumentParser()
    parser.add_argument('-base',type=str,default='out/test',help='Base file name')
    parser.add_argument('-n',type=list,default=[0,10],help='Output range')
    parser.add_argument('-t',type=list,default=[0,10],help='Time range')
    parser.add_argument('-attrs',type=list,default=['Density','Pressure','Vx1','Vx2','Vx3'],help='Attributes to load')

    args = vars(parser.parse_args())


    base = args['base'].split('/')

    dir = ''.join(base[:-1])
    if dir != '':
        dir += '/'
    base = base[-1]


    attrs = args['attrs']

    nlist = args['n']
    nlist = parse_range(args['n'],dtype=int)
    tlist = parse_range(args['t'],dtype=float)
    nlist = list(range(nlist[0], nlist[1]+1))
    num = len(nlist)
    dt = (tlist[1]-tlist[0])/(num-1)
    tlist = [i*dt for i in nlist]


    flist = ['{}_{:d}.h5'.format(base,i) for i in nlist]
    with h5py.File(dir + flist[0],'r') as f:
        dims = f['Data/Density'][...].shape
        dims = tuple([i for i in dims[::-1]])

    print('Compling ', attrs)
    lines = make_xml(flist,nlist,tlist,dims,attrs)
    outfile = dir + '{}.xdmf'.format(base)
    print("Saving to {}".format(outfile))
    with open(outfile,'w') as f:
        f.write(lines)

