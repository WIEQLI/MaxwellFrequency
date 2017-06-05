export getEdgeIntegralOfMagneticDipole

"""
        function MaxwellFrequency.Utils.getEdgeIntegralOfMagneticDipole

        s = getEdgeIntegralOfMagneticDipole(mesh,dipoleLoc,dipoleMoment)
        s = getEdgeIntegralOfMagneticDipole(mesh,dipoleLoc,dipoleMoment,normalize)

        Compute the integral of a piecewise linear edge grid basis projection
        for the modeling of a magnetic dipole. This function can be used to evaluate the
        source and receiver term of an infinitly small closed loop in electromagnetic
        modelling.

        This function is built on getEdgeIntegralOfPolygonalChain. The dipole 
        is first located in one cell of the mesh, then the 12 edges of the cell
        are used to model 6 orthogonal loops on P faces, with 2 loops for each 
        direction. The x-component moment is distributed to the two x-facing 
        loops according to the relative x-location of the dipole within the cell.

        INPUT
        mesh: tensor mesh or octree mesh
        dipoleLoc: a (x,y,z) vector of the location of the dipole
        dipoleMoment: a (x,y,z) vector of the vector moment of the dipole
        normalize: if true, model a unit dipole by regardless of the norm of dipolo moment (boolean, optional)

        OUTPUT
        s: source vector

        NOTE (Dikun)
        Currently only work with octree mesh, because tensor mesh dones't have 'getNodeToCellCenteredMatrix'.

"""
function getEdgeIntegralOfMagneticDipole(mesh::OcTreeMesh, dipoleLoc::Array{Float64,1},dipoleMoment::Array{Float64,1}; normalize=true)

    if normalize
        moment = dipoleMoment / norm(dipoleMoment)
    else
        moment = dipoleMoment
    end

    cc = getCellCenteredGrid(mesh)
    distance = (cc[:,1]-dipoleLoc[1]).^2 + (cc[:,2]-dipoleLoc[2]).^2 + (cc[:,3]-dipoleLoc[3]).^2 
    ccInd = sortperm(distance)[1] # cell center index of the cube that encloses the loc
    nodeInd = find(getNodeToCellCenteredMatrix(mesh)[ccInd,:]) # find the neighboring nodes's indices
    nodes = getNodalGrid(mesh)
    
    
    (xmin, xmax, ymin, ymax, zmin, zmax) = ( nodes[nodeInd[1],1], nodes[nodeInd[2],1],
                                             nodes[nodeInd[1],2], nodes[nodeInd[3],2],
                                             nodes[nodeInd[1],3], nodes[nodeInd[5],3]) # cell bounds
    (dx, dy, dz) = (xmax-xmin, ymax-ymin, zmax-zmin) # cell size
    
    # build x-dipole
    (d1, d2) = ( dipoleLoc[1]-xmin, xmax-dipoleLoc[1] )
    xdipole1_poly = nodes[nodeInd[[1,3,7,5,1]],:]
    xdipole1_amp = moment[1] / (dy*dz)  *  d2/(d1+d2)
    xdipole2_poly = nodes[nodeInd[[2,4,8,6,2]],:]
    xdipole2_amp = moment[1] / (dy*dz)  *  d1/(d1+d2)
    xdipole1 = complex(getEdgeIntegralOfPolygonalChain(mesh,xdipole1_poly,normalize=false)) * xdipole1_amp
    xdipole2 = complex(getEdgeIntegralOfPolygonalChain(mesh,xdipole2_poly,normalize=false)) * xdipole2_amp
    
    # build y-dipole
    (d1, d2) = ( dipoleLoc[2]-ymin, ymax-dipoleLoc[2] )
    ydipole1_poly = nodes[nodeInd[[1,5,6,2,1]],:]
    ydipole1_amp = moment[2] / (dx*dz)  *  d2/(d1+d2)
    ydipole2_poly = nodes[nodeInd[[3,7,8,4,3]],:]
    ydipole2_amp = moment[2] / (dx*dz)  *  d1/(d1+d2)
    ydipole1 = complex(getEdgeIntegralOfPolygonalChain(mesh,ydipole1_poly,normalize=false)) * ydipole1_amp
    ydipole2 = complex(getEdgeIntegralOfPolygonalChain(mesh,ydipole2_poly,normalize=false)) * ydipole2_amp
    
    # build z-dipole
    (d1, d2) = ( dipoleLoc[3]-zmin, zmax-dipoleLoc[3] )
    zdipole1_poly = nodes[nodeInd[[1,2,4,3,1]],:]
    zdipole1_amp = moment[3] / (dx*dy)  *  d2/(d1+d2)
    zdipole2_poly = nodes[nodeInd[[5,6,8,7,5]],:]
    zdipole2_amp = moment[3] / (dx*dy)  *  d1/(d1+d2)
    zdipole1 = complex(getEdgeIntegralOfPolygonalChain(mesh,zdipole1_poly,normalize=false)) * zdipole1_amp
    zdipole2 = complex(getEdgeIntegralOfPolygonalChain(mesh,zdipole2_poly,normalize=false)) * zdipole2_amp
    
    return xdipole1 + xdipole2 + ydipole1 + ydipole2 + zdipole1 + zdipole2

end

