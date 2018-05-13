/* This file is part of the Palabos library.
 *
 * Copyright (C) 2011-2017 FlowKit Sarl
 * Route d'Oron 2
 * 1010 Lausanne, Switzerland
 * E-mail contact: contact@flowkit.com
 *
 * The most recent release of Palabos can be downloaded at
 * <http://www.palabos.org/>
 *
 * The library Palabos is free software: you can redistribute it and/or
 * modify it under the terms of the GNU Affero General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * The library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

/** \file
  * This code solves the steady flow inside an aneurysm. It introduces several
  * new concepts like Guo off lattice boundary conditions, reading of
  * surface geometry STL files, smooth grid refinement and voxelization.
  * Make sure to unpack the file aneurysm.stl.tgz before running the
  * simulation.
  **/

#include "palabos3D.h"
#include "palabos3D.hh"

#include <cstdlib>
#include <cmath>
#include <memory>

using namespace plb;
using namespace std;

typedef double T;
typedef Array<T,3> Vec3D;
#define DESCRIPTOR descriptors::D3Q19Descriptor

#define PADDING 4
static std::string outputDir("./Images/");

// This structure holds all the user-defined parameters, and some
// derived values needed for the simulation.
struct Param
{
    // xml�������
    //----------------------------------
    std::string PCFSF_Geometry_FName;

    T lx, ly, lz;                       // ������ĳߴ磬������λ����
    plint extraLayer;                   // ʹ��Χ�и��󣻽����ڿ��ӻ���Ŀ�ġ����ڷ��棬��extraLayer=0Ҳ�ǿ��Եġ�

    T nu;                               // Kinematic viscosity������λ
    T rho;                              // �����ܶȣ�����λ

    std::vector<T> uVecs;               // ����ٶȣ� ������λ����
    T uLB;                              // ���ӵ�λ�ٶ�(��ֵ����)

    plint flowDirection;                // ����ٶȷ���
    plint referenceResolution;          // �زο����ȵĸ�����

    bool performOutput;
    bool doImages;

    // ����Ĳ���
    //----------------------------------
    plint cxLB, cyLB, czLB;             // PCFSF���ĵ�λ�ã������ӵ�λ����

    T nuLB;                             // Kinematic viscosity�����ӵ�λ
    T rhoLB;                            // �����ܶȣ����ӵ�λ

    T     curUVec;                      // ��ǰ�ٶȣ����ӵ�λ
    plint curUIdx;                      // ��ǰ�ٶȵ�����
    plint resolution;                   // ��ͬϸ�ֲ�ε�ʵ�ʸ�����

    plint nx, ny, nz;                   // ��Χ�е�����ֱ���
    T     omega;                        // �ɳڲ���
    T     dx, dt;                       // ��ɢ�ռ��ʱ�䲽

    plint curIter, maxIter;             // ��������

    plint startLevel;
    plint maxLevel;
    T     epsilon;

    Box3D inlet, outlet;                // �ⲿ����߽磬�����ӵ�λ����
    Box3D lateral1, lateral2, lateral3, lateral4;

    Vec3D uBoundary;                   // ����ٶȱ߽�

    Param() {}

    Param(std::string xmlFname)
    {
        XMLreader document(xmlFname);
        document["geometry"]["mesh"].read(PCFSF_Geometry_FName);
        document["geometry"]["domain"]["x"].read(lx);
        document["geometry"]["domain"]["y"].read(ly);
        document["geometry"]["domain"]["z"].read(lz);

        document["fluid"]["kinematicViscosity"].read(nu);
        document["fluid"]["density"].read(rho);

        document["numerics"]["inletVelocity"].read(uVecs);
        document["numerics"]["flowDirection"].read(flowDirection);
        document["numerics"]["referenceResolution"].read(referenceResolution);
        document["numerics"]["uLB"].read(uLB);

        document["simulation"]["maxIter"].read(maxIter);

        document["simulation"]["maxLevel"].read(maxLevel);
        document["simulation"]["epsilon"].read(epsilon);

        document["simulation"]["performOutput"].read(performOutput);
        document["simulation"]["doImages"].read(doImages);
        plbIOError(flowDirection<0 || flowDirection>2, "Inlet direction of PCFSF must be 0 (x), 1 (y), or 2 (z).");
    }

    void computeLBparameters(plint level)
    {
        // �ֱ��ʲ������1��ÿ�����귽��ķֱ���������1��
        //   ���ݶ��壬����``referenceResolution''������ϸ�ֲ��Ϊ0ʱ�ķֱ���
        resolution = referenceResolution * util::twoToThePower(level);

        curUVec = uVecs[0];
        dx = ly / (resolution - 1.0);
        dt = (uLB/curUVec) * dx;
        nuLB = nu * dt/(dx*dx);
        omega = 1.0/(DESCRIPTOR<T>::invCs2*nuLB+0.5);

        if (flowDirection == 0) {
            nx = util::roundToInt(lx/dx) + 1;
            ny = util::roundToInt(ly/dx) + 1;
            nz = util::roundToInt(lz/dx) + 1;
            uBoundary = Vec3D( uLB, (T)0.0, (T)0.0);
        } else if (flowDirection == 1) {
            nx = util::roundToInt(lx/dx) + 1;
            ny = util::roundToInt(ly/dx) + 1;
            nz = util::roundToInt(lz/dx) + 1;
            uBoundary = Vec3D( (T)0.0, uLB, (T)0.0);
        } else {
            nx = util::roundToInt(lx/dx) + 1;
            ny = util::roundToInt(ly/dx) + 1;
            nz = util::roundToInt(lz/dx) + 1;
            uBoundary = Vec3D( (T)0.0, (T)0.0, uLB);
        }

        cxLB = util::roundToInt(lx/(2*dx));
        cyLB = util::roundToInt(ly/(2*dx));
        czLB = util::roundToInt(lz/(2*dx));

        curIter = 0;
        rhoLB = 1.0;

        defineOuterDomain();
    }

    void defineOuterDomain()
    {
        if (flowDirection == 0) {
            inlet    = Box3D(0,    0,    1,    ny-2, 1, nz-2);
            outlet   = Box3D(nx-1, nx-1, 1,    ny-2, 1, nz-2);
            lateral1 = Box3D(0,    nx-1, 0,    0,    0,    nz-1);
            lateral2 = Box3D(0,    nx-1, ny-1, ny-1, 0,    nz-1);
            lateral3 = Box3D(0,    nx-1, 1,    ny-2, 0,    0);
            lateral4 = Box3D(0,    nx-1, 1,    ny-2, nz-1, nz-1);
        } else if (flowDirection == 1) {
            inlet    = Box3D(1,    nx-2, 0,    0,    1,    nz-2);
            outlet   = Box3D(1,    nx-2, ny-1, ny-1, 1,    nz-2);
            lateral1 = Box3D(0,    0,    0,    ny-1, 0,    nz-1);
            lateral2 = Box3D(nx-1, nx-1, 0,    ny-1, 0,    nz-1);
            lateral3 = Box3D(1,    nx-2, 0,    ny-1, 0,    0);
            lateral4 = Box3D(1,    nx-2, 0,    ny-1, nz-1, nz-1);
        } else {
            inlet    = Box3D(1,    nx-2, 1,    ny-2, 0,       0);
            outlet   = Box3D(1,    nx-2, 1,    ny-2, nz-1, nz-1);
            lateral1 = Box3D(0,    0,    0,    ny-1, 0,    nz-1);
            lateral2 = Box3D(nx-1, nx-1, 0,    ny-1, 0,    nz-1);
            lateral3 = Box3D(1,    nx-2, 0,    0,    0,    nz-1);
            lateral4 = Box3D(1,    nx-2, ny-1, ny-1, 0,    nz-1);
        }
    }

    Box3D boundingBox() const
    {
        return Box3D(0, nx-1, 0, ny-1, 0, nz-1);
    }

    T getInletVelocity(plint& _idx)
    {
        curUIdx = _idx;
        curUVec = uVecs[_idx];
        return curUVec;
    }
};

Param param;

// �ú���������һ�����ܶȣ��Լ�һ������ٶȡ����ڳ�ʼ������Ⱥ��
class VelocityFunction {
public:
    VelocityFunction(Param& _param) : param(_param) {}

    void operator() (plint iX, plint iY, plint iZ, T& density, Vec3D& velocity) const
    {
        density = (T)1.0;

        bool isInlet = false;
        if( 0 == param.flowDirection )
        {
            if( iX == param.inlet.x0 &&
                (iY >= param.inlet.y0 && iY <= param.inlet.y1) &&
                (iZ >= param.inlet.z0 && iZ <= param.inlet.z1) )
                isInlet = true;
        }
        else if( 1 == param.flowDirection )
        {
            if( iY == param.inlet.y0 &&
                (iX >= param.inlet.x0 && iX <= param.inlet.x1) &&
                (iZ >= param.inlet.z0 && iZ <= param.inlet.z1) )
                isInlet = true;
        }
        else
        {
            if( iZ == param.inlet.z0 &&
                (iX >= param.inlet.x0 && iX <= param.inlet.x1) &&
                (iY >= param.inlet.y0 && iY <= param.inlet.y1) )
                isInlet = true;
        }

        if( isInlet )
            velocity = param.uBoundary;
        else
            velocity = Vec3D((T)0.0, (T)0.0, (T)0.0);
    }

private:
    Param param;
};

// ʵ����outer domain�ı߽�����
void outerDomainBoundaries(Param& param, std::auto_ptr<MultiBlockLattice3D<T,DESCRIPTOR> >& lattice, std::auto_ptr< MultiScalarField3D<T> >& rhoBar,
                           std::auto_ptr< MultiTensorField3D<T,3> >& j, MultiScalarField3D<int>& voxelMatrix,
                           std::auto_ptr< OnLatticeBoundaryCondition3D<T,DESCRIPTOR> >& bc)
{
    pcout << "Define boundary conditions." << std::endl;

    lattice->periodicity().toggleAll(false);
    rhoBar->periodicity().toggleAll(false);
    j->periodicity().toggleAll(false);

    bc->setVelocityConditionOnBlockBoundaries(*lattice);
    setBoundaryVelocity(*lattice, lattice->getBoundingBox(), Vec3D((T)0.,(T)0.,(T)0.) );

    Box3D globalDomain(lattice->getBoundingBox());
    if (param.flowDirection == 0) {
        bc->addVelocityBoundary0N(param.inlet, *lattice);
        bc->addPressureBoundary0P(param.outlet, *lattice);
    } else if (param.flowDirection == 1) {
        bc->addVelocityBoundary1N(param.inlet, *lattice);
        bc->addPressureBoundary1P(param.outlet, *lattice);
    } else {
        bc->addVelocityBoundary2N(param.inlet, *lattice);
        bc->addPressureBoundary2P(param.outlet, *lattice);
    }

    setBoundaryVelocity(*lattice, param.inlet, param.uBoundary);
    setBoundaryDensity(*lattice, param.outlet, param.rhoLB);

    // PCFSF���ڲ�������no-dynamics��ʲôҲ������
    defineDynamics(*lattice, voxelMatrix, new NoDynamics<T,DESCRIPTOR>(), voxelFlag::inside);

    initializeAtEquilibrium( *lattice, lattice->getBoundingBox(), VelocityFunction(param) );
}


void computePermeability( std::auto_ptr<MultiBlockLattice3D<T,DESCRIPTOR> >& lattice, Param& param, T deltaP, plb_ofstream& permeablityFile )
{
    pcout << "Computing the permeability." << std::endl;

    plint n = 0;
    T inletVecLB = 0.0;
    T length = 0.0;
    if( param.flowDirection == 0 )
    {
        n = param.nx;
        length = param.lx;
        inletVecLB = param.uBoundary[0];
    }
    else if( param.flowDirection == 1 )
    {
        n = param.ny;
        length = param.ly;
        inletVecLB = param.uBoundary[1];
    }
    else
    {
        n = param.nz;
        length = param.lz;
        inletVecLB = param.uBoundary[2];
    }

    // �������������ƽ���ٶ�
    T scale = param.dx/param.dt;
    T meanVecLB = computeAverage(*computeVelocityComponent(*lattice, lattice->getBoundingBox(), param.flowDirection));
    T meanVec = meanVecLB * scale;

    T delaPLB = deltaP / (param.rho * util::sqr(param.dx/param.dt) * DESCRIPTOR<T>::cs2);
    T pGradLB = delaPLB/(T)(n-1);
    T pGrad = deltaP * param.rho / length;

    T kLB = param.nuLB*meanVecLB / pGradLB;
    T k   = param.nu * meanVec / pGrad;

    if (param.performOutput) {
        pcout << "Average velocity in lattice units  = " << meanVecLB  << std::endl;
        pcout << "Gradient pressure in lattice units = " << pGradLB    << std::endl;
        pcout << "Permeability in lattice units      = " << kLB        << std::endl;
        pcout << "Reynolds number in lattice units   = " << inletVecLB*(n-1) / param.nuLB << std::endl;
        pcout << "Reynolds number in physical units  = " << param.curUVec*length / param.nu    << std::endl;

        permeablityFile << "Average velocity in lattice units   = " << meanVecLB     << std::endl;
        permeablityFile << "Average velocity in physical units  = " << meanVec       << std::endl;
        permeablityFile << "Inlet velocity in lattice units     = " << inletVecLB    << std::endl;
        permeablityFile << "Inlet velocity in physical units    = " << param.curUVec << std::endl;
        permeablityFile << "Gradient pressure in lattice units  = " << pGradLB       << std::endl;
        permeablityFile << "Gradient pressure in physical units = " << pGrad         << std::endl;
        permeablityFile << "Permeability in lattice units       = " << kLB           << std::endl;
        permeablityFile << "Permeability in physical units      = " << k             << std::endl;
        permeablityFile << "Reynolds number in lattice units    = " << inletVecLB*(n-1) / param.nuLB   << std::endl;
        permeablityFile << "Reynolds number in physical units   = " << param.curUVec*length / param.nu << std::endl;
    }
}


// Write VTK file for the flow around the PCFSF, to be viewed with Paraview.
void writeVTK( OffLatticeBoundaryCondition3D<T,DESCRIPTOR,Vec3D>& bc, Param& param, plint level )
{
    VtkImageOutput3D<T> vtkVolumeOut("l"+util::val2str(level)+"_v"+util::val2str(param.curUIdx)+createFileName("_volume", param.curIter, PADDING));

    vtkVolumeOut.writeData<float>( *bc.computeVelocityNorm(param.boundingBox()), "velocityNorm", param.dx/param.dt );
    vtkVolumeOut.writeData<3,float>( *bc.computeVelocity(param.boundingBox()), "velocity", param.dx/param.dt );
    vtkVolumeOut.writeData<float>( *bc.computePressure(param.boundingBox()), "pressure", util::sqr(param.dx/param.dt) );

    VtkImageOutput3D<T> vtkSliceOut("l"+util::val2str(level)+"_v"+util::val2str(param.curUIdx)+createFileName("_slice", param.curIter, PADDING));

    Box3D xSlice(param.cxLB, param.cxLB, 1,          param.ny-2, 1,          param.nz-2);
    Box3D ySlice(1,          param.nx-2, param.cyLB, param.cyLB, 1,          param.nz-2);
    Box3D zSlice(1,          param.nx-2, 1,          param.ny-2, param.czLB, param.czLB);

    vtkSliceOut.writeData<float>( *bc.computeVelocityNorm(xSlice), "velocityNorm", param.dx/param.dt );
    vtkSliceOut.writeData<3,float>( *bc.computeVelocity(ySlice), "velocity", param.dx/param.dt );
    vtkSliceOut.writeData<float>( *bc.computePressure(zSlice), "pressure", util::sqr(param.dx/param.dt)*param.rho );
}


// Write PPM images on slices.
void writeGifs(OffLatticeBoundaryCondition3D<T,DESCRIPTOR,Vec3D>& bc, Param& param, plint level)
{
    Box3D xSlice(param.cxLB, param.cxLB, 1,          param.ny-2, 1,          param.nz-2);
    Box3D ySlice(1,          param.nx-2, param.cyLB, param.cyLB, 1,          param.nz-2);
    Box3D zSlice(1,          param.nx-2, 1,          param.ny-2, param.czLB, param.czLB);

    const plint imSize = 600;
    ImageWriter<T> writer("leeloo");

    writer.writeScaledGif("l"+util::val2str(level)+"_v"+util::val2str(param.curUIdx) + createFileName("_vnorm_xslice", param.curIter, PADDING),
                          *bc.computeVelocityNorm(xSlice), imSize, imSize);
    writer.writeScaledGif("l"+util::val2str(level)+"_v"+util::val2str(param.curUIdx) + createFileName("_vnorm_yslice", param.curIter, PADDING),
                          *bc.computeVelocityNorm(ySlice), imSize, imSize);
    writer.writeScaledGif("l"+util::val2str(level)+"_v"+util::val2str(param.curUIdx) + createFileName("_vnorm_zslice", param.curIter, PADDING),
                          *bc.computeVelocityNorm(zSlice), imSize, imSize);
}


// �ú���׼����ִ��ʵ�ʵķ���
std::auto_ptr<MultiBlockLattice3D<T,DESCRIPTOR> > run ( Param& param, plint level, plb_ofstream& permeablityFile, MultiBlockLattice3D<T,DESCRIPTOR>* iniVal=0 )
{
    // ���ȼ�����Ӳ���
    param.computeLBparameters(level);

    if (param.performOutput) {
        pcout << "tau = " << 1.0/param.omega << std::endl;
        pcout << "dx  = " <<    param.dx     << std::endl;
        pcout << "dt  = " <<    param.dt     << std::endl;
        pcout << "Lattice viscosity nu    = " << param.nuLB         << std::endl;
        pcout << "Number of iterations in an integral time scale: " << (plint) (1.0/param.dt) << std::endl;
    }

    /*
     * ��ȡPCFSF�ļ���
     */

    // ��PCFSF�ı��漸�Σ���ASCII�������STL�ļ���ʾ����������ŵ��������μ�����ɵ����ݽṹ��
    //   ����DBL��ʾ����˫���ȣ�һ���Ƽ�����ѡ��
    pcout << std::endl << "Reading STL data for the PCFSF geometry." << std::endl;
    Vec3D centerLB(param.cxLB, param.cyLB, param.czLB);
    TriangleSet<T> triangleSet(param.PCFSF_Geometry_FName, DBL);

    // ��PCFSF�ŵ��������еĺ���λ�á��˴������Ȼ��PCFSF�İ�Χ�У�Ȼ���ֶ������伸������
    // �ھ���PCFSF���ε�STL�ļ��������ĵ㣬��(0, 0, 0)������£����뽫����"PCFSFCenter"
    // �ֶ����õ�(0, 0, 0)
    Cuboid<T> bCuboid = triangleSet.getBoundingCuboid();
    Vec3D PCFSFCenter = (T) 0.5 * (bCuboid.lowerLeftCorner + bCuboid.upperRightCorner);
    triangleSet.translate(-PCFSFCenter);
    triangleSet.scale(1.0/param.dx); // ���ڿ�ʼ���ø��ӵ�λ...
    triangleSet.translate(centerLB);
    triangleSet.writeBinarySTL(outputDir+"PCFSF_LB.stl");

    // ���¼����ǵ��ʹ��롣�佫�û�������PCFSF���漸��ת��Ϊpalabos�ڲ�ʹ�õġ�
    //   ��Ϊ��Ч�����ݽṹ��������ʹ��TriangleBoundary3D�ṹ��ָ�����ʵı߽�����
    pcout << std::endl << "Constructing TriangleBoundary3D at level " << level << std::endl;
    plint margin = 1;      // �ϰ�����Χ�����cells�Ķ���ԣ��
    plint borderWidth = 1; // ��ΪGuo�߽������ڵ�cell�������á�Ҫ��margin>=borderWidth
    plint blockSize = 0;   // �����0���ʾ����ʹ��ϡ���ʾ

    DEFscaledMesh<T> defMesh(triangleSet, 0, param.flowDirection, margin, Dot3D(0, 0, 0));
    TriangleBoundary3D<T> boundary(defMesh);
    //boundary.getMesh().inflate();

    /*
     * ���ػ�����
     */

    // PCSFS��������ⲿ(���ڲ��෴)�������⡣Ϊ�ˣ�����ʶ�����ڼ�����֮�ڵ�lattice�ڵ�
    //   ���������ⲿ�Ľڵ�������֡���ͨ����������ػ����̴���
    pcout << std::endl << "Running VoxelizedDomain3D at level " << level << std::endl;
    plint extendedEnvelopeWidth = 2;   // ���Ƶ�off-lattice�߽���������ΪGuo��off lattice�߽�������Ҫ2-cells�ڽ�����
    const int flowType = voxelFlag::outside;
    VoxelizedDomain3D<T> voxelizedDomain ( boundary, flowType, param.boundingBox(), borderWidth, extendedEnvelopeWidth, blockSize );

    if (param.performOutput) {
        pcout << getMultiBlockInfo(voxelizedDomain.getVoxelMatrix()) << std::endl;
    }

    /*
     * ����lattice���ܶȺͶ���blocks.
     */

    pcout << "Generating the lattice, the rhoBar and j fields." << std::endl;
    plint envelopeWidth = 1;  // ���ڱ�׼��BGK dynamics��
    std::auto_ptr<MultiBlockLattice3D<T,DESCRIPTOR> > lattice = generateMultiBlockLattice<T,DESCRIPTOR> ( voxelizedDomain.getVoxelMatrix(),
                                                                envelopeWidth, new IncBGKdynamics<T,DESCRIPTOR>(param.omega) ); // ��ģ�����ٶȵ��ڶ���
    lattice->toggleInternalStatistics(false);

    // �ܶȳ�rhoBar�Ͷ�����j����ײ��ʵ��outflow�߽�����ʱ��Ҫʹ��
    std::auto_ptr< MultiScalarField3D<T> > rhoBar( generateMultiScalarField<T>((MultiBlock3D&) *lattice, envelopeWidth).release() );
    rhoBar->toggleInternalStatistics(false);

    std::auto_ptr< MultiTensorField3D<T,3> > j( generateMultiTensorField<T,3>((MultiBlock3D&) *lattice, envelopeWidth).release() );
    j->toggleInternalStatistics(false);

    std::vector<MultiBlock3D*> lattice_rho_bar_j_arg;
    lattice_rho_bar_j_arg.push_back(lattice.get());
    lattice_rho_bar_j_arg.push_back(rhoBar.get());
    lattice_rho_bar_j_arg.push_back(j.get());
    integrateProcessingFunctional( new ExternalRhoJcollideAndStream3D<T,DESCRIPTOR>(), lattice->getBoundingBox(), lattice_rho_bar_j_arg, 0);
    integrateProcessingFunctional( new BoxRhoBarJfunctional3D<T,DESCRIPTOR>(), lattice->getBoundingBox(), lattice_rho_bar_j_arg, 2);
    // �ܶȳ�rhoBar�Ͷ�����j���ڲ��3���㣬��Ϊ�߽������ڲ��0����

    /*
     * ��PCFSF������off-lattice�߽�������outer-domain�߽�����
     */
    pcout << "Generating boundary conditions." << std::endl;
    std::auto_ptr< OffLatticeBoundaryCondition3D<T,DESCRIPTOR,Vec3D> > boundaryCondition;

    BoundaryProfiles3D<T,Vec3D> profiles;
    bool useAllDirections=true;
    OffLatticeModel3D<T,Vec3D>* offLatticeModel=0;

    bool velIsJ = false;
    profiles.setWallProfile(new NoSlipProfile3D<T>);
    offLatticeModel = new GuoOffLatticeModel3D<T,DESCRIPTOR> (
            new TriangleFlowShape3D<T,Array<T,3> >(voxelizedDomain.getBoundary(), profiles),
            flowType, useAllDirections );
    offLatticeModel->setVelIsJ(velIsJ);
    boundaryCondition = std::auto_ptr< OffLatticeBoundaryCondition3D<T,DESCRIPTOR,Vec3D> >( new OffLatticeBoundaryCondition3D<T,DESCRIPTOR,Vec3D>(offLatticeModel, voxelizedDomain, *lattice) );

    boundaryCondition->insert();

    /*
     * ����PCFSF�ⲿ������ı߽������ͳ�ʼ����
     */
    // �߽������㷨��outer domain.
    for( plint idx=0; idx < plint(param.uVecs.size()); idx++ )
    {
        permeablityFile << std::endl << "Running simulation for the inlet velocity: " << param.getInletVelocity(idx) << std::endl;
        permeablityFile << "===========================================================" << std::endl;

        std::auto_ptr< OnLatticeBoundaryCondition3D<T,DESCRIPTOR> > outerBoundaryCondition( createLocalBoundaryCondition3D<T,DESCRIPTOR>() );
        outerDomainBoundaries(param, lattice, rhoBar, j, voxelizedDomain.getVoxelMatrix(), outerBoundaryCondition);

        // �����ܶȳ�rhoBar�Ͷ�����j
        applyProcessingFunctional(new BoxRhoBarJfunctional3D<T,DESCRIPTOR>(), lattice->getBoundingBox(), lattice_rho_bar_j_arg);
        if(iniVal) {
            Box3D toDomain(lattice->getBoundingBox());
            Box3D fromDomain(toDomain.shift(margin,margin,margin)); // ����������ʱ��margin�ڳߴ緽�淭�����˴�ͨ��һ��shift��ȡ����һЧ��
            copy(*iniVal, fromDomain, *lattice, toDomain, modif::staticVariables);
            computeRhoBarJ(*lattice, *rhoBar, *j, lattice->getBoundingBox());
            boundaryCondition->apply(lattice_rho_bar_j_arg);
        }

        // ��һ��ѡ�����(ѹ��)����ʱ����Ҫ���ValueTracer���Ա�ó������ض�������ϸ�ֲ�Σ��Ѿ��ﵽ�ȶ�̬�Ľ��ۣ���ֹͣ����
        plint convergeIter=20;
        util::ValueTracer<T> energeTracer(0.05*convergeIter, param.resolution, param.epsilon);//0.05*convergenceIter
        global::timer("iteration").restart();
        lattice->resetTime(param.curIter*param.dt);
        bool checkForErrors = true;

        T inletPressure = param.rhoLB;
        T outletPressure = param.rhoLB;
        T pressureDrop = inletPressure - outletPressure;
        T pressureScale = param.rho * util::sqr(param.dx/param.dt) * DESCRIPTOR<T>::cs2;

        plb_ofstream iterFile((outputDir+"IterationInfo_p"+util::val2str(idx)+".txt").c_str());
        iterFile << "The iteration information for the inlet velocity: " << param.curUVec << std::endl;
        iterFile << "==========================================================================" << std::endl;
        iterFile << "Time\t\t" << "Iteration\t\t" << "Pressure drop" << std::endl;

        // ��ײ��������.
        while( !energeTracer.hasConverged() && param.curIter < param.maxIter )//
        {
            if (param.curIter % 50==0 && param.performOutput) {
                inletPressure = pressureScale * computeAverageDensity(*lattice, param.inlet);
                outletPressure = pressureScale * computeAverageDensity(*lattice, param.outlet);
                pressureDrop = inletPressure - outletPressure;

                pcout << "T= " << param.curIter*param.dt << "; " << "Iteration " << param.curIter
                << "; Inlet pressure: " << inletPressure << "; " << "Pressure drop: "  << pressureDrop << std::endl;

                iterFile << "Time " << param.curIter*param.dt << "\t\t" << param.curIter << "\t\t" << pressureDrop << std::endl;

                pcout << "Writing VTK ..." << endl;
                writeVTK(*boundaryCondition, param, level);

                pcout << "Writing Gif ..." << endl;
                writeGifs(*boundaryCondition, param, level);

                pcout << "Average Energy: " << boundaryCondition->computeAverageEnergy()*util::sqr(param.dx/param.dt);
            }

            if( param.curIter % convergeIter==0 )
                energeTracer.takeValue(computeAverageEnergy(*lattice),true);//getStoredAverageEnergy(*lattice)

            lattice->executeInternalProcessors();
            lattice->incrementTime();

            if (checkForErrors) {
                abortIfErrorsOccurred();
                checkForErrors = false;
            }

            ++param.curIter;
        }

        inletPressure = pressureScale * computeAverageDensity(*lattice, param.inlet);
        outletPressure = pressureScale * computeAverageDensity(*lattice, param.outlet);
        pressureDrop = inletPressure - outletPressure;

        computePermeability(lattice, param, pressureDrop, permeablityFile);

        pcout << "Elapsed time: " << global::timer("iteration").stop() << std::endl;
        pcout << "Total elapsed time: " << global::timer("global").getTime() << std::endl;

        pcout << "End of simulation at iteration " << param.curIter << std::endl;
        permeablityFile << "End of simulation at iteration " << param.curIter << std::endl;

        permeablityFile << "Elapsed time: " << global::timer("iteration").stop() << std::endl;
        permeablityFile << "Total elapsed time: " << global::timer("global").getTime() << std::endl;

        param.curIter = 0;

        iterFile.close();
    }

    return lattice;
}


int main(int argc, char* argv[])
{
    //std::cout.precision(10);
    //std::scientific(std::cout);

    plbInit(&argc, &argv);
    global::directories().setOutputDir("./Images/");
    global::IOpolicy().activateParallelIO(false);

    string paramXmlFileName;
    try {
        global::argv(1).read(paramXmlFileName);
    }
    catch (PlbIOException& exception) {
        pcout << "Wrong parameters; the syntax is: " << (std::string)global::argv(0) << " parameter-input-file.xml" << std::endl;
        return -1;
    }

    plb_ofstream permeablityFile((outputDir+"Permeability.txt").c_str());

    // ��ȡ����XML�����ļ�(����Ҳ�����˴�����ע��).
    try {
        param = Param(paramXmlFileName);
    }
    catch (PlbIOException& exception) {
        pcout << "Error in input file " << paramXmlFileName << ": " << exception.what() << std::endl;
        return -1;
    }

    global::timer("global").start();
    plint iniLevel=0;
    std::auto_ptr< MultiBlockLattice3D<T,DESCRIPTOR> > iniConditionLattice(0);
    // ���´�������ƽ������ϸ��ֱ�������ĸ��ƽ����ʾϸ�ֲ������1������������ÿ�������Ϸ�����
    //   ������ϸ����dx��dt��Ҫ�ı䡣dt�ǰ�dx^2 (��ɢ��Ϊ)�仯
    try {
        for (plint level=iniLevel; level<=param.maxLevel; ++level) {
            pcout << std::endl << "Running new simulation at level " << level << std::endl;
            permeablityFile << std::endl << "Running new simulation at level " << level << std::endl;
            std::auto_ptr< MultiBlockLattice3D<T,DESCRIPTOR> > convergedLattice ( run(param, level, permeablityFile, iniConditionLattice.get()) );
            if (level != param.maxLevel) {
                plint dxScale = -1;
                plint dtScale = -2;

                // ǰһ�������������ķ�������һ�������Σ��ڽ����ʵ��Ĳ�ֵ֮�󣩷���ĳ�ʼ����
                iniConditionLattice = std::auto_ptr<MultiBlockLattice3D<T,DESCRIPTOR> > (refine(*convergedLattice, dxScale, dtScale, new IncBGKdynamics<T,DESCRIPTOR>(1.)) );
            }
        }

        permeablityFile.close();
    }
    catch(PlbException& exception) {
        pcout << exception.what() << std::endl;
        return -1;
    }
}
