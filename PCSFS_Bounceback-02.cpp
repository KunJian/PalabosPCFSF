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

using namespace plb;
using namespace std;

typedef double T;
typedef Array<T,3> Vec3D;
#define DESCRIPTOR descriptors::D3Q19Descriptor

#define PADDING 5
static std::string outputDir("./Images/");

// This structure holds all the user-defined parameters, and some
// derived values needed for the simulation.
struct Param
{
    // xml�������
    //----------------------------------
    std::string PCFSF_Geometry_FName;

    T lx, ly, lz;                       // ������ĳߴ磬����λ
    plint extraLayer;                   // ʹ��Χ�и��󣻽����ڿ��ӻ���Ŀ�ġ����ڷ��棬��extraLayer=0Ҳ�ǿ��Եġ�

    T nu;                               // Kinematic viscosity������λ
    T rho;                              // �����ܶȣ�����λ

    std::vector<T> deltaP;              // ѹ�������ӵ�λ
    plint flowDirection;                // ����ٶȷ���
    plint referenceResolution;          // �زο����ȵĸ�����

    plint curIter, maxIter;             // ��ǰ������������

    plint maxLevel;
    T     epsilon;

    bool performOutput;
    bool doImages;

    // ����Ĳ���
    //----------------------------------
    plint cxLB, cyLB, czLB;             // PCFSF���ĵ�λ�ã����ӵ�λ
    plint nx, ny, nz;                   // ��Χ�е�����ֱ���
    Box3D inlet, outlet;                // �ⲿ����߽磬���ӵ�λ
    Box3D lateral1, lateral2, lateral3, lateral4;
    T     nuLB;                         // Kinematic viscosity�����ӵ�λ

    T     curDeltaP;                    // ѹ�������ӵ�λ
    plint curPIdx;                      // ��ǰѹ��������
    plint resolution;                   // ��ͬϸ�ֲ�ε�ʵ�ʸ�����
    T     omega;                        // �ɳڲ���
    T     dx, dt;                       // ��ɢ�ռ��ʱ�䲽

    plint startLevel;

    Param() {}

    Param(std::string xmlFname, plb_ofstream& permeablityFile)
    {
        XMLreader document(xmlFname);
        document["geometry"]["mesh"].read(PCFSF_Geometry_FName);
        document["geometry"]["domain"]["x"].read(lx);
        document["geometry"]["domain"]["y"].read(ly);
        document["geometry"]["domain"]["z"].read(lz);
        document["geometry"]["extraLayer"].read(extraLayer);

        document["fluid"]["kinematicViscosity"].read(nu);
        document["fluid"]["density"].read(rho);

        document["numerics"]["pressureDrop"].read(deltaP);
        document["numerics"]["flowDirection"].read(flowDirection);
        document["numerics"]["referenceResolution"].read(referenceResolution);

        document["simulation"]["maxIter"].read(maxIter);

        document["simulation"]["maxLevel"].read(maxLevel);
        document["simulation"]["epsilon"].read(epsilon);

        document["simulation"]["performOutput"].read(performOutput);
        document["simulation"]["doImages"].read(doImages);

        plbIOError(flowDirection<0 || flowDirection>2, "Inlet direction of PCFSF must be 0 (x), 1 (y), or 2 (z).");
    }

    void computeLBparameters(plint level, plb_ofstream& permeablityFile)
    {
        omega = 1.0;
        nuLB  = ((T)1/omega- (T)0.5)/DESCRIPTOR<T>::invCs2;

        // �ֱ��ʲ������1��ÿ�����귽��ķֱ���������1��
        //   ���ݶ��壬����``referenceResolution''������ϸ�ֲ��Ϊ0ʱ�ķֱ���
        resolution = referenceResolution * util::twoToThePower(level);

        dx = ly / (resolution - 1.0);
        dt = nuLB/nu * (dx*dx);

        nx = util::roundToInt(lx/dx) + 1;
        ny = util::roundToInt(ly/dx) + 1;
        nz = util::roundToInt(lz/dx) + 1;

        cxLB = util::roundToInt(lx/(2*dx));
        cyLB = util::roundToInt(ly/(2*dx));
        czLB = util::roundToInt(lz/(2*dx));

        curIter = 0;

        defineOuterDomain();

        if (performOutput) {
            pcout << "tau = " << 1.0/omega << std::endl;
            pcout << "dx  = " <<    dx     << std::endl;
            pcout << "dt  = " <<    dt     << std::endl;
            pcout << "Physical viscosity nu   = " << nu      << std::endl;
            pcout << "Lattice viscosity nuLB  = " << nuLB    << std::endl;
            pcout << "Max number of iterations: " << maxIter << std::endl;
            pcout << "Number of iterations in an integral time scale: " << (plint) (1.0/dt) << std::endl;

            permeablityFile << "Max number of iterations: " << maxIter << std::endl;
            permeablityFile << "Number of iterations in an integral time scale: " << (plint) (1.0/dt) << std::endl;
            permeablityFile << "tau = " << 1.0/omega << std::endl;
            permeablityFile << "dx  = " <<    dx     << std::endl;
            permeablityFile << "dt  = " <<    dt     << std::endl;
            permeablityFile << "Physical viscosity nu   = " << nu   << std::endl;
            permeablityFile << "Lattice viscosity nuLB  = " << nuLB << std::endl;
        }
    }

    void defineOuterDomain()
    {
        if (flowDirection == 0) {
            inlet    = Box3D(0,    0,    1,    ny-2, 1,    nz-2);
            outlet   = Box3D(nx-1, nx-1, 1,    ny-2, 1,    nz-2);
            lateral1 = Box3D(0,    nx-1, 0,    0,    0,    nz-1);
            lateral2 = Box3D(0,    nx-1, ny-1, ny-1, 0,    nz-1);
            lateral3 = Box3D(0,    nx-1, 1,    ny-2, 0,    0);
            lateral4 = Box3D(0,    nx-1, 1,    ny-2, nz-1, nz-1);
        } else if (flowDirection == 1) {
            inlet    = Box3D(1,    nx-2, 0,    0,    1,    nz-2);
            outlet   = Box3D(1,    nx-2, ny-1, ny-1, 1,    nz-2);
            lateral1 = Box3D(0,    nx-1, 0,    ny-1, 0,    0);
            lateral2 = Box3D(0,    nx-1, 0,    ny-1, nz-1, nz-1);
            lateral3 = Box3D(0,    0,    0,    ny-1, 1,    nz-2);
            lateral4 = Box3D(nx-1, nx-1, 0,    ny-1, 1,    nz-2);
        } else {
            inlet    = Box3D(1,    nx-2, 1,    ny-2, 0,    0);
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

    T getPressureDrop(plint& _idx)
    {
        curPIdx = _idx;
        curDeltaP = deltaP[_idx];
        return curDeltaP;
    }
};

Param param;

// �ú���������һ�����ٶȣ��Լ�һ�������������������Լ�С��ѹ�������ڳ�ʼ������Ⱥ��
class PressureGradient {
public:
    PressureGradient(Param& _param) : param(_param)
    { }
    void operator() (plint iX, plint iY, plint iZ, T& density, Vec3D& velocity) const
    {
        velocity.resetToZero();

        if( 0 == param.flowDirection )
            density = (T)1 - param.curDeltaP*DESCRIPTOR<T>::invCs2 / (T)(param.nx-1) * (T)iX;
        else if( 1 == param.flowDirection )
            density = (T)1 - param.curDeltaP*DESCRIPTOR<T>::invCs2 / (T)(param.ny-1) * (T)iY;
        else
            density = (T)1 - param.curDeltaP*DESCRIPTOR<T>::invCs2 / (T)(param.nz-1) * (T)iZ;
    }

private:
    Param param;
};


// ʵ����outer domain�ı߽�����
void outerDomainBoundaries(Param& param, std::auto_ptr<MultiBlockLattice3D<T,DESCRIPTOR> >& lattice, MultiScalarField3D<int>& voxelMatrix,
                           std::auto_ptr< OnLatticeBoundaryCondition3D<T,DESCRIPTOR> >& bc)
{
    /*
     * ���ó�ʼ����
     */
    // ��ʼ���������еط����ǳ�ѹ����velocity-at-infinity
    pcout << "Define boundary conditions." << std::endl;

    lattice->periodicity().toggleAll(false);

    if( param.flowDirection == 0 )
    {
        bc->addPressureBoundary0N(param.inlet, *lattice);
        setBoundaryDensity(*lattice, param.inlet, (T) 1.);

        bc->addPressureBoundary0P(param.outlet, *lattice);
        setBoundaryDensity(*lattice, param.outlet, (T) 1. - param.curDeltaP*DESCRIPTOR<T>::invCs2);
    }
    else if( param.flowDirection == 1 )
    {
        bc->addPressureBoundary1N(param.inlet, *lattice);
        setBoundaryDensity(*lattice, param.inlet, (T) 1.);

        bc->addPressureBoundary1P(param.outlet, *lattice);
        setBoundaryDensity(*lattice, param.outlet, (T) 1. - param.curDeltaP*DESCRIPTOR<T>::invCs2);
    }
    else
    {
        bc->addPressureBoundary2N(param.inlet, *lattice);
        setBoundaryDensity(*lattice, param.inlet, (T) 1.);

        bc->addPressureBoundary2P(param.outlet, *lattice);
        setBoundaryDensity(*lattice, param.outlet, (T) 1. - param.curDeltaP*DESCRIPTOR<T>::invCs2);
    }

    // PCFSF����߽硢�Լ���������ĸ����棬����bounce-back�߽�����
    defineDynamics(*lattice, voxelMatrix, new BounceBack<T,DESCRIPTOR>(), voxelFlag::outerBorder);
    defineDynamics(*lattice, param.lateral1, new BounceBack<T,DESCRIPTOR>((T)1.));
    defineDynamics(*lattice, param.lateral2, new BounceBack<T,DESCRIPTOR>((T)1.));
    defineDynamics(*lattice, param.lateral3, new BounceBack<T,DESCRIPTOR>((T)1.));
    defineDynamics(*lattice, param.lateral4, new BounceBack<T,DESCRIPTOR>((T)1.));

    // PCFSF���ڲ�������no-dynamics��ʲôҲ������
    defineDynamics(*lattice, voxelMatrix, new NoDynamics<T,DESCRIPTOR>(), voxelFlag::inside);

    pcout << "Initilization of rho and u." << std::endl;
    initializeAtEquilibrium( *lattice, lattice->getBoundingBox(), PressureGradient(param) );

    lattice->initialize();
}


void computePermeability(std::auto_ptr<MultiBlockLattice3D<T,DESCRIPTOR> >& lattice, Param& param, plb_ofstream& permeablityFile)
{
    pcout << "Computing the permeability." << std::endl;

    // �������������ƽ���ٶ�
    T meanVecLB = computeAverage(*computeVelocityComponent(*lattice, lattice->getBoundingBox(), param.flowDirection));
    T inletVecLB = computeAverage(*computeVelocityComponent(*lattice, param.inlet, param.flowDirection));
    T scale = param.dx/param.dt;

    T meanVec = meanVecLB * scale;
    T inletVec =  inletVecLB * scale;

    plint n = 0;
    T length = 0.0;
    if( param.flowDirection == 0 )
    {
        n = param.nx;
        length = param.lx;
    }
    else if( param.flowDirection == 1 )
    {
        n = param.ny;
        length = param.ly;
    }
    else
    {
        n = param.nz;
        length = param.lz;
    }

    T pGradLB = param.curDeltaP/(T)(n-1);
    T pGrad = param.curDeltaP * scale*scale * param.rho / length;

    T kLB = param.nuLB*meanVecLB / pGradLB;
    T k   = param.nu * inletVec / pGrad;

    if (param.performOutput) {
        pcout << "Average velocity in lattice units  = " << meanVecLB  << std::endl;
        pcout << "Gradient pressure in lattice units = " << pGradLB    << std::endl;
        pcout << "Permeability in lattice units      = " << kLB        << std::endl;
        pcout << "Reynolds number in lattice units   = " << inletVecLB*(n-1) / param.nuLB << std::endl;
        pcout << "Reynolds number in physical units  = " << inletVec*length / param.nu    << std::endl;

        permeablityFile << "Average velocity in lattice units   = " << meanVecLB   << std::endl;
        permeablityFile << "Average velocity in physical units  = " << meanVec     << std::endl;
        permeablityFile << "Inlet velocity in lattice units     = " << inletVecLB  << std::endl;
        permeablityFile << "Inlet velocity in physical units    = " << inletVec    << std::endl;
        permeablityFile << "Gradient pressure in lattice units  = " << pGradLB     << std::endl;
        permeablityFile << "Gradient pressure in physical units = " << pGrad       << std::endl;
        permeablityFile << "Permeability in lattice units       = " << kLB         << std::endl;
        permeablityFile << "Permeability in physical units      = " << k           << std::endl;
        permeablityFile << "Reynolds number in lattice units    = " << inletVecLB*(n-1) / param.nuLB << std::endl;
        permeablityFile << "Reynolds number in physical units   = " << inletVec*length / param.nu    << std::endl;
    }
}


// ����Χ��PCFSF��������VTK�ļ�������Paraview�鿴
void writeVTK(MultiBlockLattice3D<T,DESCRIPTOR>& lattice, Param& param, plint level)
{
    VtkImageOutput3D<T> vtkVolumeOut("l"+util::val2str(level)+"_p"+util::val2str(param.curPIdx)+createFileName("_volume", param.curIter, PADDING));

    vtkVolumeOut.writeData<float>( *computeVelocityNorm(lattice, lattice.getBoundingBox()), "velocityNorm", param.dx/param.dt );
    vtkVolumeOut.writeData<3,float>(*computeVelocity(lattice, lattice.getBoundingBox()), "velocity", param.dx/param.dt );
    vtkVolumeOut.writeData<float>( *computeDensity(lattice, lattice.getBoundingBox()), "pressure", util::sqr(param.dx/param.dt)*param.rho );

    VtkImageOutput3D<T> vtkSliceOut("l"+util::val2str(level)+"_p"+util::val2str(param.curPIdx)+createFileName("_slice", param.curIter, PADDING));

    Box3D xSlice(param.cxLB, param.cxLB, 1,          param.ny-2, 1,          param.nz-2);
    Box3D ySlice(1,          param.nx-2, param.cyLB, param.cyLB, 1,          param.nz-2);
    Box3D zSlice(1,          param.nx-2, 1,          param.ny-2, param.czLB, param.czLB);

    vtkSliceOut.writeData<float>( *computeVelocityNorm(lattice, xSlice), "velocityNorm", param.dx/param.dt );
    vtkSliceOut.writeData<3,float>(*computeVelocity(lattice, ySlice), "velocity", param.dx/param.dt );
    vtkSliceOut.writeData<float>( *computeDensity(lattice, zSlice), "pressure", util::sqr(param.dx/param.dt)*param.rho );
}


// ����������Ƭ��Gifͼ
void writeGifs(MultiBlockLattice3D<T,DESCRIPTOR>& lattice, Param& param, plint level)
{
    Box3D xSlice(param.cxLB, param.cxLB, 1,          param.ny-2, 1,          param.nz-2);
    Box3D ySlice(1,          param.nx-2, param.cyLB, param.cyLB, 1,          param.nz-2);
    Box3D zSlice(1,          param.nx-2, 1,          param.ny-2, param.czLB, param.czLB);

    const plint imSize = 600;
    ImageWriter<T> writer("leeloo");

    writer.writeScaledGif("l"+util::val2str(level)+"_p"+util::val2str(param.curPIdx) + createFileName("_vnorm_xslice", param.curIter, PADDING),
                          *computeVelocityNorm(lattice, xSlice), imSize, imSize);
    writer.writeScaledGif("l"+util::val2str(level)+"_p"+util::val2str(param.curPIdx) + createFileName("_vnorm_yslice", param.curIter, PADDING),
                          *computeVelocityNorm(lattice, ySlice), imSize, imSize);
    writer.writeScaledGif("l"+util::val2str(level)+"_p"+util::val2str(param.curPIdx) + createFileName("_vnorm_zslice", param.curIter, PADDING),
                          *computeVelocityNorm(lattice, zSlice), imSize, imSize);
}


// �ú���׼����ִ��ʵ�ʵķ���
std::auto_ptr<MultiBlockLattice3D<T,DESCRIPTOR> > run ( Param& param, plint level, plb_ofstream& permeablityFile, MultiBlockLattice3D<T,DESCRIPTOR>* iniVal=0 )
{
    // ���ȼ�����Ӳ���
    param.computeLBparameters(level, permeablityFile);

    /*
     * ��ȡPCFSF�ļ���
     */

    // ��PCFSF�ı��漸�Σ���ASCII�������STL�ļ���ʾ����������ŵ��������μ�����ɵ����ݽṹ��
    //   ����DBL��ʾ����˫���ȣ�һ���Ƽ�����ѡ��
    pcout << std::endl << "Reading STL data for the PCFSF geometry." << std::endl;

    TriangleSet<T> triangleSet(param.PCFSF_Geometry_FName, DBL);
    triangleSet.writeBinarySTL(outputDir+"PCFSF_LB.stl");

    // ���¼����ǵ��ʹ��롣�佫�û�������PCFSF���漸��ת��Ϊpalabos�ڲ�ʹ�õġ�
    //   ��Ϊ��Ч�����ݽṹ��������ʹ��TriangleBoundary3D�ṹ��ָ�����ʵı߽�����
    pcout << std::endl << "Constructing TriangleBoundary3D at level " << level << std::endl;
    plint margin = 1;      // �ϰ�����Χ�����cells�Ķ���ԣ��
    plint borderWidth = 1; // ��ΪGuo�߽������ڵ�cell�������á�Ҫ��margin>=borderWidth
    plint blockSize = 0;   // �����0���ʾ����ʹ��ϡ���ʾ

    DEFscaledMesh<T>* defMesh = new DEFscaledMesh<T>(triangleSet, param.resolution, param.flowDirection, margin, param.extraLayer);
    TriangleBoundary3D<T> boundary(*defMesh);
    delete defMesh;
    boundary.getMesh().inflate();

    /*
     * ���ػ�����
     */

    // PCSFS��������ⲿ(���ڲ��෴)�������⡣Ϊ�ˣ�����ʶ�����ڼ�����֮�ڵ�lattice�ڵ�
    //   ���������ⲿ�Ľڵ�������֡���ͨ����������ػ����̴���
    pcout << std::endl << "Running VoxelizedDomain3D at level " << level << std::endl;
    plint envelopeWidth = 2;   // ���Ƶ�off-lattice�߽���������ΪGuo��off lattice�߽�������Ҫ2-cells�ڽ�����
    const int flowType = voxelFlag::outside;
    VoxelizedDomain3D<T> voxelizedDomain ( boundary, flowType, param.extraLayer, borderWidth, envelopeWidth, blockSize );//param.boundingBox()
    if (param.performOutput) {
        pcout << getMultiBlockInfo(voxelizedDomain.getVoxelMatrix()) << std::endl;
    }

    /*
     * ����lattice���ܶȺͶ���blocks.
     */

    pcout << "Generating the lattice." << std::endl;
    std::auto_ptr<MultiBlockLattice3D<T,DESCRIPTOR> > lattice = generateMultiBlockLattice<T,DESCRIPTOR> ( voxelizedDomain.getVoxelMatrix(),
                                                                envelopeWidth, new IncBGKdynamics<T,DESCRIPTOR>(param.omega) ); // ��ģ�����ٶȵ��ڶ���
    lattice->toggleInternalStatistics(false);

    /*
     * ��PCFSF������Bounceback�߽�������outer-domain�߽�����
     */

    pcout << "Generating boundary conditions." << std::endl;
    // �߽������㷨��outer domain.
    for( plint idx=0; idx < plint(param.deltaP.size()); idx++ )
    {
        permeablityFile << std::endl << "Running simulation for the pressure drop: " << param.getPressureDrop(idx) << std::endl;
        permeablityFile << "===========================================================" << std::endl;

        std::auto_ptr< OnLatticeBoundaryCondition3D<T,DESCRIPTOR> > outerBoundaryCondition( createLocalBoundaryCondition3D<T,DESCRIPTOR>() );
        outerDomainBoundaries(param, lattice, voxelizedDomain.getVoxelMatrix(), outerBoundaryCondition);

        if(iniVal) {
            Box3D toDomain(lattice->getBoundingBox());
            Box3D fromDomain(toDomain.shift(margin,margin,margin)); // ����������ʱ��margin�ڳߴ緽�淭�����˴�ͨ��һ��shift��ȡ����һЧ��
            copy(*iniVal, fromDomain, *lattice, toDomain, modif::staticVariables);
        }

        // ��һ��ѡ�����(ѹ��)����ʱ����Ҫ���ValueTracer���Ա�ó������ض�������ϸ�ֲ�Σ��Ѿ��ﵽ�ȶ�̬�Ľ��ۣ���ֹͣ����
        plint convergeIter=20;
        util::ValueTracer<T> converge(0.05*convergeIter, param.resolution, param.epsilon);
        global::timer("iteration").restart();
        lattice->resetTime(param.curIter);
        bool checkForErrors = true;

        pcout << "Saving a " << lattice->getNx() << " by " << lattice->getNy() << " by " << lattice->getNz() << " lattice." << std::endl;
        global::timer("io").start();
        parallelIO::save(*lattice, "checkpoint", false);
        pcout << "Total time for i/o: " << global::timer("io").getTime() << std::endl;

        plb_ofstream iterFile((outputDir+"IterationInfo_p"+util::val2str(idx)+".txt").c_str());
        iterFile << "The iteration information for the input pressure drop: " << param.curDeltaP << std::endl;
        iterFile << "==========================================================================" << std::endl;
        iterFile << "Time\t\t" << "Iteration\t\t" << "Average velocity" << std::endl;

        // ��ײ��������.
        while(!converge.hasConverged() && param.curIter < param.maxIter)
        {
            if (param.curIter % 50==0 && param.performOutput) {
                pcout << "T= " << param.curIter*param.dt << "; " << "Iteration " << param.curIter << std::endl;
                pcout << "Average velocity: "
                    << computeAverage(*computeVelocityComponent(*lattice, lattice->getBoundingBox(), param.flowDirection))*param.dx/param.dt
                    << std::endl;

                iterFile << "Time " << param.curIter*param.dt << "\t\t" << param.curIter << "\t\t"
                    << computeAverage(*computeVelocityComponent(*lattice, lattice->getBoundingBox(), param.flowDirection))*param.dx/param.dt
                    << std::endl;

                pcout << "Writing VTK ..."  << std::endl;
                writeVTK(*lattice, param, level);

                pcout << "Writing Gif image ..." << std::endl;
                writeGifs(*lattice, param, level);
            }

            if (param.curIter % convergeIter==0) {
                converge.takeValue(computeAverage(*computeVelocityComponent(*lattice, lattice->getBoundingBox(), param.flowDirection)),true);//getStoredAverageEnergy(*lattice)
            }

            lattice->collideAndStream();

            if (checkForErrors) {
                abortIfErrorsOccurred();
                checkForErrors = false;
            }
            ++param.curIter;
        }

        computePermeability(lattice, param, permeablityFile);

        pcout << "End of simulation at iteration " << param.curIter << std::endl;
        permeablityFile << "End of simulation at iteration " << param.curIter << std::endl;

        pcout << "Elapsed time: " << global::timer("iteration").stop() << std::endl;
        pcout << "Total elapsed time: " << global::timer("global").getTime() << std::endl;

        permeablityFile << "Elapsed time: " << global::timer("iteration").stop() << std::endl;
        permeablityFile << "Total elapsed time: " << global::timer("global").getTime() << std::endl;

        param.curIter = 0;

        iterFile.close();
    }

    return lattice;
}


int main(int argc, char* argv[])
{
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
        param = Param(paramXmlFileName, permeablityFile);
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
