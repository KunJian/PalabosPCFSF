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
    // xml输入参数
    //----------------------------------
    std::string PCFSFDataFile;

    T lx, ly, lz;                       // 计算域的尺寸，物理单位

    T nu;                               // Kinematic viscosity，物理单位
    T rho;                              // 流体密度，物理单位

    std::vector<T> deltaP;              // 压降，格子单位
    plint flowDirection;                // 入口速度方向
    plint referenceResolution;          // 沿参考长度的格子数

    plint curIter, maxIter;             // 当前和最大迭代次数

    plint maxLevel;
    T     epsilon;

    bool performOutput;
    bool doImages;

    // 计算的参数
    //----------------------------------
    plint cxLB, cyLB, czLB;             // PCFSF中心的位置，格子单位
    plint nx, ny, nz;                   // 包围盒的网格分辨率
    Box3D inlet, outlet;                // 外部区域边界，格子单位
    Box3D lateral1, lateral2, lateral3, lateral4;
    T     nuLB;                         // Kinematic viscosity，格子单位

    T     curDeltaP;                    // 压降，格子单位
    plint curPIdx;                      // 当前压降的索引
    plint resolution;                   // 不同细分层次的实际格子数
    T     omega;                        // 松弛参数
    T     dx, dt;                       // 离散空间和时间步

    Param() {}

    Param(std::string xmlFname, plb_ofstream& permeablityFile)
    {
        XMLreader document(xmlFname);
        document["geometry"]["data"].read(PCFSFDataFile);
        document["geometry"]["domain"]["x"].read(lx);
        document["geometry"]["domain"]["y"].read(ly);
        document["geometry"]["domain"]["z"].read(lz);

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

        // 分辨率层次增加1，每个坐标方向的分辨率则增加1倍
        //   根据定义，参数``referenceResolution''是网格细分层次为0时的分辨率
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

// 该函数对象定义一个零速度，以及一个沿流体流动方向线性减小的压力。用于初始化粒子群。
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

void readGeometry(Param& param, MultiScalarField3D<int>& geometry)
{
    const plint nx = geometry.getNx();
    const plint ny = geometry.getNy();
    const plint nz = geometry.getNz();

    //Box3D sliceBox(0, nx-1, 0, ny-1, 0, 0);
    Box3D sliceBox(0, 0, 0, ny-1, 0, nz-1);
    std::auto_ptr<MultiScalarField3D<int> > slice = generateMultiScalarField<int>(geometry, sliceBox);

    plb_ifstream geometryFile( param.PCFSFDataFile.c_str() );

    //for (plint iZ=0; iZ<nz-1; ++iZ) {
    for (plint iX=0; iX<nx-1; iX++) {
        if (!geometryFile.is_open()) {
            pcout << "Error: could not open geometry file " << param.PCFSFDataFile << std::endl;
            exit(EXIT_FAILURE);
        }
        geometryFile >> *slice;
        //copy(*slice, slice->getBoundingBox(), geometry, Box3D(0, nx-1, 0, ny-1, iZ, iZ));
        copy(*slice, slice->getBoundingBox(), geometry, Box3D(iX, iX, 0, ny-1, 0, nz-1));
    }

    {
        VtkImageOutput3D<T> vtkOut("PCFSFGeom", param.dx);
        vtkOut.writeData<float>(*copyConvert<int,T>(geometry, geometry.getBoundingBox()), "tag", 1.0);
    }

    {
        std::auto_ptr<MultiScalarField3D<T> > floatTags = copyConvert<int,T>(geometry, geometry.getBoundingBox());
        std::vector<T> isoLevels;
        isoLevels.push_back(0.5);
        typedef TriangleSet<T>::Triangle Triangle;
        std::vector<Triangle> triangles;
        Box3D domain = floatTags->getBoundingBox().enlarge(-1);
        domain.x0++;
        domain.x1--;
        isoSurfaceMarchingCube(triangles, *floatTags, isoLevels, domain);
        TriangleSet<T> set(triangles);
        set.writeBinarySTL(outputDir + "PCFSF_LBGeom.stl");
    }
}

// 实例化outer domain的边界条件
void outerDomainBoundaries(MultiBlockLattice3D<T,DESCRIPTOR>& lattice, MultiScalarField3D<int>& geometry,
                           std::auto_ptr< OnLatticeBoundaryCondition3D<T,DESCRIPTOR> >& bc, Param& param)
{
    /*
     * 设置初始条件
     */
    // 初始条件：所有地方都是常压力和velocity-at-infinity
    pcout << "Define boundary conditions." << std::endl;
    if( param.flowDirection == 0 )
    {
        bc->addPressureBoundary0N(param.inlet, lattice);
        setBoundaryDensity(lattice, param.inlet, (T) 1.);

        bc->addPressureBoundary0P(param.outlet, lattice);
        setBoundaryDensity(lattice, param.outlet, (T) 1. - param.curDeltaP*DESCRIPTOR<T>::invCs2);
    }
    else if( param.flowDirection == 1 )
    {
        bc->addPressureBoundary1N(param.inlet, lattice);
        setBoundaryDensity(lattice, param.inlet, (T) 1.);

        bc->addPressureBoundary1P(param.outlet, lattice);
        setBoundaryDensity(lattice, param.outlet, (T) 1. - param.curDeltaP*DESCRIPTOR<T>::invCs2);
    }
    else
    {
        bc->addPressureBoundary2N(param.inlet, lattice);
        setBoundaryDensity(lattice, param.inlet, (T) 1.);

        bc->addPressureBoundary2P(param.outlet, lattice);
        setBoundaryDensity(lattice, param.outlet, (T) 1. - param.curDeltaP*DESCRIPTOR<T>::invCs2);
    }

    // PCFSF的外边界(标记为1)、以及流体域的四个壁面，采用bounce-back边界条件
    defineDynamics(lattice, geometry, new BounceBack<T,DESCRIPTOR>(), 1);
    defineDynamics(lattice, param.lateral1, new BounceBack<T,DESCRIPTOR>((T)1.));
    defineDynamics(lattice, param.lateral2, new BounceBack<T,DESCRIPTOR>((T)1.));
    defineDynamics(lattice, param.lateral3, new BounceBack<T,DESCRIPTOR>((T)1.));
    defineDynamics(lattice, param.lateral4, new BounceBack<T,DESCRIPTOR>((T)1.));

    // PCFSF的内部(标记为2)，采用no-dynamics（什么也不做）
    defineDynamics(lattice, geometry, new NoDynamics<T,DESCRIPTOR>(), 2);

    pcout << "Initilization of rho and u." << std::endl;
    initializeAtEquilibrium( lattice, lattice.getBoundingBox(), PressureGradient(param) );

    lattice.initialize();
}

void computePermeability(MultiBlockLattice3D<T,DESCRIPTOR>& lattice, Param& param, plb_ofstream& permeablityFile)
{
    pcout << "Computing the permeability." << std::endl;

    // 计算流动方向的平均速度
    T meanVecLB = computeAverage(*computeVelocityComponent(lattice, lattice.getBoundingBox(), param.flowDirection));
    T inletVecLB = computeAverage(*computeVelocityComponent(lattice, param.inlet, param.flowDirection));
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

// 保存围绕PCFSF的流动的VTK文件，可用Paraview查看
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
    vtkSliceOut.writeData<float>( *computeVelocityNorm(lattice, ySlice), "velocityNorm", param.dx/param.dt );
    vtkSliceOut.writeData<float>( *computeVelocityNorm(lattice, zSlice), "velocityNorm", param.dx/param.dt );

    vtkSliceOut.writeData<3,float>(*computeVelocity(lattice, xSlice), "velocity", param.dx/param.dt );
    vtkSliceOut.writeData<3,float>(*computeVelocity(lattice, ySlice), "velocity", param.dx/param.dt );
    vtkSliceOut.writeData<3,float>(*computeVelocity(lattice, zSlice), "velocity", param.dx/param.dt );

    vtkSliceOut.writeData<float>( *computeDensity(lattice, xSlice), "pressure", util::sqr(param.dx/param.dt)*param.rho );
    vtkSliceOut.writeData<float>( *computeDensity(lattice, ySlice), "pressure", util::sqr(param.dx/param.dt)*param.rho );
    vtkSliceOut.writeData<float>( *computeDensity(lattice, zSlice), "pressure", util::sqr(param.dx/param.dt)*param.rho );
}

// 保存中心切片的Gif图
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

// 该函数准备并执行实际的仿真
MultiBlockLattice3D<T,DESCRIPTOR> run ( Param& param, plint level, plb_ofstream& permeablityFile, MultiBlockLattice3D<T,DESCRIPTOR>* iniVal=0 )
{
    // 计算格子参数
    param.computeLBparameters(level, permeablityFile);

    /*
     * 生成lattice
     */
    pcout << "Generating the lattice." << std::endl;
    MultiBlockLattice3D<T,DESCRIPTOR> lattice(param.nx, param.ny, param.nz, new IncBGKdynamics<T,DESCRIPTOR>(param.omega)); // 该模型中速度等于动量
    lattice.toggleInternalStatistics(false);
    lattice.periodicity().toggleAll(false);

    pcout << "Reading the geometry file." << std::endl;
    MultiScalarField3D<int> geometry(param.nx, param.ny, param.nz);
    readGeometry(param, geometry);

    /*
     * 在PCFSF上生成Bounceback边界条件和outer-domain边界条件
     */
    pcout << "Generating boundary conditions." << std::endl;
    for( plint idx=0; idx < plint(param.deltaP.size()); idx++ )
    {
        permeablityFile << std::endl << "Running simulation for the pressure drop: " << param.getPressureDrop(idx) << std::endl;
        permeablityFile << "===========================================================" << std::endl;

        // 边界条件算法
        std::auto_ptr< OnLatticeBoundaryCondition3D<T,DESCRIPTOR> > outerBoundaryCondition( createLocalBoundaryCondition3D<T,DESCRIPTOR>() );
        outerDomainBoundaries(lattice, geometry, outerBoundaryCondition, param);

        if(iniVal) {
            Box3D toDomain(lattice.getBoundingBox());

            plint margin = 1;      // 障碍物周围分配的cells的额外裕度
            Box3D fromDomain(toDomain.shift(margin,margin,margin)); // 在重新缩放时，margin在尺寸方面翻倍，此处通过一个shift来取消这一效果
            copy(*iniVal, fromDomain, lattice, toDomain, modif::staticVariables);
        }

        // 当一个选择的量(压力)收敛时，需要检查ValueTracer，以便得出对于特定的网格细分层次，已经达到稳定态的结论，并停止仿真
        plint convergeIter=20;
        util::ValueTracer<T> converge(0.05*convergeIter, param.resolution, param.epsilon);
        global::timer("iteration").restart();
        lattice.resetTime(param.curIter);
        bool checkForErrors = true;

        pcout << "Saving a " << lattice.getNx() << " by " << lattice.getNy() << " by " << lattice.getNz() << " lattice." << std::endl;
        global::timer("io").start();
        parallelIO::save(lattice, "checkpoint", false);
        pcout << "Total time for i/o: " << global::timer("io").getTime() << std::endl;

        plb_ofstream iterFile((outputDir+"IterationInfo_p"+util::val2str(idx)+".txt").c_str());
        iterFile << "The iteration information for the input pressure drop: " << param.curDeltaP << std::endl;
        iterFile << "==========================================================================" << std::endl;
        iterFile << "Time\t\t" << "Iteration\t\t" << "Average velocity" << std::endl;

        // 碰撞和流迭代.
        while(!converge.hasConverged() && param.curIter < param.maxIter)
        {
            if (param.curIter % 100==0 && param.performOutput) {
                 iterFile << "Time " << param.curIter*param.dt << "\t\t" << param.curIter << "\t\t"
                    << computeAverage(*computeVelocityComponent(lattice, lattice.getBoundingBox(), param.flowDirection))*param.dx/param.dt
                    << std::endl;
            }

            if (param.curIter % 200==0) {
                pcout << "Writing Gif image ..." << std::endl;
                writeGifs(lattice, param, level);
            }

            if (param.curIter % 500==0) {
                pcout << "Writing VTK ..."  << std::endl;
                writeVTK(lattice, param, level);
            }

            if (param.curIter % convergeIter==0) {
                converge.takeValue(computeAverage(*computeVelocityComponent(lattice, lattice.getBoundingBox(), param.flowDirection)),true);//getStoredAverageEnergy(*lattice)
            }

            lattice.collideAndStream();

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

    // 读取参数XML输入文件(里面也包含了大量的注释).
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
    // 以下代码结合了平滑网格细分直到收敛的概念。平滑表示细分层次增加1，整个网格在每个方向上翻倍。
    //   当网格细化后，dx和dt都要改变。dt是按dx^2 (扩散行为)变化
    try {
        for (plint level=iniLevel; level<=param.maxLevel; ++level) {
            pcout << std::endl << "Running new simulation at level " << level << std::endl;
            permeablityFile << std::endl << "Running new simulation at level " << level << std::endl;

            MultiBlockLattice3D<T,DESCRIPTOR> convergedLattice( run(param, level, permeablityFile, iniConditionLattice.get()) );
            if (level != param.maxLevel) {
                plint dxScale = -1;
                plint dtScale = -2;

                // 前一个网格层次收敛的仿真是下一个网格层次（在进行适当的差值之后）仿真的初始条件
                iniConditionLattice = std::auto_ptr<MultiBlockLattice3D<T,DESCRIPTOR> >(refine(convergedLattice, dxScale, dtScale,
                                                                                               new IncBGKdynamics<T,DESCRIPTOR>(param.omega)) );
            }
        }

        permeablityFile.close();
    }
    catch(PlbException& exception) {
        pcout << exception.what() << std::endl;
        return -1;
    }
}
