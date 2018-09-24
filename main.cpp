
#include <stdio.h>
#include <iostream>

// Geometry utilities
#include "pzgmesh.h"
#include "TPZGmshReader.h"
#include "TPZVTKGeoMesh.h"

// Material utilities
#include "pzbndcond.h"
#include "TPZMatElastoPlastic2D.h"
#include "TPZMatElastoPlastic.h"
#include "TPZElastoPlasticMem.h"
// Elasticity
#include "TPZElasticCriterion.h"
// Plasticity
#include "TPZPlasticStepPV.h"
#include "TPZSandlerExtended.h"
#include "TPZYCMohrCoulombPV.h"

// Analysis utilities
#include "pzanalysis.h"
#include "pzstepsolver.h"

// Structure matrix utilities
#include "TPZSSpStructMatrix.h"
#include "pzskylstrmatrix.h"




/// Read the geometry from a gmsh script
TPZGeoMesh * ReadGeometry();

/// Print the geometry partition
void PrintGeometry(TPZGeoMesh * gmesh);

/// Create computational mesh for computing footing deformation
TPZCompMesh * CMeshFooting2D(TPZGeoMesh * gmesh, int p_order);

/// Analysis that solve the Non-Linear System (NLS)
TPZAnalysis * Analysis(TPZCompMesh * cmesh);

/// Non-linear solver
bool FindRoot(TPZAnalysis * analysis);

/// Post-process displacement, strain and stress
void PostProcess(TPZAnalysis *analysis, std::string plotfile);

/// Apply loading ramp
void LoadingRamp(REAL pseudo_t, TPZCompMesh * cmesh);


/// Global variables
enum EMat_Id { ERock = 1, EBottomBC = 2, ELateralBC = 3, ETopBC = 4, ETopNullBC = 5};
enum EBC_Type { ETn = 1, Eu_null = 3, Eu = 7};

int main(int argc, char *argv[]){
    
    int p_order = 1;
    TPZGeoMesh * gmesh = ReadGeometry();
    PrintGeometry(gmesh);
    
    TPZCompMesh *cmesh      = CMeshFooting2D(gmesh, p_order);
    TPZAnalysis *analysis   = Analysis(cmesh);
    std::string plotfile = "footing.vtk";
    
    // For a single step
    REAL dt = 0.1;
    int n_steps = 10;
    for (int it = 1; it <= n_steps; it++) {
        REAL t = it*dt;
        LoadingRamp(t,cmesh);
        FindRoot(analysis);
        PostProcess(analysis,plotfile);
    }
    
    std::cout << "Execution complete." << std::endl;
	return 0;
}

TPZGeoMesh * ReadGeometry(){
    
    TPZGeoMesh * gmesh = new TPZGeoMesh;
    std::string grid("footing_problem.msh");
    
    TPZGmshReader Geometry;
    REAL s = 1.0;
    Geometry.SetfDimensionlessL(s);
    gmesh = Geometry.GeometricGmshMesh(grid);
    const std::string name("Footing section");
    gmesh->SetName(name);
    
    return gmesh;
}

void PrintGeometry(TPZGeoMesh * gmesh){
    std::stringstream text_name;
    std::stringstream vtk_name;
    text_name   << "geometry" << ".txt";
    vtk_name    << "geometry" << ".vtk";
    std::ofstream textfile(text_name.str().c_str());
    gmesh->Print(textfile);
    std::ofstream vtkfile(vtk_name.str().c_str());
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, vtkfile, true);
}

TPZCompMesh * CMeshFooting2D(TPZGeoMesh * gmesh, int p_order){
    
    unsigned int dim  = gmesh->Dimension();
    const std::string name("ElastoPlastic Footing Problem ");

    // Setting up attributes
    TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
    cmesh->SetName(name);
    cmesh->SetDefaultOrder(p_order);
    cmesh->SetDimModel(dim);
    
    /// ElastoPlastic Material using Mohr Coulomb
    // Elastic predictor
    TPZElasticResponse ER;
    REAL G = 12000.0;
    REAL nu = 0.3;
    REAL E = 2.0*G*(1+nu);
    
    // Mohr Coulomb data
    REAL mc_cohesion    = 25.0*1.0e10;
    REAL mc_phi         = (90.0*M_PI/180);//(90.0*M_PI/180);//(10.0*M_PI/180);
    REAL mc_psi         = mc_phi; // because MS do not understand
    
    TPZPlasticStepPV<TPZYCMohrCoulombPV, TPZElasticResponse> LEMC;
    ER.SetUp(E, nu);
    LEMC.SetElasticResponse(ER);
    LEMC.fYC.SetUp(mc_phi, mc_psi, mc_cohesion, ER);
    int PlaneStrain = 1;
    
    TPZElasticCriterion MatEla;
    MatEla.SetElasticResponse(ER);
    
    TPZMatElastoPlastic2D < TPZElasticCriterion, TPZElastoPlasticMem > * material = new TPZMatElastoPlastic2D < TPZElasticCriterion, TPZElastoPlasticMem >(ERock,PlaneStrain);
    material->SetPlasticity(MatEla);
    cmesh->InsertMaterialObject(material);
    
//    TPZMatElastoPlastic2D < TPZPlasticStepPV<TPZYCMohrCoulombPV, TPZElasticResponse>, TPZElastoPlasticMem > * material = new TPZMatElastoPlastic2D < TPZPlasticStepPV<TPZYCMohrCoulombPV, TPZElasticResponse>, TPZElastoPlasticMem >(ERock,PlaneStrain);
//    material->SetPlasticity(LEMC);
//    material->SetPlasticity(LEMC);
//    cmesh->InsertMaterialObject(material);
    
    
    TPZFMatrix<STATE> val1(2,2,0.), val2(2,1,0.);
    
    val2(0,0) = 0;
    val2(1,0) = 1;
    TPZBndCond * bc_bottom = material->CreateBC(material, EBottomBC, Eu_null, val1, val2);
    
    val2(0,0) = 1;
    val2(1,0) = 0;
    TPZBndCond * bc_lateral = material->CreateBC(material, ELateralBC, Eu_null, val1, val2);
    
    val2(0,0) = 0;
    val2(1,0) = 0;
    val1(0,0) = 0;
    val1(1,1) = 1;
    TPZBndCond * bc_top = material->CreateBC(material, ETopBC, Eu, val1, val2);
    
    val2(0,0) = 0;
    val2(1,0) = 0;
    TPZBndCond * bc_top_null = material->CreateBC(material, ETopNullBC, ETn, val1, val2);
    
    cmesh->InsertMaterialObject(bc_bottom);
    cmesh->InsertMaterialObject(bc_lateral);
    cmesh->InsertMaterialObject(bc_top);
//    cmesh->InsertMaterialObject(bc_top_null);
    
    cmesh->SetAllCreateFunctionsContinuousWithMem();
    cmesh->AutoBuild();
    
#ifdef PZDEBUG
    std::ofstream out("cmesh.txt");
    cmesh->Print(out);
#endif
    
    return cmesh;
}

TPZAnalysis * Analysis(TPZCompMesh * cmesh){
    
    int numofThreads = 0;
    TPZAnalysis * analysis = new TPZAnalysis(cmesh,true);
//    TPZSkylineStructMatrix matrix(cmesh);
    TPZSymetricSpStructMatrix matrix(cmesh);
    TPZStepSolver<STATE> step;
    step.SetDirect(ELDLt);
    matrix.SetNumThreads(numofThreads);
    analysis->SetStructuralMatrix(matrix);
    analysis->SetSolver(step);
    return analysis;
}

void PostProcess(TPZAnalysis *analysis, std::string plotfile)
{
    const int dim = analysis->Mesh()->Dimension();
    int div = 0;
    TPZStack<std::string> scalnames, vecnames;
//    scalnames.Push("SigmaX");
//    scalnames.Push("SigmaY");
//    scalnames.Push("SigmaZ");
    vecnames.Push("Displacement");
    analysis->DefineGraphMesh(dim, scalnames, vecnames, plotfile);
    analysis->PostProcess(div);
}

void LoadingRamp(REAL pseudo_t, TPZCompMesh * cmesh){
    
    if (!cmesh) {
        DebugStop();
    }
    
    TPZMaterial *mat = cmesh->FindMaterial(ERock);
    if (!mat) {
        DebugStop();
    }
    
    /// Compute the footing lenght
    REAL footing_lenght = 0;
    {
        TPZGeoMesh * gmesh = cmesh->Reference();
        if (!gmesh) {
            DebugStop();
        }
        int n_el = gmesh ->NElements();
        for (int iel = 0; iel < n_el; iel++) {
            TPZGeoEl * gel = gmesh->Element(iel);
            if (!gel) {
                DebugStop();
            }
            
            if (gel->MaterialId() != ETopBC) {
                continue;
            }
            REAL gel_length = gel->SideArea(gel->NSides() - 1);
            footing_lenght += gel_length;
        }
    }
    
    /// Apply loading
    REAL max_uy = 5.0;
    REAL min_uy = 0.0;
    
    /// Compute current displacement
    REAL uy = footing_lenght*(max_uy - min_uy)*pseudo_t/100.0;
    
    /// Apply current displacement
    TPZFMatrix<STATE> val2(2,1,0.);
    val2(1,0) = -uy;
    TPZBndCond * bc_top = NULL;
    bc_top = dynamic_cast<TPZBndCond *> (cmesh->FindMaterial(ETopBC));
    if (!bc_top || bc_top->Material() != mat) {
        DebugStop();
    } else {
        bc_top->Val2() = val2;
    }
    
    
}

bool FindRoot(TPZAnalysis *analysis){
    
    TPZFMatrix<STATE> x(analysis->Solution()), dx;
//    x.Zero();
    REAL tol = 1.0e2;
    int n_it = 20;
    bool stop_criterion_Q = false;
    REAL norm_res;
    for (int i = 1; i <= n_it; i++) {
        analysis->Assemble();
//        analysis->Rhs() *= -1.0;
        analysis->Solve();
        dx = analysis->Solution();
        x += dx;
        analysis->LoadSolution(x);
        analysis->AssembleResidual();
        norm_res = Norm(analysis->Rhs());
        stop_criterion_Q = norm_res < tol;
        if (stop_criterion_Q) {
            std::cout << "Nonlinear process converged with residue norm = " << norm_res << std::endl;
            std::cout << "Number of iterations = " << i << std::endl;
            break;
        }
    }
    
    if (stop_criterion_Q == false) {
        std::cout << "Nonlinear process not converged with residue norm = " << norm_res << std::endl;
    }
    return stop_criterion_Q;
}
