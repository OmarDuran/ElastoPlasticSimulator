
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
#include "pzpostprocanalysis.h"

// Structure matrix utilities
#include "TPZSSpStructMatrix.h"
#include "pzskylstrmatrix.h"
#include "pzfstrmatrix.h"




/// Read the geometry from a gmsh script
TPZGeoMesh * ReadGeometry();

/// Print the geometry partition
void PrintGeometry(TPZGeoMesh * gmesh);

/// Create computational mesh for computing footing deformation
TPZCompMesh * CMeshFooting2D(TPZGeoMesh * gmesh, int p_order);

/// Analysis that solve the Non-Linear System (NLS)
TPZAnalysis * Analysis(TPZCompMesh * cmesh, int n_threads);

/// Load the solution inside memory once the solution is converged
void AcceptSolution(TPZCompMesh * cmesh, TPZAnalysis * analysis);

/// Non-linear solver
bool FindRoot(TPZAnalysis * analysis);

/// Post-process displacement, strain and stress
void PostProcess(TPZAnalysis *analysis, TPZStack<std::string> & scalnames, std::string plotfile);

/// Apply loading ramp
void LoadingRamp(REAL pseudo_t, TPZCompMesh * cmesh);


/// Global variables
enum EMat_Id { ERock = 1, EBottomBC = 2, ELateralBC = 3, ETopBC = 4, ETopNullBC = 5};
enum EBC_Type { ETn = 1, Eu_null = 3, Eu = 7};

int main(int argc, char *argv[]){
    
    int p_order = 2;
    int n_threads = 8;
    TPZGeoMesh * gmesh = ReadGeometry();
    PrintGeometry(gmesh);
    
    TPZCompMesh *cmesh      = CMeshFooting2D(gmesh, p_order);
    TPZAnalysis *analysis   = Analysis(cmesh,n_threads);
    
    /// Creation of a post-process analysis
    TPZPostProcAnalysis * post_processor = new TPZPostProcAnalysis;
    post_processor->SetCompMesh(cmesh);
    
    int n_regions = 1;
    TPZManVector<int,10> post_mat_id(n_regions);
    post_mat_id[0] = ERock;
    
    TPZStack<std::string> var_names;
    var_names.Push("DisplacementMemX");
    var_names.Push("DisplacementMemY");
    var_names.Push("XStress");
    var_names.Push("YStress");
    var_names.Push("ZStress");
    var_names.Push("VolTotalStrain");
    var_names.Push("PlasticSqJ2");
    var_names.Push("I1Stress");
    
    post_processor->SetPostProcessVariables(post_mat_id, var_names);
    TPZFStructMatrix structmatrix(post_processor->Mesh());
    structmatrix.SetNumThreads(n_threads);
    post_processor->SetStructuralMatrix(structmatrix);
    
    std::string plotfile = "footing.vtk";
    
    // For a single step
    REAL dt = 0.1;
    int n_steps = 10;
    for (int it = 1; it <= n_steps; it++) {
        REAL t = it*dt;
        LoadingRamp(t,cmesh);
        FindRoot(analysis);
        AcceptSolution(cmesh, analysis);
        post_processor->TransferSolution();
        PostProcess(post_processor,var_names,plotfile);
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
    
    // Mohr Coulomb data
    REAL mc_cohesion    = 10.0;
    REAL mc_phi         = (20.0*M_PI/180);
    REAL mc_psi         = mc_phi;
    
    /// ElastoPlastic Material using Mohr Coulomb
    // Elastic predictor
    TPZElasticResponse ER;
    REAL G = 400*mc_cohesion;
    REAL nu = 0.3;
    REAL E = 2.0*G*(1+nu);
    

    
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
//    cmesh->InsertMaterialObject(material);
    
    
    TPZFMatrix<STATE> val1(2,2,0.), val2(2,1,0.);
    
    val2(0,0) = 0;
    val2(1,0) = 1;
    TPZBndCond * bc_bottom = material->CreateBC(material, EBottomBC, Eu_null, val1, val2);
    
    val2(0,0) = 1;
    val2(1,0) = 0;
    TPZBndCond * bc_lateral = material->CreateBC(material, ELateralBC, Eu_null, val1, val2);
    
//    val2(0,0) = 0;
//    val2(1,0) = 0;
//    val1(0,0) = 0;
//    val1(1,1) = 1;
    val2(0,0) = 0;
    val2(1,0) = 0;
    TPZBndCond * bc_top = material->CreateBC(material, ETopBC, ETn, val1, val2);
    
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

TPZAnalysis * Analysis(TPZCompMesh * cmesh, int n_threads){
    
    TPZAnalysis * analysis = new TPZAnalysis(cmesh,true);
//    TPZSkylineStructMatrix matrix(cmesh);
    TPZSymetricSpStructMatrix matrix(cmesh);
    TPZStepSolver<STATE> step;
    step.SetDirect(ELDLt);
    matrix.SetNumThreads(n_threads);
    analysis->SetStructuralMatrix(matrix);
    analysis->SetSolver(step);
    return analysis;
}

void PostProcess(TPZAnalysis *analysis, TPZStack<std::string> & scalnames, std::string plotfile)
{
    const int dim = analysis->Mesh()->Dimension();
    int div = 1;
    TPZStack<std::string> vecnames;
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
    REAL max_uy = 100.0;
    REAL min_uy = 0.0;
    
    /// Compute current displacement
//    REAL uy = (footing_lenght*(max_uy - min_uy)*pseudo_t)/100.0;
    REAL uy = (footing_lenght*(max_uy - min_uy)*pseudo_t);
    
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
    x.Zero();
    REAL tol = 1.0e-4;
    int n_it = 50;
    bool stop_criterion_Q = false;
    REAL norm_res;
    analysis->Assemble();
    for (int i = 1; i <= n_it; i++) {
//        analysis->Solver().Matrix()->SetIsDecomposed(0);
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
        analysis->Assemble();
    }
    
    if (stop_criterion_Q == false) {
        std::cout << "Nonlinear process not converged with residue norm = " << norm_res << std::endl;
    }
    return stop_criterion_Q;
}

void AcceptSolution(TPZCompMesh * cmesh, TPZAnalysis * analysis){
    {
        TPZMaterial *mat = cmesh->FindMaterial(ERock);
        if (!mat) {
            DebugStop();
        }
        
        bool accept_solution_Q = true;
        {
            std::map<int, TPZMaterial *> & refMatVec = cmesh->MaterialVec();
            std::map<int, TPZMaterial * >::iterator mit;
            TPZMatWithMem<TPZElastoPlasticMem> * pMatWithMem;
            
            for(mit=refMatVec.begin(); mit!= refMatVec.end(); mit++)
            {
                pMatWithMem = dynamic_cast<TPZMatWithMem<TPZElastoPlasticMem> *>( mit->second );
                if(pMatWithMem != NULL)
                {
                    pMatWithMem->SetUpdateMem(accept_solution_Q);
                }
            }
        }
        
        analysis->AssembleResidual();
        
        accept_solution_Q = false;
        {
            std::map<int, TPZMaterial *> & refMatVec = cmesh->MaterialVec();
            std::map<int, TPZMaterial * >::iterator mit;
            TPZMatWithMem<TPZElastoPlasticMem> * pMatWithMem;
            
            for(mit=refMatVec.begin(); mit!= refMatVec.end(); mit++)
            {
                pMatWithMem = dynamic_cast<TPZMatWithMem<TPZElastoPlasticMem> *>( mit->second );
                if(pMatWithMem != NULL)
                {
                    pMatWithMem->SetUpdateMem(accept_solution_Q);
                }
            }
        }
    }
}
