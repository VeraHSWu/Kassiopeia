/*
 * KKrylovChargeDensitySolver.hh
 *
 *  Created on: 12 Aug 2015
 *      Author: wolfgang
 */

#ifndef KEMFIELD_SOURCE_2_0_CHARGEDENSITYSOLVERS_ELECTRIC_INCLUDE_KKRYLOVCHARGEDENSITYSOLVER_HH_
#define KEMFIELD_SOURCE_2_0_CHARGEDENSITYSOLVERS_ELECTRIC_INCLUDE_KKRYLOVCHARGEDENSITYSOLVER_HH_

#include "KChargeDensitySolver.hh"
#include "KElectrostaticBasis.hh"
#include "KKrylovSolverConfiguration.hh"

namespace KEMField {

template< typename ValueType >
class KBoundaryMatrixGenerator; // forward declaration

class KKrylovChargeDensitySolver : public KChargeDensitySolver {
public:
	typedef KElectrostaticBasis::ValueType ValueType;
	typedef KBoundaryMatrixGenerator<ValueType> MatrixGenerator;

	KKrylovChargeDensitySolver();
	virtual ~KKrylovChargeDensitySolver();

	void SetMatrixGenerator(KSmartPointer<MatrixGenerator > matrixGen);
	KSmartPointer<const MatrixGenerator> GetMatrixGenerator() const;
	void SetPreconditionerGenerator(KSmartPointer<MatrixGenerator > preconGen);

	void SetIterationsBetweenRestart(unsigned int iterationsBetweenRestart) {
		fKrylovConfig.SetIterationsBetweenRestart(iterationsBetweenRestart);
	}

	void SetMaxIterations(unsigned int maxIterations) {
		fKrylovConfig.SetMaxIterations(maxIterations);
	}

	void SetSolverName(const std::string& solverName) {
		fKrylovConfig.SetSolverName(solverName);
	}

	void SetStepsBetweenCheckpoints(unsigned int stepsBetweenCheckpoints) {
		fKrylovConfig.SetStepsBetweenCheckpoints(stepsBetweenCheckpoints);
	}

	void SetStepsBetweenTimeChecks(unsigned int stepsBetweenTimeChecks) {
		fKrylovConfig.SetStepsBetweenTimeChecks(stepsBetweenTimeChecks);
	}

	void SetTimeLimitSeconds(double timeLimitSeconds) {
		fKrylovConfig.SetTimeLimitSeconds(timeLimitSeconds);
	}

	void SetTolerance(double tolerance) {
		fKrylovConfig.SetTolerance(tolerance);
	}

	void SetUseCheckpoints(bool useCheckpoints) {
		fKrylovConfig.SetUseCheckpoints(useCheckpoints);
	}

	void SetUseDisplay(bool useDisplay) {
		fKrylovConfig.SetUseDisplay(useDisplay);
	}

	void SetUsePlot(bool usePlot) {
		fKrylovConfig.SetUsePlot(usePlot);
	}

	void SetUseTimer(bool useTimer) {
		fKrylovConfig.SetUseTimer(useTimer);
	}

private:
	void InitializeCore(KSurfaceContainer& container);
	void ComputeSolution(KSurfaceContainer& container);

	KSmartPointer<MatrixGenerator> fMatrixGenerator;
	KSmartPointer<MatrixGenerator> fPreconditionerGenerator;

	KKrylovSolverConfiguration fKrylovConfig;
};

} /* namespace KEMField */

#endif /* KEMFIELD_SOURCE_2_0_CHARGEDENSITYSOLVERS_ELECTRIC_INCLUDE_KKRYLOVCHARGEDENSITYSOLVER_HH_ */