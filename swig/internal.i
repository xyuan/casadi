%exception  casadi::Adaptor< Derived, Solver >::addOptions() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Adaptor< DleToDple , DpleInternal  >::addOptions() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Adaptor< DleToLrDle , LrDleInternal  >::addOptions() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Adaptor< DpleToLrDple , LrDpleInternal  >::addOptions() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Adaptor< LpToQp , QpSolverInternal  >::addOptions() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Adaptor< LrDleToDle , DleInternal  >::addOptions() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Adaptor< LrDpleToDple , DpleInternal  >::addOptions() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Adaptor< LrDpleToLrDle , LrDleInternal  >::addOptions() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Adaptor< QcqpToSocp , SocpSolverInternal  >::addOptions() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Adaptor< QpToImplicit , NlpSolverInternal  >::addOptions() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Adaptor< QpToNlp , NlpSolverInternal  >::addOptions() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Adaptor< QpToQcqp , QcqpSolverInternal  >::addOptions() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Adaptor< SdqpToSdp , SdpSolverInternal  >::addOptions() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Adaptor< SocpToSdp , SdpSolverInternal  >::addOptions() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Assertion::clone() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Assertion::eval(const cpv_MX &input, const pv_MX &output) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Assertion::evalAdj(const std::vector< pv_MX > &adjSeed, const std::vector< pv_MX > &adjSens) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Assertion::evalFwd(const std::vector< cpv_MX > &fwdSeed, const std::vector< pv_MX > &fwdSens) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Assertion::getOp() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Assertion::printPart(std::ostream &stream, int part) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::BinaryMX< ScX, ScY >::clone() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::BinaryMX< ScX, ScY >::eval(const cpv_MX &arg, const pv_MX &res) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::BinaryMX< ScX, ScY >::evalAdj(const std::vector< pv_MX > &aseed, const std::vector< pv_MX > &asens) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::BinaryMX< ScX, ScY >::evalFwd(const std::vector< cpv_MX > &fseed, const std::vector< pv_MX > &fsens) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::BinaryMX< ScX, ScY >::generate(std::ostream &stream, const std::vector< int > &arg, const std::vector< int > &res, CodeGenerator &gen) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::BinaryMX< ScX, ScY >::getBinary(int op, const MX &y, bool scX, bool scY) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::BinaryMX< ScX, ScY >::getOp() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::BinaryMX< ScX, ScY >::getUnary(int op) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::BinaryMX< ScX, ScY >::isBinaryOp() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::BinaryMX< ScX, ScY >::numInplace() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::BinaryMX< ScX, ScY >::printPart(std::ostream &stream, int part) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::BinarySX::dep(int i) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::BinarySX::dep(int i) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::BinarySX::getOp() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::BinarySX::hasDep() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::BinarySX::isSmooth() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::BinarySX::ndep() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::BinarySX::print(std::ostream &stream, long &remaining_calls) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CallFunction::clone() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CallFunction::eval(const cpv_MX &arg, const pv_MX &res) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CallFunction::evalAdj(const std::vector< pv_MX > &adjSeed, const std::vector< pv_MX > &adjSens) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CallFunction::evalFwd(const std::vector< cpv_MX > &fwdSeed, const std::vector< pv_MX > &fwdSens) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CallFunction::generate(std::ostream &stream, const std::vector< int > &arg, const std::vector< int > &res, CodeGenerator &gen) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CallFunction::getFunction() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CallFunction::getFunctionInput() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CallFunction::getFunctionOutput() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CallFunction::getOp() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CallFunction::nTmp(size_t &ni, size_t &nr) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CallFunction::nout() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CallFunction::printPart(std::ostream &stream, int part) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CallFunction::sparsity(int oind) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CodeGenerator::addAuxiliary(Auxiliary f) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CodeGenerator::addDependency(const Function &f) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CodeGenerator::addInclude(const std::string &new_include, bool relative_path=false) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CodeGenerator::addSparsity(const Sparsity &sp) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CodeGenerator::casadi_dot(int n, const std::string &x, int inc_x, const std::string &y, int inc_y) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CodeGenerator::copyVector(std::ostream &s, const std::string &arg, std::size_t n, const std::string &res, const std::string &it="i", bool only_if_exists=false) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CodeGenerator::flush(std::ostream &s) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CodeGenerator::getConstant(const std::vector< double > &v, bool allow_adding=false) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CodeGenerator::getConstant(const std::vector< int > &v, bool allow_adding=false) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CodeGenerator::getDependency(const Function &f) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CodeGenerator::getSparsity(const Sparsity &sp) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CollocationIntegrator::calculateInitialConditions() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CollocationIntegrator::calculateInitialConditionsB() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CollocationIntegrator::clone() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CollocationIntegrator::create(const Function &f, const Function &g) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CollocationIntegrator::init() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CollocationIntegrator::setupFG() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Concat::generate(std::ostream &stream, const std::vector< int > &arg, const std::vector< int > &res, CodeGenerator &gen) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Concat::getGetNonzeros(const Sparsity &sp, const std::vector< int > &nz) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Constant< Value >::clone() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Constant< Value >::generate(std::ostream &stream, const std::vector< int > &arg, const std::vector< int > &res, CodeGenerator &gen) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Constant< Value >::getBinary(int op, const MX &y, bool ScX, bool ScY) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Constant< Value >::getGetNonzeros(const Sparsity &sp, const std::vector< int > &nz) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Constant< Value >::getHorzcat(const std::vector< MX > &x) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Constant< Value >::getMatrixValue() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Constant< Value >::getReshape(const Sparsity &sp) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Constant< Value >::getSetNonzeros(const MX &y, const std::vector< int > &nz) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Constant< Value >::getSetSparse(const Sparsity &sp) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Constant< Value >::getTranspose() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Constant< Value >::getUnary(int op) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Constant< Value >::getValue() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Constant< Value >::getVertcat(const std::vector< MX > &x) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Constant< Value >::isIdentity() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Constant< Value >::isOne() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Constant< Value >::isValue(double val) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Constant< Value >::isZero() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Constant< Value >::printPart(std::ostream &stream, int part) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ConstantDMatrix::clone() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ConstantDMatrix::generate(std::ostream &stream, const std::vector< int > &arg, const std::vector< int > &res, CodeGenerator &gen) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ConstantDMatrix::getMatrixValue() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ConstantDMatrix::getValue() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ConstantDMatrix::isIdentity() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ConstantDMatrix::isMinusOne() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ConstantDMatrix::isOne() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ConstantDMatrix::isZero() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ConstantDMatrix::printPart(std::ostream &stream, int part) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ConstantMX::clone() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ConstantMX::eval(const cpv_MX &input, const pv_MX &output) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ConstantMX::evalAdj(const std::vector< pv_MX > &adjSeed, const std::vector< pv_MX > &adjSens) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ConstantMX::evalFwd(const std::vector< cpv_MX > &fwdSeed, const std::vector< pv_MX > &fwdSens) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ConstantMX::getInnerProd(const MX &y) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ConstantMX::getMatrixValue() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ConstantMX::getOp() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ConstantMX::getValue() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ConstantSX::getOp() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ConstantSX::getValue() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ConstantSX::isConstant() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CplexInterface::clone() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CplexInterface::evaluate() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CplexInterface::freeCplex() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CplexInterface::init() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CsparseInterface::clone() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CsparseInterface::init() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CsparseInterface::prepare() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CvodesInterface::clone() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CvodesInterface::create(const Function &f, const Function &g) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CvodesInterface::freeCVodes() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CvodesInterface::getJac() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CvodesInterface::getJacB() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CvodesInterface::getJacGen() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CvodesInterface::getJacGenB() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CvodesInterface::init() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CvodesInterface::initAdj() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CvodesInterface::integrate(double t_out) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CvodesInterface::integrateB(double t_out) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CvodesInterface::printStats(std::ostream &stream) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CvodesInterface::reset() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CvodesInterface::resetB() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CvodesInterface::setStopTime(double tf) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DenseIO< Derived >::init() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DenseIO< Derived >::inputD(int i) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DenseIO< Derived >::inputD(int i) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DenseIO< Derived >::outputD(int i) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DenseIO< Derived >::outputD(int i) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DenseIO< Derived >::readInputs() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DenseIO< Derived >::writeOutputs() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DenseMultiplication::clone() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DenseMultiplication::generate(std::ostream &stream, const std::vector< int > &arg, const std::vector< int > &res, CodeGenerator &gen) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DenseTranspose::clone() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DenseTranspose::generate(std::ostream &stream, const std::vector< int > &arg, const std::vector< int > &res, CodeGenerator &gen) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DenseTranspose::nTmp(size_t &ni, size_t &nr) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Determinant::clone() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Determinant::eval(const cpv_MX &input, const pv_MX &output) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Determinant::evalAdj(const std::vector< pv_MX > &adjSeed, const std::vector< pv_MX > &adjSens) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Determinant::evalFwd(const std::vector< cpv_MX > &fwdSeed, const std::vector< pv_MX > &fwdSens) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Determinant::getOp() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Determinant::printPart(std::ostream &stream, int part) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Diagcat::clone() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Diagcat::eval(const cpv_MX &arg, const pv_MX &res) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Diagcat::evalAdj(const std::vector< pv_MX > &adjSeed, const std::vector< pv_MX > &adjSens) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Diagcat::evalFwd(const std::vector< cpv_MX > &fwdSeed, const std::vector< pv_MX > &fwdSens) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Diagcat::getOp() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Diagcat::printPart(std::ostream &stream, int part) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Diagsplit::clone() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Diagsplit::eval(const cpv_MX &input, const pv_MX &output) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Diagsplit::evalAdj(const std::vector< pv_MX > &adjSeed, const std::vector< pv_MX > &adjSens) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Diagsplit::evalFwd(const std::vector< cpv_MX > &fwdSeed, const std::vector< pv_MX > &fwdSens) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Diagsplit::getDiagcat(const std::vector< MX > &x) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Diagsplit::getOp() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Diagsplit::printPart(std::ostream &stream, int part) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DleToDple::clone() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DleToDple::create(const DleStructure &st) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DleToDple::evaluate() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DleToDple::getDerForward(int nfwd) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DleToDple::getDerReverse(int nadj) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DleToDple::hasDerForward() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DleToDple::hasDerReverse() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DleToDple::init() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DleToDple::printStats(std::ostream &stream) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DleToLrDle::clone() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DleToLrDle::create(const DleStructure &st) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DleToLrDle::evaluate() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DleToLrDle::getDerForward(int nfwd) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DleToLrDle::getDerReverse(int nadj) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DleToLrDle::hasDerForward() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DleToLrDle::hasDerReverse() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DleToLrDle::init() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DleToLrDle::printStats(std::ostream &stream) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DpleToLrDple::clone() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DpleToLrDple::create(const DpleStructure &st) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DpleToLrDple::evaluate() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DpleToLrDple::getDerForward(int nfwd) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DpleToLrDple::getDerReverse(int nadj) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DpleToLrDple::hasDerForward() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DpleToLrDple::hasDerReverse() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DpleToLrDple::init() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DpleToLrDple::printStats(std::ostream &stream) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DsdpInterface::clone() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DsdpInterface::evaluate() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DsdpInterface::init() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FixedStepIntegrator::calculateInitialConditions() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FixedStepIntegrator::calculateInitialConditionsB() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FixedStepIntegrator::clone() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FixedStepIntegrator::create(const Function &f, const Function &g) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FixedStepIntegrator::getExplicit() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FixedStepIntegrator::getExplicitB() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FixedStepIntegrator::init() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FixedStepIntegrator::integrate(double t_out) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FixedStepIntegrator::integrateB(double t_out) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FixedStepIntegrator::reset() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FixedStepIntegrator::resetB() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FixedStepIntegrator::setupFG() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::callDerivative(const DMatrixVector &arg, DMatrixVector &output_res, const DMatrixVectorVector &fseed, DMatrixVectorVector &output_fsens, const DMatrixVectorVector &aseed, DMatrixVectorVector &output_asens, bool always_inline=false, bool never_inline=false) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::callDerivative(const MXVector &arg, MXVector &output_res, const MXVectorVector &fseed, MXVectorVector &output_fsens, const MXVectorVector &aseed, MXVectorVector &output_asens, bool always_inline=false, bool never_inline=false) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::callDerivative(const SXVector &arg, SXVector &output_res, const SXVectorVector &fseed, SXVectorVector &output_fsens, const SXVectorVector &aseed, SXVectorVector &output_asens, bool always_inline=false, bool never_inline=false) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::checkAdjSeed(const std::vector< std::vector< M > > &aseed) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::checkArg(const std::vector< M > &arg) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::checkFwdSeed(const std::vector< std::vector< M > > &fseed) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::checkInputs() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::checkRes(const std::vector< M > &res) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::inputScheme() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::inputScheme() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::input_struct() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::input_struct() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::matchingAdjSeed(const std::vector< std::vector< M > > &aseed) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::matchingArg(const std::vector< M > &arg) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::matchingFwdSeed(const std::vector< std::vector< M > > &fseed) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::matchingRes(const std::vector< M > &arg) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::outputScheme() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::outputScheme() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::output_struct() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::output_struct() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::replaceAdjSeed(const std::vector< std::vector< M > > &aseed) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::replaceArg(const std::vector< M > &arg) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::replaceFwdSeed(const std::vector< std::vector< M > > &fseed) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::replaceRes(const std::vector< M > &res) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::spCanEvaluate(bool fwd) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::spEvaluate(bool fwd) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::spInit(bool fwd) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::adWeight() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::adWeightSp() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::adjViaJac(int nadj) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::call(const DMatrixVector &arg, DMatrixVector &res, bool always_inline, bool never_inline) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::call(const MXVector &arg, MXVector &res, bool always_inline, bool never_inline) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::call(const SXVector &arg, SXVector &res, bool always_inline, bool never_inline) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::callForward(const std::vector< DMatrix > &arg, const std::vector< DMatrix > &res, const std::vector< std::vector< DMatrix > > &fseed, std::vector< std::vector< DMatrix > > &fsens, bool always_inline, bool never_inline) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::callForward(const std::vector< MX > &arg, const std::vector< MX > &res, const std::vector< std::vector< MX > > &fseed, std::vector< std::vector< MX > > &fsens, bool always_inline, bool never_inline) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::callForward(const std::vector< SX > &arg, const std::vector< SX > &res, const std::vector< std::vector< SX > > &fseed, std::vector< std::vector< SX > > &fsens, bool always_inline, bool never_inline) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::callReverse(const std::vector< DMatrix > &arg, const std::vector< DMatrix > &res, const std::vector< std::vector< DMatrix > > &aseed, std::vector< std::vector< DMatrix > > &asens, bool always_inline, bool never_inline) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::callReverse(const std::vector< MX > &arg, const std::vector< MX > &res, const std::vector< std::vector< MX > > &aseed, std::vector< std::vector< MX > > &asens, bool always_inline, bool never_inline) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::callReverse(const std::vector< SX > &arg, const std::vector< SX > &res, const std::vector< std::vector< SX > > &aseed, std::vector< std::vector< SX > > &asens, bool always_inline, bool never_inline) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::checkInputs() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::createCall(const std::vector< MX > &arg) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::derForward(int nfwd) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::derReverse(int nadj) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::dynamicCompilation(Function f, std::string fname, std::string fdescr, std::string compiler) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::evalMX(const std::vector< MX > &arg, std::vector< MX > &res) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::evalSX(const std::vector< SX > &arg, std::vector< SX > &res) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::fullJacobian() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::fwdViaJac(int nfwd) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::generate(std::ostream &stream, const std::vector< int > &arg, const std::vector< int > &res, CodeGenerator &gen) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::generateBody(std::ostream &stream, const std::string &type, CodeGenerator &gen) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::generateCode(std::ostream &cfile, bool generate_main) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::generateDeclarations(std::ostream &stream, const std::string &type, CodeGenerator &gen) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::generateFunction(std::ostream &stream, const std::string &fname, const std::string &type, CodeGenerator &gen) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::generateIO(CodeGenerator &gen) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::getAdaptorSolverName() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::getDerForward(int nfwd) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::getDerReverse(int nadj) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::getFullJacobian() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::getGradient(int iind, int oind) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::getHessian(int iind, int oind) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::getJacSparsity(int iind, int oind, bool symmetric) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::getJacSparsityHierarchical(int iind, int oind) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::getJacSparsityHierarchicalSymm(int iind, int oind) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::getJacSparsityPlain(int iind, int oind) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::getJacobian(int iind, int oind, bool compact, bool symmetric) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::getNumInputElements() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::getNumInputNonzeros() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::getNumOutputElements() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::getNumOutputNonzeros() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::getNumericJacobian(int iind, int oind, bool compact, bool symmetric) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::getPartition(int iind, int oind, Sparsity &D1, Sparsity &D2, bool compact, bool symmetric) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::getSanitizedName() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::getStat(const std::string &name) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::getStats() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::getTangent(int iind, int oind) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::gradient(int iind, int oind) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::hasDerForward() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::hasDerReverse() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::hasDerivative() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::hessian(int iind, int oind) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::inputNoCheck(int iind=0) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::inputNoCheck(int iind=0) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::inputScheme() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::inputScheme() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::input_struct() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::input_struct() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::jacSparsity(int iind, int oind, bool compact, bool symmetric) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::jacobian(int iind, int oind, bool compact, bool symmetric) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::log(const std::string &fcn, const std::string &msg) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::log(const std::string &msg) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::monitored(const std::string &mod) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::nTmp(size_t &ni, size_t &nr) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::outputNoCheck(int oind=0) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::outputNoCheck(int oind=0) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::outputScheme() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::outputScheme() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::output_struct() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::output_struct() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::print(std::ostream &stream) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::repr(std::ostream &stream) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::setDerForward(const Function &fcn, int nfwd) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::setDerReverse(const Function &fcn, int nadj) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::setJacSparsity(const Sparsity &sp, int iind, int oind, bool compact) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::setJacobian(const Function &jac, int iind, int oind, bool compact) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::spCanEvaluate(bool fwd) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::spEvaluate(bool fwd) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::spEvaluateViaJacSparsity(bool fwd) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::spInit(bool fwd) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::symbolicAdjSeed(int nadj, const std::vector< MatType > &v) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::symbolicFwdSeed(int nfwd, const std::vector< MatType > &v) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::symbolicInput() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::symbolicInputSX() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::symbolicOutput() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::symbolicOutput(const std::vector< MX > &arg) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::tangent(int iind, int oind) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::verbose() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::wrapMXFunction() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericExpression< SXElement  >::zz_ge(const SXElement &y) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericExpression< SXElement  >::zz_gt(const SXElement &y) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< MX  >::colind() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< MX  >::row() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< MatType >::colind() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< MatType >::row() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< MatType >::sparsityRef() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< Matrix< DataType >  >::colind() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< Matrix< DataType >  >::row() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericType::is_a() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericType::toDictionary() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericType::toDictionary() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericType::toDoubleVector() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericType::toDoubleVector() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericType::toIntVector() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericType::toIntVector() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericType::toIntVectorVector() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericType::toIntVectorVector() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GetNonzeros::eval(const cpv_MX &input, const pv_MX &output) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GetNonzeros::evalAdj(const std::vector< pv_MX > &adjSeed, const std::vector< pv_MX > &adjSens) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GetNonzeros::evalFwd(const std::vector< cpv_MX > &fwdSeed, const std::vector< pv_MX > &fwdSens) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GetNonzeros::getAll() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GetNonzeros::getGetNonzeros(const Sparsity &sp, const std::vector< int > &nz) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GetNonzeros::getOp() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GetNonzeros::mapping() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GetNonzerosSlice2::clone() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GetNonzerosSlice2::generate(std::ostream &stream, const std::vector< int > &arg, const std::vector< int > &res, CodeGenerator &gen) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GetNonzerosSlice2::getAll() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GetNonzerosSlice2::printPart(std::ostream &stream, int part) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GetNonzerosSlice::clone() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GetNonzerosSlice::generate(std::ostream &stream, const std::vector< int > &arg, const std::vector< int > &res, CodeGenerator &gen) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GetNonzerosSlice::getAll() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GetNonzerosSlice::isIdentity() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GetNonzerosSlice::printPart(std::ostream &stream, int part) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GetNonzerosSlice::simplifyMe(MX &ex) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GetNonzerosVector::clone() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GetNonzerosVector::generate(std::ostream &stream, const std::vector< int > &arg, const std::vector< int > &res, CodeGenerator &gen) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GetNonzerosVector::getAll() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GetNonzerosVector::printPart(std::ostream &stream, int part) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Horzcat::clone() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Horzcat::eval(const cpv_MX &arg, const pv_MX &res) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Horzcat::evalAdj(const std::vector< pv_MX > &adjSeed, const std::vector< pv_MX > &adjSens) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Horzcat::evalFwd(const std::vector< cpv_MX > &fwdSeed, const std::vector< pv_MX > &fwdSens) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Horzcat::getOp() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Horzcat::printPart(std::ostream &stream, int part) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Horzsplit::clone() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Horzsplit::eval(const cpv_MX &input, const pv_MX &output) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Horzsplit::evalAdj(const std::vector< pv_MX > &adjSeed, const std::vector< pv_MX > &adjSens) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Horzsplit::evalFwd(const std::vector< cpv_MX > &fwdSeed, const std::vector< pv_MX > &fwdSens) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Horzsplit::getHorzcat(const std::vector< MX > &x) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Horzsplit::getOp() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Horzsplit::printPart(std::ostream &stream, int part) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IOInterface< Derived >::getInput(T val, const std::string &iname) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IOInterface< Derived >::getInput(T val, int iind=0) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IOInterface< Derived >::getOutput(T val, const std::string &oname) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IOInterface< Derived >::getOutput(T val, int oind=0) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IOInterface< Derived >::inputS(int i) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IOInterface< Derived >::inputS(int i) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IOInterface< Derived >::inputSchemeEntry(const std::string &name) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IOInterface< Derived >::outputS(int i) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IOInterface< Derived >::outputS(int i) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IOInterface< Derived >::outputSchemeEntry(const std::string &name) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IOInterface< Derived >::schemeEntry(const casadi::IOScheme &scheme, const std::string &name, bool input) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IOInterface< Function  >::getInput(T val, const std::string &iname) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IOInterface< Function  >::getInput(T val, int iind=0) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IOInterface< Function  >::getOutput(T val, const std::string &oname) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IOInterface< Function  >::getOutput(T val, int oind=0) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IOInterface< Function  >::inputS(int i) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IOInterface< Function  >::inputS(int i) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IOInterface< Function  >::inputSchemeEntry(const std::string &name) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IOInterface< Function  >::outputS(int i) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IOInterface< Function  >::outputS(int i) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IOInterface< Function  >::outputSchemeEntry(const std::string &name) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IOInterface< Function  >::schemeEntry(const casadi::IOScheme &scheme, const std::string &name, bool input) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IOInterface< FunctionInternal  >::getInput(T val, const std::string &iname) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IOInterface< FunctionInternal  >::getInput(T val, int iind=0) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IOInterface< FunctionInternal  >::getInput(const std::string &iname) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IOInterface< FunctionInternal  >::getInput(int iind=0) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IOInterface< FunctionInternal  >::getInputScheme() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IOInterface< FunctionInternal  >::getNumInputs() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IOInterface< FunctionInternal  >::getNumOutputs() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IOInterface< FunctionInternal  >::getOutput(T val, const std::string &oname) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IOInterface< FunctionInternal  >::getOutput(T val, int oind=0) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IOInterface< FunctionInternal  >::getOutput(const std::string &oname) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IOInterface< FunctionInternal  >::getOutput(int oind=0) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IOInterface< FunctionInternal  >::getOutputScheme() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IOInterface< FunctionInternal  >::input(const std::string &iname) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IOInterface< FunctionInternal  >::input(const std::string &iname) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IOInterface< FunctionInternal  >::input(int iind=0) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IOInterface< FunctionInternal  >::input(int iind=0) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IOInterface< FunctionInternal  >::inputS(int i) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IOInterface< FunctionInternal  >::inputS(int i) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IOInterface< FunctionInternal  >::inputSchemeEntry(const std::string &name) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IOInterface< FunctionInternal  >::output(const std::string &oname) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IOInterface< FunctionInternal  >::output(const std::string &oname) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IOInterface< FunctionInternal  >::output(int oind=0) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IOInterface< FunctionInternal  >::output(int oind=0) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IOInterface< FunctionInternal  >::outputS(int i) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IOInterface< FunctionInternal  >::outputS(int i) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IOInterface< FunctionInternal  >::outputSchemeEntry(const std::string &name) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IOInterface< FunctionInternal  >::schemeEntry(const casadi::IOScheme &scheme, const std::string &name, bool input) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IOInterface< FunctionInternal  >::setInput(T val, const std::string &iname) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IOInterface< FunctionInternal  >::setInput(T val, int iind=0) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IOInterface< FunctionInternal  >::setInputScheme(const casadi::IOScheme &scheme) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IOInterface< FunctionInternal  >::setNumInputs(int num_in) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IOInterface< FunctionInternal  >::setNumOutputs(int num_out) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IOInterface< FunctionInternal  >::setOutput(T val, const std::string &oname) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IOInterface< FunctionInternal  >::setOutput(T val, int oind=0) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IOInterface< FunctionInternal  >::setOutputScheme(const casadi::IOScheme &scheme) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IOScheme::print(std::ostream &stream=std::cout) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IOScheme::repr(std::ostream &stream) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IOSchemeVector< M  >::print(std::ostream &stream=CASADI_COUT, bool trailing_newline=true) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IOSchemeVector< M  >::repr(std::ostream &stream=CASADI_COUT, bool trailing_newline=true) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IOSchemeVector< M  >::vector() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IdasInterface::clone() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IdasInterface::correctInitialConditions() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IdasInterface::create(const Function &f, const Function &g) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IdasInterface::freeIDAS() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IdasInterface::getJac() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IdasInterface::getJacB() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IdasInterface::getJacGen() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IdasInterface::getJacGenB() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IdasInterface::init() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IdasInterface::initAdj() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IdasInterface::initBandedLinearSolver() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IdasInterface::initBandedLinearSolverB() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IdasInterface::initDenseLinearSolver() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IdasInterface::initDenseLinearSolverB() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IdasInterface::initIterativeLinearSolver() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IdasInterface::initIterativeLinearSolverB() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IdasInterface::initTaping() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IdasInterface::initUserDefinedLinearSolver() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IdasInterface::initUserDefinedLinearSolverB() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IdasInterface::integrate(double t_out) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IdasInterface::integrateB(double t_out) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IdasInterface::printStats(std::ostream &stream) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IdasInterface::reset() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IdasInterface::resetB() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IdasInterface::setStopTime(double tf) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ImplicitFixedStepIntegrator::clone() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ImplicitFixedStepIntegrator::create(const Function &f, const Function &g) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ImplicitFixedStepIntegrator::getExplicit() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ImplicitFixedStepIntegrator::getExplicitB() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ImplicitFixedStepIntegrator::init() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ImplicitFunctionInternal::callForward(const std::vector< MX > &arg, const std::vector< MX > &res, const std::vector< std::vector< MX > > &fseed, std::vector< std::vector< MX > > &fsens, bool always_inline, bool never_inline) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ImplicitFunctionInternal::callReverse(const std::vector< MX > &arg, const std::vector< MX > &res, const std::vector< std::vector< MX > > &aseed, std::vector< std::vector< MX > > &asens, bool always_inline, bool never_inline) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ImplicitFunctionInternal::evaluate() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ImplicitFunctionInternal::getDerForward(int nfwd) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ImplicitFunctionInternal::getDerReverse(int nadj) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ImplicitFunctionInternal::hasDerForward() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ImplicitFunctionInternal::hasDerReverse() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ImplicitFunctionInternal::spCanEvaluate(bool fwd) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ImplicitFunctionInternal::spEvaluate(bool fwd) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::InfSX::getValue() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::InfSX::isInf() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::InnerProd::clone() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::InnerProd::eval(const cpv_MX &input, const pv_MX &output) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::InnerProd::evalAdj(const std::vector< pv_MX > &adjSeed, const std::vector< pv_MX > &adjSens) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::InnerProd::evalFwd(const std::vector< cpv_MX > &fwdSeed, const std::vector< pv_MX > &fwdSens) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::InnerProd::generate(std::ostream &stream, const std::vector< int > &arg, const std::vector< int > &res, CodeGenerator &gen) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::InnerProd::getOp() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::InnerProd::printPart(std::ostream &stream, int part) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IntegerSX::getIntValue() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IntegerSX::getValue() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IntegerSX::isInteger() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IntegratorInternal::clone() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IntegratorInternal::create(const Function &f, const Function &g) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IntegratorInternal::evaluate() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IntegratorInternal::getAugOffset(int nfwd, int nadj) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IntegratorInternal::getAugmented(int nfwd, int nadj, AugOffset &offset) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IntegratorInternal::getDerForward(int nfwd) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IntegratorInternal::getDerReverse(int nadj) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IntegratorInternal::hasDerForward() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IntegratorInternal::hasDerReverse() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IntegratorInternal::integrate(double t_out) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IntegratorInternal::integrateB(double t_out) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IntegratorInternal::p() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IntegratorInternal::printStats(std::ostream &stream) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IntegratorInternal::qf() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IntegratorInternal::resetB() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IntegratorInternal::rp() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IntegratorInternal::rqf() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IntegratorInternal::rx0() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IntegratorInternal::rxf() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IntegratorInternal::rz0() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IntegratorInternal::rzf() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IntegratorInternal::setDerivativeOptions(Integrator &integrator, const AugOffset &offset) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IntegratorInternal::setStopTime(double tf) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IntegratorInternal::spCanEvaluate(bool fwd) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IntegratorInternal::spEvaluate(bool fwd) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IntegratorInternal::spJacF() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IntegratorInternal::spJacG() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IntegratorInternal::x0() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IntegratorInternal::xf() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IntegratorInternal::z0() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IntegratorInternal::zf() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Inverse::clone() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Inverse::eval(const cpv_MX &input, const pv_MX &output) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Inverse::evalAdj(const std::vector< pv_MX > &adjSeed, const std::vector< pv_MX > &adjSens) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Inverse::evalFwd(const std::vector< cpv_MX > &fwdSeed, const std::vector< pv_MX > &fwdSens) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Inverse::getOp() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Inverse::printPart(std::ostream &stream, int part) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IpoptInterface::clone() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IpoptInterface::evaluate() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IpoptInterface::finalize_metadata(int n, const std::map< std::string, std::vector< std::string > > &var_string_md, const std::map< std::string, std::vector< int > > &var_integer_md, const std::map< std::string, std::vector< double > > &var_numeric_md, int m, const std::map< std::string, std::vector< std::string > > &con_string_md, const std::map< std::string, std::vector< int > > &con_integer_md, const std::map< std::string, std::vector< double > > &con_numeric_md) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IpoptInterface::freeIpopt() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IpoptInterface::getReducedHessian() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IpoptInterface::get_nlp_info(int &n, int &m, int &nnz_jac_g, int &nnz_h_lag) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IpoptInterface::get_number_of_nonlinear_variables() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IpoptInterface::get_var_con_metadata(int n, std::map< std::string, std::vector< std::string > > &var_string_md, std::map< std::string, std::vector< int > > &var_integer_md, std::map< std::string, std::vector< double > > &var_numeric_md, int m, std::map< std::string, std::vector< std::string > > &con_string_md, std::map< std::string, std::vector< int > > &con_integer_md, std::map< std::string, std::vector< double > > &con_numeric_md) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IpoptInterface::init() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IpoptInterface::setQPOptions() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IpoptUserClass::finalize_metadata(Index n, const StringMetaDataMapType &var_string_md, const IntegerMetaDataMapType &var_integer_md, const NumericMetaDataMapType &var_numeric_md, Index m, const StringMetaDataMapType &con_string_md, const IntegerMetaDataMapType &con_integer_md, const NumericMetaDataMapType &con_numeric_md) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IpoptUserClass::get_nlp_info(Index &n, Index &m, Index &nnz_jac_g, Index &nnz_h_lag, IndexStyleEnum &index_style) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IpoptUserClass::get_number_of_nonlinear_variables() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IpoptUserClass::get_var_con_metadata(Index n, StringMetaDataMapType &var_string_md, IntegerMetaDataMapType &var_integer_md, NumericMetaDataMapType &var_numeric_md, Index m, StringMetaDataMapType &con_string_md, IntegerMetaDataMapType &con_integer_md, NumericMetaDataMapType &con_numeric_md) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::KinsolInterface::bjac(long N, long mupper, long mlower, N_Vector u, N_Vector fu, DlsMat J, N_Vector tmp1, N_Vector tmp2) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::KinsolInterface::clone() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::KinsolInterface::create(const Function &f) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::KinsolInterface::djac(long N, N_Vector u, N_Vector fu, DlsMat J, N_Vector tmp1, N_Vector tmp2) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::KinsolInterface::func(N_Vector u, N_Vector fval) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::KinsolInterface::init() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::KinsolInterface::kinsol_error(const std::string &module, int flag, bool fatal=true) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::KinsolInterface::lsetup(KINMem kin_mem) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::KinsolInterface::psetup(N_Vector u, N_Vector uscale, N_Vector fval, N_Vector fscale, N_Vector tmp1, N_Vector tmp2) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::KinsolInterface::psolve(N_Vector u, N_Vector uscale, N_Vector fval, N_Vector fscale, N_Vector v, N_Vector tmp) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::KinsolInterface::solveNonLinear() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::KnitroInterface::clone() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::KnitroInterface::evaluate() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::KnitroInterface::init() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::LapackLuDense::clone() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::LapackLuDense::init() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::LapackLuDense::prepare() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::LapackQrDense::clone() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::LapackQrDense::init() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::LapackQrDense::prepare() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::LinearSolver::spSolve(DMatrix &X, const DMatrix &B, bool transpose=false) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::LinearSolverInternal::callForwardLinsol(const std::vector< MX > &arg, const std::vector< MX > &res, const std::vector< std::vector< MX > > &fseed, std::vector< std::vector< MX > > &fsens, bool tr) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::LinearSolverInternal::callReverseLinsol(const std::vector< MX > &arg, const std::vector< MX > &res, const std::vector< std::vector< MX > > &aseed, std::vector< std::vector< MX > > &asens, bool tr) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::LinearSolverInternal::colind() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::LinearSolverInternal::evaluate() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::LinearSolverInternal::getFactorization(bool transpose) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::LinearSolverInternal::getFactorizationSparsity(bool transpose) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::LinearSolverInternal::ncol() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::LinearSolverInternal::nnz() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::LinearSolverInternal::nrow() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::LinearSolverInternal::row() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::LinearSolverInternal::solve(bool transpose) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::LinearSolverInternal::solve(const MX &A, const MX &B, bool transpose) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::LinearSolverInternal::spSolve(DMatrix &X, const DMatrix &B, bool transpose) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::LpSolverInternal::checkInputs() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::LpSolverInternal::solve() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::LpToQp::clone() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::LpToQp::evaluate() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::LpToQp::init() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::LrDleInternal::create(const LrDleStructure &st) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::LrDleToDle::clone() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::LrDleToDle::create(const LrDleStructure &st) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::LrDleToDle::evaluate() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::LrDleToDle::getDerForward(int nfwd) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::LrDleToDle::getDerReverse(int nadj) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::LrDleToDle::hasDerForward() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::LrDleToDle::hasDerReverse() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::LrDleToDle::init() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::LrDleToDle::printStats(std::ostream &stream) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::LrDpleToDple::clone() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::LrDpleToDple::create(const LrDpleStructure &st) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::LrDpleToDple::evaluate() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::LrDpleToDple::getDerForward(int nfwd) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::LrDpleToDple::getDerReverse(int nadj) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::LrDpleToDple::hasDerForward() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::LrDpleToDple::hasDerReverse() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::LrDpleToDple::init() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::LrDpleToDple::printStats(std::ostream &stream) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::LrDpleToLrDle::clone() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::LrDpleToLrDle::create(const LrDleStructure &st, const std::vector< int > &Hs) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::LrDpleToLrDle::evaluate() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::LrDpleToLrDle::getDerForward(int nfwd) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::LrDpleToLrDle::getDerReverse(int nadj) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::LrDpleToLrDle::hasDerForward() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::LrDpleToLrDle::hasDerReverse() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::LrDpleToLrDle::init() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::LrDpleToLrDle::printStats(std::ostream &stream) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::at(int k) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::at(int k) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::getTemp() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::setTemp(int t) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::sparsity() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::sparsityRef() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MXFunction::algorithm() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MXFunction::generateLiftingFunctions(MXFunction &vdef_fcn, MXFunction &vinit_fcn) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MXNode::addDependency(const MX &dep) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MXNode::assign(const MX &d, const std::vector< int > &inz, bool add=false) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MXNode::assign(const MX &d, const std::vector< int > &inz, const std::vector< int > &onz, bool add=false) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MXNode::clone() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MXNode::dep(int ind=0) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MXNode::dep(int ind=0) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MXNode::eval(const cpv_MX &arg, const pv_MX &res) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MXNode::evalAdj(const std::vector< pv_MX > &aseed, const std::vector< pv_MX > &fsens) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MXNode::evalFwd(const std::vector< cpv_MX > &fseed, const std::vector< pv_MX > &fsens) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MXNode::generate(std::ostream &stream, const std::vector< int > &arg, const std::vector< int > &res, CodeGenerator &gen) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MXNode::getAddNonzeros(const MX &y, const std::vector< int > &nz) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MXNode::getAssertion(const MX &y, const std::string &fail_message="") const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MXNode::getAssign(const MX &y, const Slice &i, const Slice &j) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MXNode::getBinary(int op, const MX &y, bool scX, bool scY) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MXNode::getBinarySwitch(int op, const MX &y) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MXNode::getDeterminant() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MXNode::getDiagcat(const std::vector< MX > &x) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MXNode::getDiagsplit(const std::vector< int > &offset1, const std::vector< int > &offset2) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MXNode::getFunction() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MXNode::getFunction() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MXNode::getFunctionInput() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MXNode::getFunctionOutput() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MXNode::getGetNonzeros(const Sparsity &sp, const std::vector< int > &nz) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MXNode::getHorzcat(const std::vector< MX > &x) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MXNode::getHorzsplit(const std::vector< int > &output_offset) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MXNode::getInnerProd(const MX &y) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MXNode::getInverse() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MXNode::getMatrixValue() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MXNode::getMultiplication(const MX &y, const MX &z) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MXNode::getName() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MXNode::getNorm1() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MXNode::getNorm2() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MXNode::getNormF() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MXNode::getNormInf() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MXNode::getOp() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MXNode::getOutput(int oind) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MXNode::getRef(const Slice &i, const Slice &j) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MXNode::getReshape(const Sparsity &sp) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MXNode::getSetNonzeros(const MX &y, const std::vector< int > &nz) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MXNode::getSetSparse(const Sparsity &sp) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MXNode::getSolve(const MX &r, bool tr, const LinearSolver &linear_solver) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MXNode::getTranspose() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MXNode::getUnary(int op) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MXNode::getValue() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MXNode::getVertcat(const std::vector< MX > &x) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MXNode::getVertsplit(const std::vector< int > &output_offset) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MXNode::hasDep() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MXNode::isBinaryOp() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MXNode::isIdentity() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MXNode::isMultipleOutput() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MXNode::isNonLinear() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MXNode::isOne() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MXNode::isOutputNode() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MXNode::isUnaryOp() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MXNode::isValue(double val) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MXNode::isZero() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MXNode::mapping() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MXNode::nTmp(size_t &ni, size_t &nr) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MXNode::ndep() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MXNode::nnz() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MXNode::nout() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MXNode::numInplace() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MXNode::numel() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MXNode::print(std::ostream &stream) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MXNode::print(std::ostream &stream, long &remaining_calls) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MXNode::printPart(std::ostream &stream, int part) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MXNode::repr(std::ostream &stream) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MXNode::setDependencies(const MX &dep) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MXNode::setDependencies(const MX &dep1, const MX &dep2) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MXNode::setDependencies(const MX &dep1, const MX &dep2, const MX &dep3) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MXNode::setDependencies(const std::vector< MX > &dep) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MXNode::setSparsity(const Sparsity &sparsity) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MXNode::shape() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MXNode::simplifyMe(MX &ex) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MXNode::size1() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MXNode::size2() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MXNode::sparsity() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MXNode::sparsity(int oind) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::T() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::append(const Matrix< DataType > &y) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::appendColumns(const Matrix< DataType > &y) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::binary(int op, const Matrix< DataType > &x, const Matrix< DataType > &y) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::borBV(const Matrix< DataType > &val) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::clear() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::data() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::data() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::elem(int rr, int cc=0) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::elem(int rr, int cc=0) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::enlarge(int nrow, int ncol, const std::vector< int > &rr, const std::vector< int > &cc, bool ind1=false) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::erase(const std::vector< int > &rr, bool ind1=false) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::erase(const std::vector< int > &rr, const std::vector< int > &cc, bool ind1=false) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::get(Matrix< DataType > &output_m, bool ind1, const Matrix< int > &rr) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::get(Matrix< DataType > &output_m, bool ind1, const Matrix< int > &rr, const Matrix< int > &cc) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::get(Matrix< DataType > &output_m, bool ind1, const Matrix< int > &rr, const Slice &cc) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::get(Matrix< DataType > &output_m, bool ind1, const Slice &rr) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::get(Matrix< DataType > &output_m, bool ind1, const Slice &rr, const Matrix< int > &cc) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::get(Matrix< DataType > &output_m, bool ind1, const Slice &rr, const Slice &cc) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::get(Matrix< DataType > &output_m, bool ind1, const Sparsity &sp) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::get(Matrix< DataType > &val) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::get(double &val) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::get(std::vector< double > &output_m, bool tr=false) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::getDep(int ch=0) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::getElementHash() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::getIntValue() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::getNZ(Matrix< DataType > &output_m, bool ind1, const Matrix< int > &k) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::getNZ(Matrix< DataType > &output_m, bool ind1, const Slice &k) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::getNZ(double &val) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::getNZ(std::vector< double > &output_m) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::getName() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::getNdeps() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::getSym(std::vector< double > &output_m) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::getValue() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::getValue(int k) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::hasNonStructuralZeros() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::inf(const Sparsity &sp) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::inf(const std::pair< int, int > &rc) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::inf(int nrow=1, int ncol=1) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::isCommutative() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::isConstant() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::isIdentity() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::isInteger() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::isLeaf() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::isMinusOne() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::isOne() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::isRegular() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::isSlice(bool ind1=false) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::isSmooth() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::isSymbolic() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::isSymbolicSparse() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::isZero() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::makeDense(const DataType &val=0) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::makeSparse(double tol=0) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::matrix_matrix(int op, const Matrix< DataType > &x, const Matrix< DataType > &y) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::matrix_scalar(int op, const Matrix< DataType > &x, const Matrix< DataType > &y) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::nan(const Sparsity &sp) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::nan(const std::pair< int, int > &rc) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::nan(int nrow=1, int ncol=1) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::nonzeros() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::nonzeros_int() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::print(std::ostream &stream=CASADI_COUT, bool trailing_newline=true) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::printDense(std::ostream &stream=CASADI_COUT, bool trailing_newline=true) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::printScalar(std::ostream &stream=CASADI_COUT, bool trailing_newline=true) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::printSparse(std::ostream &stream=CASADI_COUT, bool trailing_newline=true) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::printVector(std::ostream &stream=CASADI_COUT, bool trailing_newline=true) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::printme(const Matrix< DataType > &y) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::remove(const std::vector< int > &rr, const std::vector< int > &cc) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::repr(std::ostream &stream=CASADI_COUT, bool trailing_newline=true) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::reserve(int nnz) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::reserve(int nnz, int ncol) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::resize(int nrow, int ncol) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::sanityCheck(bool complete=false) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::scalar_matrix(int op, const Matrix< DataType > &x, const Matrix< DataType > &y) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::set(const Matrix< DataType > &m, bool ind1, const Matrix< int > &rr) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::set(const Matrix< DataType > &m, bool ind1, const Matrix< int > &rr, const Matrix< int > &cc) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::set(const Matrix< DataType > &m, bool ind1, const Matrix< int > &rr, const Slice &cc) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::set(const Matrix< DataType > &m, bool ind1, const Slice &rr) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::set(const Matrix< DataType > &m, bool ind1, const Slice &rr, const Matrix< int > &cc) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::set(const Matrix< DataType > &m, bool ind1, const Slice &rr, const Slice &cc) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::set(const Matrix< DataType > &m, bool ind1, const Sparsity &sp) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::set(const Matrix< DataType > &val) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::set(const std::vector< double > &val, bool tr=false) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::set(double val) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::setAll(const DataType &val) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::setBV(const Matrix< DataType > &val) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::setNZ(const Matrix< DataType > &m, bool ind1, const Matrix< int > &k) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::setNZ(const Matrix< DataType > &m, bool ind1, const Slice &k) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::setNZ(const std::vector< double > &val) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::setNZ(double val) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::setSparse(const Sparsity &sp, bool intersect=false) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::setSym(const std::vector< double > &val) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::setValue(double m) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::setValue(double m, int k) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::setZero() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::setZeroBV() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::sparsityRef() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::toScalar() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::toSlice(bool ind1=false) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::triplet(const std::vector< int > &row, const std::vector< int > &col, const Matrix< DataType > &d) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::triplet(const std::vector< int > &row, const std::vector< int > &col, const Matrix< DataType > &d, const std::pair< int, int > &rc) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::triplet(const std::vector< int > &row, const std::vector< int > &col, const Matrix< DataType > &d, int nrow, int ncol) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::unary(int op, const Matrix< DataType > &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::zz_abs() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::zz_acos() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::zz_acosh() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::zz_adj() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::zz_all() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::zz_and(const Matrix< DataType > &y) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::zz_any() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::zz_asin() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::zz_asinh() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::zz_atan() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::zz_atan2(const Matrix< DataType > &y) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::zz_atanh() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::zz_blockcat(const std::vector< std::vector< Matrix< DataType > > > &v) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::zz_ceil() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::zz_cofactor(int i, int j) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::zz_cos() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::zz_cosh() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::zz_countNodes() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::zz_dependsOn(const Matrix< DataType > &arg) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::zz_det() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::zz_diag() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::zz_diagcat(const std::vector< Matrix< DataType > > &A) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::zz_diagsplit(const std::vector< int > &offset1, const std::vector< int > &offset2) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::zz_eig_symbolic() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::zz_eq(const Matrix< DataType > &y) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::zz_erf() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::zz_erfinv() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::zz_exp() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::zz_expand(Matrix< DataType > &weights, Matrix< DataType > &terms) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::zz_extractShared(std::vector< Matrix< DataType > > &ex, std::vector< Matrix< DataType > > &v, std::vector< Matrix< DataType > > &vdef, const std::string &v_prefix="v_", const std::string &v_suffix="") {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::zz_floor() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::zz_gauss_quadrature(const Matrix< DataType > &x, const Matrix< DataType > &a, const Matrix< DataType > &b, int order, const Matrix< DataType > &w) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::zz_getMinor(int i, int j) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::zz_getOperatorRepresentation(const std::vector< std::string > &args) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::zz_getSymbols() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::zz_getSymbols(const std::vector< Matrix< DataType > > &e) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::zz_gradient(const Matrix< DataType > &arg) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::zz_heaviside() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::zz_hessian(const Matrix< DataType > &arg) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::zz_hessian(const Matrix< DataType > &arg, Matrix< DataType > &H, Matrix< DataType > &g) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::zz_horzcat(const std::vector< Matrix< DataType > > &v) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::zz_horzsplit(const std::vector< int > &offset) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::zz_if_else(const Matrix< DataType > &if_true, const Matrix< DataType > &if_false) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::zz_if_else_zero(const Matrix< DataType > &y) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::zz_inner_prod(const Matrix< DataType > &y) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::zz_inv() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::zz_isEqual(const Matrix< DataType > &ex2, int depth=0) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::zz_jacobian(const Matrix< DataType > &arg) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::zz_jacobianTimesVector(const Matrix< DataType > &arg, const Matrix< DataType > &v, bool transpose_jacobian=false) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::zz_kron(const Matrix< DataType > &b) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::zz_le(const Matrix< DataType > &y) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::zz_log() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::zz_log10() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::zz_lt(const Matrix< DataType > &y) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::zz_max(const Matrix< DataType > &y) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::zz_min(const Matrix< DataType > &y) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::zz_minus(const Matrix< DataType > &y) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::zz_mod(const Matrix< DataType > &y) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::zz_mpower(const Matrix< DataType > &y) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::zz_mtaylor(const Matrix< DataType > &x, const Matrix< DataType > &a, int order, const std::vector< int > &order_contributions) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::zz_mtaylor(const Matrix< DataType > &x, const Matrix< DataType > &a, int order=1) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::zz_mtimes(const Matrix< DataType > &y) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::zz_mtimes(const Matrix< DataType > &y, const Matrix< DataType > &z) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::zz_ne(const Matrix< DataType > &y) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::zz_norm_1() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::zz_norm_2() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::zz_norm_F() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::zz_norm_inf() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::zz_norm_inf_mul(const Matrix< DataType > &y) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::zz_not() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::zz_nullspace() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::zz_or(const Matrix< DataType > &y) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::zz_outer_prod(const Matrix< DataType > &y) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::zz_pinv() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::zz_plus(const Matrix< DataType > &y) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::zz_poly_coeff(const Matrix< DataType > &x) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::zz_poly_roots() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::zz_polyval(const Matrix< DataType > &x) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::zz_power(const Matrix< DataType > &y) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::zz_printCompact(std::ostream &stream=CASADI_COUT) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::zz_project(const Sparsity &sparsity) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::zz_pw_const(const Matrix< DataType > &tval, const Matrix< DataType > &val) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::zz_pw_lin(const Matrix< DataType > &tval, const Matrix< DataType > &val) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::zz_qr(Matrix< DataType > &Q, Matrix< DataType > &R) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::zz_quad_form(const Matrix< DataType > &A) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::zz_ramp() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::zz_rdivide(const Matrix< DataType > &y) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::zz_rectangle() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::zz_reshape(const Sparsity &sp) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::zz_reshape(int nrow, int ncol) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::zz_sign() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::zz_simplify() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::zz_sin() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::zz_sinh() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::zz_solve(const Matrix< DataType > &b) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::zz_sparsify(double tol=0) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::zz_spy() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::zz_sqrt() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::zz_substitute(const Matrix< DataType > &v, const Matrix< DataType > &vdef) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::zz_substitute(const std::vector< Matrix< DataType > > &ex, const std::vector< Matrix< DataType > > &v, const std::vector< Matrix< DataType > > &vdef) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::zz_substituteInPlace(const std::vector< Matrix< DataType > > &v, std::vector< Matrix< DataType > > &vdef, std::vector< Matrix< DataType > > &ex, bool reverse=false) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::zz_sumAll() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::zz_sumCols() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::zz_sumRows() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::zz_tan() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::zz_tangent(const Matrix< DataType > &arg) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::zz_tanh() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::zz_taylor(const Matrix< DataType > &x, const Matrix< DataType > &a=0, int order=1) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::zz_times(const Matrix< DataType > &y) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::zz_trace() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::zz_triangle() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::zz_unite(const Matrix< DataType > &B) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::zz_vecNZ() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::zz_vertcat(const std::vector< Matrix< DataType > > &v) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::zz_vertsplit(const std::vector< int > &offset) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< T >::at(int k) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< T >::at(int k) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< T >::back() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< T >::back() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< T >::begin() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< T >::begin() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< T >::end() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< T >::end() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< T >::front() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< T >::front() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< T >::getBV(Matrix< DataType > &val) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< T >::getElement(int rr, int cc=0) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< T >::ptr() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< T >::ptr() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< T >::rbegin() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< T >::rbegin() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< T >::rend() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< T >::rend() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< T >::sparsity() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MinusInfSX::getValue() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MinusInfSX::isMinusInf() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MinusOneSX::getIntValue() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MinusOneSX::getValue() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MinusOneSX::isInteger() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MinusOneSX::isMinusOne() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MultipleOutput::getOutput(int oind) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MultipleOutput::isMultipleOutput() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MultipleOutput::nout() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MultipleOutput::sparsity(int oind) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Multiplication::clone() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Multiplication::eval(const cpv_MX &arg, const pv_MX &res) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Multiplication::evalAdj(const std::vector< pv_MX > &aseed, const std::vector< pv_MX > &asens) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Multiplication::evalFwd(const std::vector< cpv_MX > &fseed, const std::vector< pv_MX > &fsens) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Multiplication::generate(std::ostream &stream, const std::vector< int > &arg, const std::vector< int > &res, CodeGenerator &gen) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Multiplication::getOp() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Multiplication::nTmp(size_t &ni, size_t &nr) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Multiplication::numInplace() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Multiplication::printPart(std::ostream &stream, int part) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::NanSX::getValue() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::NanSX::isNan() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Newton::clone() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Newton::create(const Function &f) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Newton::init() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Newton::solveNonLinear() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::NlpSolverInternal::checkInitialBounds() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::NlpSolverInternal::checkInputs() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::NlpSolverInternal::getGradF() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::NlpSolverInternal::getGradLag() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::NlpSolverInternal::getHessLag() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::NlpSolverInternal::getJacF() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::NlpSolverInternal::getJacG() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::NlpSolverInternal::getReducedHessian() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::NlpSolverInternal::getSpHessLag() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::NlpSolverInternal::gradF() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::NlpSolverInternal::gradLag() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::NlpSolverInternal::hessLag() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::NlpSolverInternal::jacF() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::NlpSolverInternal::jacG() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::NlpSolverInternal::reportConstraints(std::ostream &stream=std::cout) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::NlpSolverInternal::setOptionsFromFile(const std::string &file) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::NlpSolverInternal::setQPOptions() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::NlpSolverInternal::spHessLag() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Norm1::clone() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Norm1::getOp() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Norm1::printPart(std::ostream &stream, int part) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Norm2::clone() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Norm2::getOp() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Norm2::printPart(std::ostream &stream, int part) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::NormF::clone() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::NormF::eval(const cpv_MX &input, const pv_MX &output) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::NormF::evalAdj(const std::vector< pv_MX > &adjSeed, const std::vector< pv_MX > &adjSens) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::NormF::evalFwd(const std::vector< cpv_MX > &fwdSeed, const std::vector< pv_MX > &fwdSens) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::NormF::generate(std::ostream &stream, const std::vector< int > &arg, const std::vector< int > &res, CodeGenerator &gen) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::NormF::getOp() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::NormF::printPart(std::ostream &stream, int part) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::NormInf::clone() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::NormInf::getOp() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::NormInf::printPart(std::ostream &stream, int part) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::OldCollocationIntegrator::clone() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::OldCollocationIntegrator::create(const Function &f, const Function &g) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::OldCollocationIntegrator::init() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::OldCollocationIntegrator::integrate(double t_out) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::OldCollocationIntegrator::integrateB(double t_out) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::OldCollocationIntegrator::reset() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::OldCollocationIntegrator::resetB() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::OneSX::getIntValue() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::OneSX::getValue() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::OneSX::isInteger() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::OneSX::isOne() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::OoqpInterface::clone() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::OoqpInterface::evaluate() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::OoqpInterface::init() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::OptionsFunctionality::getOptionAllowedIndex(const std::string &name) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::OptionsFunctionality::getOptionEnumValue(const std::string &name) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::OptionsFunctionality::setOptionByAllowedIndex(const std::string &name, int i) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::OptionsFunctionality::setOptionByEnumValue(const std::string &name, int v) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::OptionsFunctionalityNode::addOption(const std::string &str, const opt_type &type, const GenericType &def_val, const std::string &desc, const std::string &allowed_vals, bool inherit=false) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::OptionsFunctionalityNode::addOption(const std::string &str, const opt_type &type, const GenericType &def_val=GenericType(), const std::string &desc="n/a", const std::vector< GenericType > &allowed_vals=std::vector< GenericType >(), bool inherit=false, std::vector< int > enum_values=std::vector< int >(), std::vector< std::string > enum_descr=std::vector< std::string >()) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::OptionsFunctionalityNode::copyOptions(const OptionsFunctionality &obj, bool skipUnknown=false) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::OptionsFunctionalityNode::dictionary() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::OptionsFunctionalityNode::getBestMatches(const std::string &name, std::vector< std::string > &suggestions, int amount=5) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::OptionsFunctionalityNode::getOption(const std::string &str) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::OptionsFunctionalityNode::getOptionAllowed(const std::string &str) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::OptionsFunctionalityNode::getOptionAllowedIndex(const std::string &name) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::OptionsFunctionalityNode::getOptionDefault(const std::string &str) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::OptionsFunctionalityNode::getOptionDescription(const std::string &str) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::OptionsFunctionalityNode::getOptionEnumValue(const std::string &name) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::OptionsFunctionalityNode::getOptionNames() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::OptionsFunctionalityNode::getOptionType(const std::string &str) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::OptionsFunctionalityNode::getOptionTypeName(const std::string &str) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::OptionsFunctionalityNode::hasOption(const std::string &str) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::OptionsFunctionalityNode::hasSetOption(const std::string &str) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::OptionsFunctionalityNode::print(std::ostream &stream) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::OptionsFunctionalityNode::printOption(const std::string &name, std::ostream &stream=CASADI_COUT) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::OptionsFunctionalityNode::printOptions(std::ostream &stream=CASADI_COUT) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::OptionsFunctionalityNode::repr(std::ostream &stream) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::OptionsFunctionalityNode::setOption(const Dictionary &dict, bool skipUnknown=false) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::OptionsFunctionalityNode::setOption(const std::string &str, const GenericType &val) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::OptionsFunctionalityNode::setOptionByAllowedIndex(const std::string &name, int i) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::OptionsFunctionalityNode::setOptionByEnumValue(const std::string &name, int v) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::OutputNode::clone() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::OutputNode::getFunctionInput() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::OutputNode::getFunctionOutput() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::OutputNode::getHorzcat(const std::vector< MX > &x) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::OutputNode::getOp() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::OutputNode::getVertcat(const std::vector< MX > &x) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::OutputNode::isNonLinear() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::OutputNode::isOutputNode() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::OutputNode::printPart(std::ostream &stream, int part) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::PrintableObject< SXElement  >::getDescription() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::PrintableObject< SXElement  >::getRepresentation() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ProfilingType() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ProfilingType< ProfilingData_ENTRY >() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ProfilingType< ProfilingData_EXIT >() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ProfilingType< ProfilingData_IO >() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ProfilingType< ProfilingData_NAME >() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ProfilingType< ProfilingData_SOURCE >() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ProfilingType< ProfilingData_TIMELINE >() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::QcqpSolverInternal::checkInputs() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::QcqpSolverInternal::setQPOptions() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::QcqpSolverInternal::solve() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::QcqpToSocp::clone() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::QcqpToSocp::evaluate() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::QcqpToSocp::init() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::QpSolverInternal::checkInputs() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::QpSolverInternal::generateNativeCode(std::ostream &file) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::QpSolverInternal::setLPOptions() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::QpSolverInternal::solve() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::QpToImplicit::clone() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::QpToImplicit::create(const Function &f) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::QpToImplicit::init() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::QpToImplicit::solveNonLinear() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::QpToNlp::clone() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::QpToNlp::evaluate() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::QpToNlp::init() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::QpToQcqp::clone() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::QpToQcqp::evaluate() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::QpToQcqp::init() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::QpoasesInterface::clone() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::QpoasesInterface::evaluate() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::QpoasesInterface::init() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::RealtypeSX::getIntValue() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::RealtypeSX::getValue() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::RealtypeSX::isAlmostZero(double tol) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Reshape::clone() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Reshape::eval(const cpv_MX &input, const pv_MX &output) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Reshape::evalAdj(const std::vector< pv_MX > &adjSeed, const std::vector< pv_MX > &adjSens) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Reshape::evalFwd(const std::vector< cpv_MX > &fwdSeed, const std::vector< pv_MX > &fwdSens) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Reshape::generate(std::ostream &stream, const std::vector< int > &arg, const std::vector< int > &res, CodeGenerator &gen) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Reshape::getOp() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Reshape::getReshape(const Sparsity &sp) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Reshape::getTranspose() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Reshape::numInplace() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Reshape::printPart(std::ostream &stream, int part) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::RkIntegrator::clone() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::RkIntegrator::create(const Function &f, const Function &g) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::RkIntegrator::init() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::RkIntegrator::setupFG() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXElement::assignIfDuplicate(const SXElement &scalar, int depth=1) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXElement::assignNoDelete(const SXElement &scalar) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXElement::constpow(const Matrix< SXElement > &n) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXElement::constpow(const SXElement &y) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXElement::get() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXElement::getDep(int ch=0) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXElement::getIntValue() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXElement::getName() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXElement::getNdeps() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXElement::getOp() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXElement::getTemp() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXElement::getValue() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXElement::hasDep() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXElement::inv() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXElement::isAlmostZero(double tol) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXElement::isCommutative() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXElement::isConstant() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXElement::isDoubled() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXElement::isInf() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXElement::isInteger() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXElement::isLeaf() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXElement::isMinusInf() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXElement::isMinusOne() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXElement::isNan() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXElement::isNonNegative() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXElement::isNull() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXElement::isOne() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXElement::isOp(int op) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXElement::isRegular() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXElement::isSymbolic() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXElement::isZero() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXElement::mark() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXElement::marked() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXElement::print(std::ostream &stream, long &remaining_calls) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXElement::print(std::ostream &stream=CASADI_COUT, bool trailing_newline=true) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXElement::printme(const SXElement &y) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXElement::repr(std::ostream &stream=CASADI_COUT, bool trailing_newline=true) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXElement::setTemp(int t) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXElement::sq() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXElement::zz_abs() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXElement::zz_acos() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXElement::zz_acosh() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXElement::zz_and(const SXElement &y) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXElement::zz_asin() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXElement::zz_asinh() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXElement::zz_atan() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXElement::zz_atan2(const Matrix< SXElement > &b) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXElement::zz_atan2(const SXElement &y) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXElement::zz_atanh() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXElement::zz_ceil() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXElement::zz_cos() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXElement::zz_cosh() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXElement::zz_eq(const SXElement &y) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXElement::zz_erf() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXElement::zz_erfinv() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXElement::zz_exp() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXElement::zz_floor() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXElement::zz_if_else_zero(const SXElement &y) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXElement::zz_isEqual(const SXElement &scalar, int depth=0) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXElement::zz_le(const SXElement &y) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXElement::zz_log() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXElement::zz_log10() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXElement::zz_lt(const SXElement &y) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXElement::zz_max(const Matrix< SXElement > &b) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXElement::zz_max(const SXElement &y) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXElement::zz_min(const Matrix< SXElement > &b) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXElement::zz_min(const SXElement &y) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXElement::zz_minus(const SXElement &y) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXElement::zz_mod(const SXElement &y) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXElement::zz_mpower(const SXElement &b) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXElement::zz_mul(const SXElement &y) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXElement::zz_ne(const SXElement &y) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXElement::zz_not() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXElement::zz_or(const SXElement &y) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXElement::zz_plus(const SXElement &y) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXElement::zz_power(const SXElement &b) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXElement::zz_rdivide(const SXElement &y) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXElement::zz_sign() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXElement::zz_simplify() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXElement::zz_sin() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXElement::zz_sinh() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXElement::zz_sqrt() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXElement::zz_tan() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXElement::zz_tanh() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXElement::zz_times(const SXElement &y) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXFunction::algorithm() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXNode::dep(int i) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXNode::dep(int i) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXNode::getIntValue() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXNode::getName() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXNode::getOp() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXNode::getValue() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXNode::hasDep() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXNode::init() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXNode::isAlmostZero(double tol) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXNode::isConstant() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXNode::isInf() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXNode::isInteger() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXNode::isMinusInf() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXNode::isMinusOne() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXNode::isNan() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXNode::isOne() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXNode::isSmooth() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXNode::isSymbolic() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXNode::isZero() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXNode::mark() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXNode::marked() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXNode::ndep() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXNode::print(std::ostream &stream) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXNode::print(std::ostream &stream, long &remaining_calls) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Scpgen::clone() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Scpgen::dualInfeasibility() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Scpgen::eval_exp() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Scpgen::eval_mat() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Scpgen::eval_res() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Scpgen::eval_vec() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Scpgen::evaluate() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Scpgen::getQpSolver() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Scpgen::init() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Scpgen::line_search(int &ls_iter, bool &ls_success) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Scpgen::primalInfeasibility() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Scpgen::printIteration(std::ostream &stream) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Scpgen::printIteration(std::ostream &stream, int iter, double obj, double pr_inf, double du_inf, double reg, int ls_trials, bool ls_success) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Scpgen::regularize() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Scpgen::solve_qp() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SdpSolverInternal::checkInputs() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SdpSolverInternal::printProblem(std::ostream &stream=std::cout) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SdpSolverInternal::setSOCPOptions() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SdpSolverInternal::solve() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SdqpSolverInternal::printProblem(std::ostream &stream=std::cout) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SdqpSolverInternal::setSOCQPOptions() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SdqpSolverInternal::solve() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SdqpToSdp::clone() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SdqpToSdp::evaluate() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SdqpToSdp::init() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SetNonzeros< Add >::clone() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SetNonzeros< Add >::eval(const cpv_MX &input, const pv_MX &output) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SetNonzeros< Add >::evalAdj(const std::vector< pv_MX > &adjSeed, const std::vector< pv_MX > &adjSens) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SetNonzeros< Add >::evalFwd(const std::vector< cpv_MX > &fwdSeed, const std::vector< pv_MX > &fwdSens) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SetNonzeros< Add >::getAll() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SetNonzeros< Add >::getOp() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SetNonzeros< Add >::mapping() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SetNonzeros< Add >::numInplace() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SetNonzerosSlice2< Add >::clone() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SetNonzerosSlice2< Add >::generate(std::ostream &stream, const std::vector< int > &arg, const std::vector< int > &res, CodeGenerator &gen) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SetNonzerosSlice2< Add >::getAll() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SetNonzerosSlice2< Add >::printPart(std::ostream &stream, int part) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SetNonzerosSlice< Add >::clone() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SetNonzerosSlice< Add >::generate(std::ostream &stream, const std::vector< int > &arg, const std::vector< int > &res, CodeGenerator &gen) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SetNonzerosSlice< Add >::getAll() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SetNonzerosSlice< Add >::isAssignment() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SetNonzerosSlice< Add >::printPart(std::ostream &stream, int part) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SetNonzerosSlice< Add >::simplifyMe(MX &ex) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SetNonzerosVector< Add >::clone() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SetNonzerosVector< Add >::generate(std::ostream &stream, const std::vector< int > &arg, const std::vector< int > &res, CodeGenerator &gen) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SetNonzerosVector< Add >::getAll() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SetNonzerosVector< Add >::printPart(std::ostream &stream, int part) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SetSparse::clone() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SetSparse::eval(const cpv_MX &input, const pv_MX &output) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SetSparse::evalAdj(const std::vector< pv_MX > &adjSeed, const std::vector< pv_MX > &adjSens) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SetSparse::evalFwd(const std::vector< cpv_MX > &fwdSeed, const std::vector< pv_MX > &fwdSens) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SetSparse::generate(std::ostream &stream, const std::vector< int > &arg, const std::vector< int > &res, CodeGenerator &gen) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SetSparse::getOp() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SetSparse::nTmp(size_t &ni, size_t &nr) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SetSparse::printPart(std::ostream &stream, int part) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SharedObject::assertInit() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SharedObject::clone() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SharedObject::get() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SharedObject::get() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SharedObject::getCount() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SharedObject::printPtr(std::ostream &stream=CASADI_COUT) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SharedObject::swap(SharedObject &other) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SharedObject::weak() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SharedObjectNode::assertInit() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SharedObjectNode::clone() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SharedObjectNode::getCount() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SharedObjectNode::init() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SharedObjectNode::isInit() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SharedObjectNode::print(std::ostream &stream) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SharedObjectNode::repr(std::ostream &stream) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SharedObjectNode::weak() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SimpleHomotopyNlp::clone() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SimpleHomotopyNlp::evaluate() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SimpleHomotopyNlp::init() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SnoptInterface::clone() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SnoptInterface::evaluate() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SnoptInterface::formatStatus(int status) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SnoptInterface::init() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SnoptInterface::setOptionsFromFile(const std::string &file) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SnoptInterface::setQPOptions() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SocpSolverInternal::checkInputs() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SocpSolverInternal::printProblem(std::ostream &stream=std::cout) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SocpSolverInternal::solve() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SocpToSdp::clone() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SocpToSdp::evaluate() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SocpToSdp::init() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Solve< Tr >::clone() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Solve< Tr >::eval(const cpv_MX &arg, const pv_MX &res) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Solve< Tr >::evalAdj(const std::vector< pv_MX > &adjSeed, const std::vector< pv_MX > &adjSens) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Solve< Tr >::evalFwd(const std::vector< cpv_MX > &fwdSeed, const std::vector< pv_MX > &fwdSens) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Solve< Tr >::getFunction() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Solve< Tr >::getOp() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Solve< Tr >::nTmp(size_t &ni, size_t &nr) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Solve< Tr >::numInplace() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Solve< Tr >::print(std::ostream &stream, long &remaining_calls) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Solve< Tr >::printPart(std::ostream &stream, int part) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SparseStorage< DataType >::clear() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SparseStorage< DataType >::data() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SparseStorage< DataType >::data() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SparseStorage< DataType >::elem(int rr, int cc) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SparseStorage< DataType >::hasNZ(int rr, int cc) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SparseStorage< DataType >::reserve(int nnz) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SparseStorage< DataType >::reserve(int nnz, int ncol) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SparseStorage< DataType >::resize(int nrow, int ncol) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SparseStorage< DataType >::sparsity() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SparseStorage< DataType >::sparsityRef() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SparseStorage< DataType >::toScalar() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Sparsity::colind() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Sparsity::find(bool ind1=SWIG_IND1) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Sparsity::find(std::vector< int > &loc, bool ind1=false) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Sparsity::reCache() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Sparsity::row() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Sparsity::zz_blockcat(const std::vector< std::vector< Sparsity > > &v) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Sparsity::zz_diagcat(const std::vector< Sparsity > &v) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Sparsity::zz_diagsplit(const std::vector< int > &offset1, const std::vector< int > &offset2) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Sparsity::zz_horzcat(const std::vector< Sparsity > &sp) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Sparsity::zz_horzsplit(const std::vector< int > &output_offset) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Sparsity::zz_kron(const Sparsity &b) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Sparsity::zz_mtimes(const Sparsity &Y, const Sparsity &Z) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Sparsity::zz_mtimes(const Sparsity &y) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Sparsity::zz_norm_0_mul(const Sparsity &B) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Sparsity::zz_reshape(const Sparsity &sp) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Sparsity::zz_reshape(int nrow, int ncol) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Sparsity::zz_sprank() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Sparsity::zz_tril(bool includeDiagonal=true) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Sparsity::zz_triu(bool includeDiagonal=true) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Sparsity::zz_vecNZ() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Sparsity::zz_vertcat(const std::vector< Sparsity > &sp) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Sparsity::zz_vertsplit(const std::vector< int > &output_offset) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Split::generate(std::ostream &stream, const std::vector< int > &arg, const std::vector< int > &res, CodeGenerator &gen) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Split::nout() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Split::sparsity(int oind) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SqicInterface::clone() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SqicInterface::evaluate() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SqicInterface::generateNativeCode(std::ostream &file) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SqicInterface::init() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Sqpmethod::clone() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Sqpmethod::eval_f(const std::vector< double > &x, double &f) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Sqpmethod::eval_g(const std::vector< double > &x, std::vector< double > &g) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Sqpmethod::eval_grad_f(const std::vector< double > &x, double &f, std::vector< double > &grad_f) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Sqpmethod::eval_h(const std::vector< double > &x, const std::vector< double > &lambda, double sigma, Matrix< double > &H) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Sqpmethod::eval_jac_g(const std::vector< double > &x, std::vector< double > &g, Matrix< double > &J) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Sqpmethod::evaluate() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Sqpmethod::getQpSolver() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Sqpmethod::getRegularization(const Matrix< double > &H) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Sqpmethod::init() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Sqpmethod::primalInfeasibility(const std::vector< double > &x, const std::vector< double > &lbx, const std::vector< double > &ubx, const std::vector< double > &g, const std::vector< double > &lbg, const std::vector< double > &ubg) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Sqpmethod::printIteration(std::ostream &stream) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Sqpmethod::printIteration(std::ostream &stream, int iter, double obj, double pr_inf, double du_inf, double dx_norm, double reg, int ls_trials, bool ls_success) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Sqpmethod::regularize(Matrix< double > &H, double reg) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Sqpmethod::reset_h() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Sqpmethod::solve_QP(const Matrix< double > &H, const std::vector< double > &g, const std::vector< double > &lbx, const std::vector< double > &ubx, const Matrix< double > &A, const std::vector< double > &lbA, const std::vector< double > &ubA, std::vector< double > &x_opt, std::vector< double > &lambda_x_opt, std::vector< double > &lambda_A_opt) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::StabilizedQpSolverInternal::checkInputs() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::StabilizedQpSolverInternal::setLPOptions() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::StabilizedQpSolverInternal::solve() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::StabilizedQpToQp::clone() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::StabilizedQpToQp::evaluate() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::StabilizedQpToQp::generateNativeCode(std::ostream &file) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::StabilizedQpToQp::init() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::StabilizedSqicInterface::clone() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::StabilizedSqicInterface::evaluate() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::StabilizedSqicInterface::generateNativeCode(std::ostream &file) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::StabilizedSqicInterface::init() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::StabilizedSqp::clone() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::StabilizedSqp::eval_f(const std::vector< double > &x, double &f) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::StabilizedSqp::eval_g(const std::vector< double > &x, std::vector< double > &g) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::StabilizedSqp::eval_grad_f(const std::vector< double > &x, double &f, std::vector< double > &grad_f) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::StabilizedSqp::eval_h(const std::vector< double > &x, const std::vector< double > &lambda, double sigma, Matrix< double > &H) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::StabilizedSqp::eval_jac_g(const std::vector< double > &x, std::vector< double > &g, Matrix< double > &J) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::StabilizedSqp::evaluate() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::StabilizedSqp::getRegularization(const Matrix< double > &H) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::StabilizedSqp::getStabilizedQpSolver() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::StabilizedSqp::init() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::StabilizedSqp::mat_vec(const std::vector< double > &x, const DMatrix &A, std::vector< double > &y) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::StabilizedSqp::mat_vectran(const std::vector< double > &x, const DMatrix &A, std::vector< double > &y) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::StabilizedSqp::meritfg() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::StabilizedSqp::norm1matrix(const DMatrix &A) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::StabilizedSqp::primalInfeasibility(const std::vector< double > &x, const std::vector< double > &lbx, const std::vector< double > &ubx, const std::vector< double > &g, const std::vector< double > &lbg, const std::vector< double > &ubg) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::StabilizedSqp::printIteration(std::ostream &stream) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::StabilizedSqp::printIteration(std::ostream &stream, int iter, double obj, double pr_inf, double du_inf, double dx_norm, double reg, double TRdelta, int ls_trials, bool ls_success, char info) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::StabilizedSqp::regularize(Matrix< double > &H, double reg) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::StabilizedSqp::reset_h() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::StabilizedSqp::solve_QP(const Matrix< double > &H, const std::vector< double > &g, const std::vector< double > &lbx, const std::vector< double > &ubx, const Matrix< double > &A, const std::vector< double > &lbA, const std::vector< double > &ubA, std::vector< double > &x_opt, std::vector< double > &lambda_x_opt, std::vector< double > &lambda_A_opt, double muR, const std::vector< double > &mu, const std::vector< double > &muE) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SubAssign::clone() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SubAssign::eval(const cpv_MX &input, const pv_MX &output) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SubAssign::evalAdj(const std::vector< pv_MX > &adjSeed, const std::vector< pv_MX > &adjSens) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SubAssign::evalFwd(const std::vector< cpv_MX > &fwdSeed, const std::vector< pv_MX > &fwdSens) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SubAssign::generate(std::ostream &stream, const std::vector< int > &arg, const std::vector< int > &res, CodeGenerator &gen) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SubAssign::getOp() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SubAssign::printPart(std::ostream &stream, int part) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SubRef::clone() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SubRef::eval(const cpv_MX &input, const pv_MX &output) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SubRef::evalAdj(const std::vector< pv_MX > &adjSeed, const std::vector< pv_MX > &adjSens) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SubRef::evalFwd(const std::vector< cpv_MX > &fwdSeed, const std::vector< pv_MX > &fwdSens) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SubRef::generate(std::ostream &stream, const std::vector< int > &arg, const std::vector< int > &res, CodeGenerator &gen) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SubRef::getOp() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SubRef::printPart(std::ostream &stream, int part) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SundialsInterface::getBandwidth() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SundialsInterface::getBandwidthB() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SundialsInterface::getJac() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SundialsInterface::getJacB() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SundialsInterface::init() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SundialsInterface::reset() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SundialsInterface::setStopTime(double tf) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SymbolicMX::clone() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SymbolicMX::eval(const cpv_MX &input, const pv_MX &output) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SymbolicMX::evalAdj(const std::vector< pv_MX > &adjSeed, const std::vector< pv_MX > &adjSens) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SymbolicMX::evalFwd(const std::vector< cpv_MX > &fwdSeed, const std::vector< pv_MX > &fwdSens) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SymbolicMX::getName() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SymbolicMX::getOp() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SymbolicMX::printPart(std::ostream &stream, int part) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SymbolicQr::clone() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SymbolicQr::evaluateSXGen(const SXPtrV &input, SXPtrV &output, bool tr) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SymbolicQr::generateBody(std::ostream &stream, const std::string &type, CodeGenerator &gen) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SymbolicQr::generateDeclarations(std::ostream &stream, const std::string &type, CodeGenerator &gen) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SymbolicQr::init() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SymbolicQr::prepare() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SymbolicSX::getName() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SymbolicSX::getOp() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SymbolicSX::isSymbolic() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::TinyXmlInterface::clone() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::TinyXmlInterface::parse(const std::string &filename) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Transpose::clone() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Transpose::eval(const cpv_MX &input, const pv_MX &output) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Transpose::evalAdj(const std::vector< pv_MX > &adjSeed, const std::vector< pv_MX > &adjSens) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Transpose::evalFwd(const std::vector< cpv_MX > &fwdSeed, const std::vector< pv_MX > &fwdSens) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Transpose::generate(std::ostream &stream, const std::vector< int > &arg, const std::vector< int > &res, CodeGenerator &gen) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Transpose::getOp() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Transpose::getTranspose() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Transpose::nTmp(size_t &ni, size_t &nr) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Transpose::printPart(std::ostream &stream, int part) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::UnaryMX::clone() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::UnaryMX::eval(const cpv_MX &input, const pv_MX &output) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::UnaryMX::evalAdj(const std::vector< pv_MX > &adjSeed, const std::vector< pv_MX > &adjSens) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::UnaryMX::evalFwd(const std::vector< cpv_MX > &fwdSeed, const std::vector< pv_MX > &fwdSens) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::UnaryMX::generate(std::ostream &stream, const std::vector< int > &arg, const std::vector< int > &res, CodeGenerator &gen) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::UnaryMX::getBinary(int op, const MX &y, bool scX, bool scY) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::UnaryMX::getOp() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::UnaryMX::getUnary(int op) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::UnaryMX::isUnaryOp() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::UnaryMX::numInplace() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::UnaryMX::printPart(std::ostream &stream, int part) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::UnarySX::dep(int i) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::UnarySX::dep(int i) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::UnarySX::getOp() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::UnarySX::hasDep() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::UnarySX::isSmooth() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::UnarySX::ndep() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::UnarySX::print(std::ostream &stream, long &remaining_calls) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Vertcat::clone() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Vertcat::eval(const cpv_MX &arg, const pv_MX &res) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Vertcat::evalAdj(const std::vector< pv_MX > &adjSeed, const std::vector< pv_MX > &adjSens) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Vertcat::evalFwd(const std::vector< cpv_MX > &fwdSeed, const std::vector< pv_MX > &fwdSens) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Vertcat::getOp() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Vertcat::printPart(std::ostream &stream, int part) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Vertsplit::clone() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Vertsplit::eval(const cpv_MX &input, const pv_MX &output) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Vertsplit::evalAdj(const std::vector< pv_MX > &adjSeed, const std::vector< pv_MX > &adjSens) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Vertsplit::evalFwd(const std::vector< cpv_MX > &fwdSeed, const std::vector< pv_MX > &fwdSens) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Vertsplit::getOp() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Vertsplit::getVertcat(const std::vector< MX > &x) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Vertsplit::printPart(std::ostream &stream, int part) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::WeakRef::alive() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::WeakRef::shared() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::WorhpInterface::clone() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::WorhpInterface::evaluate() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::WorhpInterface::formatStatus(int status) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::WorhpInterface::init() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::WorhpInterface::passOptions() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::WorhpInterface::reset() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::WorhpInterface::setOptionsFromFile(const std::string &file) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::WorhpInterface::setQPOptions() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Wrapper< Derived >::checkDimensions() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Wrapper< Derived >::evaluate() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Wrapper< DleToLrDle  >::checkDimensions() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Wrapper< DpleToLrDple  >::checkDimensions() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Wrapper< LrDleToDle  >::checkDimensions() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Wrapper< LrDpleToDple  >::checkDimensions() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::XmlFile::parse(const std::string &filename) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::XmlFileInternal::print(std::ostream &stream) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::XmlNode::checkName(const std::string &str) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::XmlNode::dump(std::ostream &stream, int indent=0) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::XmlNode::getAttribute(const std::string &attribute_name) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::XmlNode::getName() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::XmlNode::getText() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::XmlNode::getText(T &val) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::XmlNode::hasAttribute(const std::string &attribute_name) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::XmlNode::hasChild(const std::string &childname) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::XmlNode::readAttribute(const std::string &attribute_name, T &val, bool assert_existance=true) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::XmlNode::setAttribute(const std::string &attribute_name, const std::string &attribute) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::XmlNode::setName(const std::string &name) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::XmlNode::size() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ZeroByZero::clone() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ZeroByZero::generate(std::ostream &stream, const std::vector< int > &arg, const std::vector< int > &res, CodeGenerator &gen) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ZeroByZero::getBinary(int op, const MX &y, bool ScX, bool ScY) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ZeroByZero::getGetNonzeros(const Sparsity &sp, const std::vector< int > &nz) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ZeroByZero::getMatrixValue() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ZeroByZero::getReshape(const Sparsity &sp) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ZeroByZero::getSetNonzeros(const MX &y, const std::vector< int > &nz) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ZeroByZero::getSetSparse(const Sparsity &sp) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ZeroByZero::getTranspose() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ZeroByZero::getUnary(int op) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ZeroByZero::getValue() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ZeroByZero::printPart(std::ostream &stream, int part) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ZeroSX::getIntValue() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ZeroSX::getValue() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ZeroSX::isAlmostZero(double tol) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ZeroSX::isInteger() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ZeroSX::isZero() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::acosh(double x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::asinh(double x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::atanh(double x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::casadi_load_linearsolver_csparsecholesky() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::check_exposed(T t) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::collocationInterpolators(const std::vector< double > &tau_root, std::vector< std::vector< double > > &C, std::vector< double > &D) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::collocationPointsL(int order, const std::string &scheme) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::constpow(const T &x, const T &n) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::copysign(const T &x, const T &y) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::copysign(double x, double y) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::cumsum(const std::vector< T > &values) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::cumsum0(const std::vector< T > &values) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::deepcopy(const A &a) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::erf(double x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::erfinv(const T &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::erfinv(double x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::fabs(int x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::fmax(double x, double y) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::fmax(int x, int y) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::fmin(double x, double y) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::fmin(int x, int y) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::getDescription(const std::vector< T > &v) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::getPtr(const std::vector< T > &v) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::getPtr(std::vector< T > &v) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::getRealTime() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::getRepresentation(const std::vector< T > &v) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::get_bvec_t(const std::vector< T > &v) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::get_bvec_t(const std::vector< double > &v) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::get_bvec_t(std::vector< T > &v) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::get_bvec_t(std::vector< double > &v) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::hash_combine(std::size_t &seed, T v) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::hash_combine(std::size_t &seed, const std::vector< int > &v) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::hash_sparsity(int nrow, int ncol, const std::vector< int > &colind, const std::vector< int > &row) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::hash_value(T v) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::if_else_zero(const T &x, const T &y) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::if_else_zero(double x, double y) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::inner_prod(const std::vector< T > &a, const std::vector< T > &b) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::isEqual(double x, double y, int depth=0) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::is_a(const SharedObject &A) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::isinf(double x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::isnan(double x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::iszero(double x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::linspace(std::vector< T > &v, const F &first, const L &last) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::makeVector(int size, int ind0=-1, const T &val0=T(), int ind1=-1, const T &val1=T(), int ind2=-1, const T &val2=T(), int ind3=-1, const T &val3=T(), int ind4=-1, const T &val4=T(), int ind5=-1, const T &val5=T(), int ind6=-1, const T &val6=T(), int ind7=-1, const T &val7=T(), int ind8=-1, const T &val8=T(), int ind9=-1, const T &val9=T(), int ind10=-1, const T &val10=T(), int ind11=-1, const T &val11=T(), int ind12=-1, const T &val12=T(), int ind13=-1, const T &val13=T(), int ind14=-1, const T &val14=T(), int ind15=-1, const T &val15=T(), int ind16=-1, const T &val16=T(), int ind17=-1, const T &val17=T(), int ind18=-1, const T &val18=T(), int ind19=-1, const T &val19=T()) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::matrixName< SXElement >() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::norm_1(const std::vector< T > &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::norm_2(const std::vector< T > &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::norm_inf(const std::vector< T > &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::operation_checker(unsigned int op) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::print(const std::vector< T > &v, std::ostream &stream=CASADI_COUT) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::printme(const T &x, const T &y) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::printme(double x, double y) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::profileWrite(std::ofstream &f, const T &s) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::profileWriteBare(std::ofstream &f, const T &s) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ptrVec(const std::vector< T > &v) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ptrVec(const std::vector< std::vector< T > > &v) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ptrVec(std::vector< T > &v) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ptrVec(std::vector< std::vector< T > > &v) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::range(int start, int stop, int step=1, int len=std::numeric_limits< int >::max()) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::range(int stop) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::repr(const std::vector< T > &v, std::ostream &stream=CASADI_COUT) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::reverse(const std::vector< T > &v) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::shared_cast(SharedObject &A) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::shared_cast(const SharedObject &A) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::sign(double x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::sort(const std::vector< T > &values, std::vector< T > &sorted_values, std::vector< int > &indices, bool invert_indices=false) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::sq(const T &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::toVector(const T &v0) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::toVector(const T &v0, const T &v1) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::toVector(const T &v0, const T &v1, const T &v2) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::twice(const T &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::vector_slice(const std::vector< T > &v, const std::vector< int > &i) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  snoptProblem::setLog(isnLog snLog, isnLog2 snLog2, isqLog sqLog) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  snoptProblem::setSTOP(isnSTOP snSTOP) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  snoptProblem::solve(int starttype) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  snoptProblemA::computeJac() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  snoptProblemA::setNeA(int neA) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  snoptProblemA::setNeG(int neG) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  snoptProblemA::setObjective(int ObjRow, double ObjAdd) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  snoptProblemA::setProblemSize(int n, int neF) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  snoptProblemA::setUserFun(snFunA usrfun) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  snoptProblemA::setWorkspace() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  snoptProblemA::solve(int starttype) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  snoptProblemB::setFuncon(snConB funcon) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  snoptProblemB::setFunobj(snObjB funobj) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  snoptProblemB::solve(int starttype) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  snoptProblemC::setObjective(int iObj, double ObjAdd) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  snoptProblemC::setProblemSize(int m, int n, int nnCon, int nnJac, int nnObj) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  snoptProblemC::setUserFun(snFunC usrfun) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  snoptProblemC::setWorkspace() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  snoptProblemC::solve(int starttype) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  sqicDestroy() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  summaryOff() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  summaryOn() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::Assertion::Assertion(const MX &x, const MX &y, const std::string &s) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::BinaryMX< ScX, ScY >::BinaryMX(Operation op, const MX &x, const MX &y) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::CLEInputIOSchemeVector< M >::CLEInputIOSchemeVector(const std::vector< M > &t) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::CLEOutputIOSchemeVector< M >::CLEOutputIOSchemeVector(const std::vector< M > &t) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::CallFunction::CallFunction(const Function &fcn, std::vector< MX > arg) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::CleStructIOSchemeVector< T >::CleStructIOSchemeVector(const std::vector< M > &t) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::CollocationIntegrator::CollocationIntegrator(const Function &f, const Function &g) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::Concat::Concat(const std::vector< MX > &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::Constant< Value >::Constant(const Sparsity &sp, Value v=Value()) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::ConstantDMatrix::ConstantDMatrix(const Matrix< double > &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::ConstantMX::ConstantMX(const Sparsity &sp) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::ControlSimulatorInputIOSchemeVector< M >::ControlSimulatorInputIOSchemeVector(const std::vector< M > &t) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::ControlledDAEInputIOSchemeVector< M >::ControlledDAEInputIOSchemeVector(const std::vector< M > &t) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::CplexInterface::CplexInterface() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::CplexInterface::CplexInterface(const std::vector< Sparsity > &st) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::CsparseInterface::CsparseInterface(const CsparseInterface &linsol) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::CsparseInterface::CsparseInterface(const Sparsity &sp, int nrhs) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::CvodesInterface::CvodesInterface(const Function &f, const Function &g) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::DAEInputIOSchemeVector< M >::DAEInputIOSchemeVector(const std::vector< M > &t) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::DAEOutputIOSchemeVector< M >::DAEOutputIOSchemeVector(const std::vector< M > &t) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::DLEInputIOSchemeVector< M >::DLEInputIOSchemeVector(const std::vector< M > &t) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::DLEOutputIOSchemeVector< M >::DLEOutputIOSchemeVector(const std::vector< M > &t) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::DPLEInputIOSchemeVector< M >::DPLEInputIOSchemeVector(const std::vector< M > &t) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::DPLEOutputIOSchemeVector< M >::DPLEOutputIOSchemeVector(const std::vector< M > &t) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::DenseMultiplication::DenseMultiplication(const MX &z, const MX &x, const MX &y) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::DenseTranspose::DenseTranspose(const MX &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::Determinant::Determinant(const MX &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::Diagcat::Diagcat(const std::vector< MX > &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::Diagsplit::Diagsplit(const MX &x, const std::vector< int > &offset1, const std::vector< int > &offset2) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::DleStructIOSchemeVector< T >::DleStructIOSchemeVector(const std::vector< M > &t) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::DleToDple::DleToDple(const DleStructure &st) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::DleToLrDle::DleToLrDle(const DleStructure &st) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::DpleToLrDple::DpleToLrDple(const DpleStructure &st) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::DpleVecStructIOSchemeVector< T >::DpleVecStructIOSchemeVector(const std::vector< M > &t) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::DsdpInterface::DsdpInterface(const std::vector< Sparsity > &st) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::EmptySparsity::EmptySparsity() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::FixedStepIntegrator::FixedStepIntegrator(const Function &f, const Function &g) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::GetNonzeros::GetNonzeros(const Sparsity &sp, const MX &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::GetNonzerosSlice2::GetNonzerosSlice2(const Sparsity &sp, const MX &x, const Slice &inner, const Slice &outer) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::GetNonzerosSlice::GetNonzerosSlice(const Sparsity &sp, const MX &x, const Slice &s) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::GetNonzerosVector::GetNonzerosVector(const Sparsity &sp, const MX &x, const std::vector< int > &nz) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::GradFInputIOSchemeVector< M >::GradFInputIOSchemeVector(const std::vector< M > &t) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::GradFOutputIOSchemeVector< M >::GradFOutputIOSchemeVector(const std::vector< M > &t) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::HNLPInputIOSchemeVector< M >::HNLPInputIOSchemeVector(const std::vector< M > &t) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::HessLagInputIOSchemeVector< M >::HessLagInputIOSchemeVector(const std::vector< M > &t) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::HessLagOutputIOSchemeVector< M >::HessLagOutputIOSchemeVector(const std::vector< M > &t) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::Horzcat::Horzcat(const std::vector< MX > &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::Horzsplit::Horzsplit(const MX &x, const std::vector< int > &offset) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::IdasInterface::IdasInterface(const Function &f, const Function &g) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::ImplicitFixedStepIntegrator::ImplicitFixedStepIntegrator(const Function &f, const Function &g) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::InfSX::InfSX() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::InnerProd::InnerProd(const MX &x, const MX &y) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::IntegratorInputIOSchemeVector< M >::IntegratorInputIOSchemeVector(const std::vector< M > &t) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::IntegratorOutputIOSchemeVector< M >::IntegratorOutputIOSchemeVector(const std::vector< M > &t) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::Inverse::Inverse(const MX &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::IpoptInterface::IpoptInterface(const Function &nlp) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::JacGInputIOSchemeVector< M >::JacGInputIOSchemeVector(const std::vector< M > &t) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::JacGOutputIOSchemeVector< M >::JacGOutputIOSchemeVector(const std::vector< M > &t) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::KinsolInterface::KinsolInterface(const Function &f) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::KnitroInterface::KnitroInterface(const Function &nlp) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::LPStructIOSchemeVector< T >::LPStructIOSchemeVector(const std::vector< M > &t) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::LR_DLEInputIOSchemeVector< M >::LR_DLEInputIOSchemeVector(const std::vector< M > &t) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::LR_DLEOutputIOSchemeVector< M >::LR_DLEOutputIOSchemeVector(const std::vector< M > &t) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::LR_DPLEInputIOSchemeVector< M >::LR_DPLEInputIOSchemeVector(const std::vector< M > &t) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::LR_DPLEOutputIOSchemeVector< M >::LR_DPLEOutputIOSchemeVector(const std::vector< M > &t) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::LapackLuDense::LapackLuDense(const Sparsity &sparsity, int nrhs) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::LapackQrDense::LapackQrDense(const Sparsity &sparsity, int nrhs) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::LinearSolver::LinearSolver() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::LinsolInputIOSchemeVector< M >::LinsolInputIOSchemeVector(const std::vector< M > &t) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::LinsolOutputIOSchemeVector< M >::LinsolOutputIOSchemeVector(const std::vector< M > &t) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::LpSolverInputIOSchemeVector< M >::LpSolverInputIOSchemeVector(const std::vector< M > &t) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::LpSolverOutputIOSchemeVector< M >::LpSolverOutputIOSchemeVector(const std::vector< M > &t) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::LpToQp::LpToQp(const std::vector< Sparsity > &st) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::LrDleStructIOSchemeVector< T >::LrDleStructIOSchemeVector(const std::vector< M > &t) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::LrDleToDle::LrDleToDle(const LrDleStructure &st) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::LrDpleToDple::LrDpleToDple(const LrDpleStructure &st) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::LrDpleToLrDle::LrDpleToLrDle(const LrDleStructure &st, const std::vector< int > &Hs) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::LrDpleVecStructIOSchemeVector< T >::LrDpleVecStructIOSchemeVector(const std::vector< M > &t) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::MX::MX(const Sparsity &sp, double val, bool dummy) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::MX::MX(const Sparsity &sp, int val, bool dummy) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::MX::MX(const std::pair< int, int > &rc) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::MXFunction::MXFunction(const MX &input, const MX &output) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::MXFunction::MXFunction(const MX &input, const std::vector< MX > &output) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::MXFunction::MXFunction(const std::vector< MX > &input, const MX &output) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::MXNode::MXNode() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::Matrix< DataType >::Matrix() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::Matrix< DataType >::Matrix(const Matrix< DataType > &m) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::Matrix< DataType >::Matrix(const Sparsity &sp) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::Matrix< DataType >::Matrix(const Sparsity &sp, const DataType &val, bool dummy) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::Matrix< DataType >::Matrix(const Sparsity &sp, const Matrix< DataType > &d) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::Matrix< DataType >::Matrix(const Sparsity &sp, const std::vector< DataType > &d, bool dummy) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::Matrix< DataType >::Matrix(const std::pair< int, int > &rc) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::Matrix< DataType >::Matrix(const std::vector< DataType > &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::Matrix< DataType >::Matrix(const std::vector< std::vector< double > > &m) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::Matrix< DataType >::Matrix(double val) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::Matrix< DataType >::Matrix(int nrow, int ncol) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::MinusInfSX::MinusInfSX() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::MinusOneSX::MinusOneSX() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::MultipleOutput::MultipleOutput() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::Multiplication::Multiplication(const MX &z, const MX &x, const MX &y) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::NLPInputIOSchemeVector< M >::NLPInputIOSchemeVector(const std::vector< M > &t) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::NLPOutputIOSchemeVector< M >::NLPOutputIOSchemeVector(const std::vector< M > &t) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::NanSX::NanSX() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::Newton::Newton(const Function &f) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::NlpSolverInputIOSchemeVector< M >::NlpSolverInputIOSchemeVector(const std::vector< M > &t) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::NlpSolverOutputIOSchemeVector< M >::NlpSolverOutputIOSchemeVector(const std::vector< M > &t) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::Norm1::Norm1(const MX &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::Norm2::Norm2(const MX &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::Norm::Norm(const MX &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::NormF::NormF(const MX &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::NormInf::NormInf(const MX &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::OldCollocationIntegrator::OldCollocationIntegrator(const Function &f, const Function &g) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::OneSX::OneSX() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::OoqpInterface::OoqpInterface() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::OoqpInterface::OoqpInterface(const std::vector< Sparsity > &st) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::OptionsFunctionalityNode::OptionsFunctionalityNode() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::OutputNode::OutputNode(const MX &parent, int oind) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::QCQPStructIOSchemeVector< T >::QCQPStructIOSchemeVector(const std::vector< M > &t) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::QPStructIOSchemeVector< T >::QPStructIOSchemeVector(const std::vector< M > &t) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::QcqpSolverInputIOSchemeVector< M >::QcqpSolverInputIOSchemeVector(const std::vector< M > &t) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::QcqpSolverOutputIOSchemeVector< M >::QcqpSolverOutputIOSchemeVector(const std::vector< M > &t) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::QcqpToSocp::QcqpToSocp(const std::vector< Sparsity > &st) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::QpSolverInputIOSchemeVector< M >::QpSolverInputIOSchemeVector(const std::vector< M > &t) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::QpSolverOutputIOSchemeVector< M >::QpSolverOutputIOSchemeVector(const std::vector< M > &t) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::QpToImplicit::QpToImplicit(const Function &f) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::QpToNlp::QpToNlp() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::QpToNlp::QpToNlp(const std::vector< Sparsity > &st) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::QpToQcqp::QpToQcqp(const std::vector< Sparsity > &st) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::QpoasesInterface::QpoasesInterface() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::QpoasesInterface::QpoasesInterface(const std::vector< Sparsity > &st) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::RDAEInputIOSchemeVector< M >::RDAEInputIOSchemeVector(const std::vector< M > &t) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::RDAEOutputIOSchemeVector< M >::RDAEOutputIOSchemeVector(const std::vector< M > &t) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::Reshape::Reshape(const MX &x, Sparsity sp) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::RkIntegrator::RkIntegrator(const Function &f, const Function &g) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::RuntimeConst< T >::RuntimeConst() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::RuntimeConst< T >::RuntimeConst(T v) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::SDPInputIOSchemeVector< M >::SDPInputIOSchemeVector(const std::vector< M > &t) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::SDPOutputIOSchemeVector< M >::SDPOutputIOSchemeVector(const std::vector< M > &t) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::SDPStructIOSchemeVector< T >::SDPStructIOSchemeVector(const std::vector< M > &t) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::SDQPInputIOSchemeVector< M >::SDQPInputIOSchemeVector(const std::vector< M > &t) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::SDQPOutputIOSchemeVector< M >::SDQPOutputIOSchemeVector(const std::vector< M > &t) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::SDQPStructIOSchemeVector< T >::SDQPStructIOSchemeVector(const std::vector< M > &t) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::SOCPInputIOSchemeVector< M >::SOCPInputIOSchemeVector(const std::vector< M > &t) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::SOCPOutputIOSchemeVector< M >::SOCPOutputIOSchemeVector(const std::vector< M > &t) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::SOCPStructIOSchemeVector< T >::SOCPStructIOSchemeVector(const std::vector< M > &t) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::SXElement::SXElement() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::SXElement::SXElement(const SXElement &scalar) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::SXElement::SXElement(double val) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::SXFunction::SXFunction(const SX &arg, const SX &res) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::SXFunction::SXFunction(const SX &arg, const std::vector< SX > &res) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::SXFunction::SXFunction(const SX &arg, const std::vector< std::vector< SXElement > > &res) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::SXFunction::SXFunction(const std::vector< SX > &arg, const SX &res) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::SXFunction::SXFunction(const std::vector< std::vector< SXElement > > &arg, const SX &res) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::SXFunction::SXFunction(const std::vector< std::vector< SXElement > > &arg, const std::vector< std::vector< SXElement > > &res) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::SXNode::SXNode() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::ScalarSparseSparsity::ScalarSparseSparsity() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::ScalarSparsity::ScalarSparsity() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::Scpgen::Scpgen(const Function &nlp) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::SdqpToSdp::SdqpToSdp(const std::vector< Sparsity > &st) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::SetNonzeros< Add >::SetNonzeros(const MX &y, const MX &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::SetNonzerosSlice2< Add >::SetNonzerosSlice2(const MX &y, const MX &x, const Slice &inner, const Slice &outer) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::SetNonzerosSlice< Add >::SetNonzerosSlice(const MX &y, const MX &x, const Slice &s) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::SetNonzerosVector< Add >::SetNonzerosVector(const MX &y, const MX &x, const std::vector< int > &nz) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::SetSparse::SetSparse(const MX &x, const Sparsity &sp) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::SharedObject::SharedObject() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::SharedObject::SharedObject(const SharedObject &ref) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::SharedObjectNode::SharedObjectNode() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::SharedObjectNode::SharedObjectNode(const SharedObjectNode &node) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::SimpleHomotopyNlp::SimpleHomotopyNlp(const Function &hnlp) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::SnoptInterface::SnoptInterface(const Function &nlp) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::SocpToSdp::SocpToSdp(const std::vector< Sparsity > &st) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::Solve< Tr >::Solve(const MX &r, const MX &A, const LinearSolver &linear_solver) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::SparseStorage< DataType >::SparseStorage() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::SparseStorage< DataType >::SparseStorage(const SparseStorage< DataType > &m) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::SparseStorage< DataType >::SparseStorage(const Sparsity &sparsity, const DataType &val=DataType(0)) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::Sparsity::Sparsity(const std::pair< int, int > &rc) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::Split::Split(const MX &x, const std::vector< int > &offset) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::SqicInterface::SqicInterface() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::SqicInterface::SqicInterface(const std::vector< Sparsity > &st) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::Sqpmethod::Sqpmethod(const Function &nlp) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::StabilizedQpSolverInputIOSchemeVector< M >::StabilizedQpSolverInputIOSchemeVector(const std::vector< M > &t) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::StabilizedQpToQp::StabilizedQpToQp(const std::vector< Sparsity > &st) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::StabilizedSqicInterface::StabilizedSqicInterface() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::StabilizedSqicInterface::StabilizedSqicInterface(const std::vector< Sparsity > &st) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::StabilizedSqp::StabilizedSqp(const Function &nlp) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::SubAssign::SubAssign(const MX &x, const MX &y, const Slice &i, const Slice &j) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::SubRef::SubRef(const MX &x, const Slice &i, const Slice &j) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::SundialsInterface::SundialsInterface(const Function &f, const Function &g) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::SymbolicMX::SymbolicMX(const std::string &name, const Sparsity &sp) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::SymbolicMX::SymbolicMX(const std::string &name, int nrow=1, int ncol=1) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::SymbolicQr::SymbolicQr(const Sparsity &sparsity, int nrhs) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::SymbolicSX::SymbolicSX(const std::string &name) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::TinyXmlInterface::TinyXmlInterface() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::Transpose::Transpose(const MX &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::UnaryMX::UnaryMX(Operation op, MX x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::Vertcat::Vertcat(const std::vector< MX > &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::Vertsplit::Vertsplit(const MX &x, const std::vector< int > &offset) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::WeakRef::WeakRef(SharedObject shared) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::WeakRef::WeakRef(int dummy=0) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::WorhpInterface::WorhpInterface(const Function &nlp) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::XmlNode::XmlNode() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::ZeroSX::ZeroSX() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::sortCompare< T >::sortCompare(const std::vector< T > &v) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception snoptProblemA::snoptProblemA() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception snoptProblemB::snoptProblemB() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception snoptProblemC::snoptProblemC() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}