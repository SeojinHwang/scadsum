// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#ifndef RCPP_lassosum_RCPPEXPORTS_H_GEN_
#define RCPP_lassosum_RCPPEXPORTS_H_GEN_

#include <RcppArmadillo.h>
#include <Rcpp.h>

namespace lassosum {

    using namespace Rcpp;

    namespace {
        void validateSignature(const char* sig) {
            Rcpp::Function require = Rcpp::Environment::base_env()["require"];
            require("lassosum", Rcpp::Named("quietly") = true);
            typedef int(*Ptr_validate)(const char*);
            static Ptr_validate p_validate = (Ptr_validate)
                R_GetCCallable("lassosum", "_lassosum_RcppExport_validate");
            if (!p_validate(sig)) {
                throw Rcpp::function_not_exported(
                    "C++ function with signature '" + std::string(sig) + "' not found in lassosum");
            }
        }
    }

    inline int countlines(const char* fileName) {
        typedef SEXP(*Ptr_countlines)(SEXP);
        static Ptr_countlines p_countlines = NULL;
        if (p_countlines == NULL) {
            validateSignature("int(*countlines)(const char*)");
            p_countlines = (Ptr_countlines)R_GetCCallable("lassosum", "_lassosum_countlines");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_countlines(Shield<SEXP>(Rcpp::wrap(fileName)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<int >(rcpp_result_gen);
    }

    inline arma::mat multiBed3(const std::string fileName, int N, int P, const arma::mat input, arma::Col<int> col_skip_pos, arma::Col<int> col_skip, arma::Col<int> keepbytes, arma::Col<int> keepoffset, const int trace) {
        typedef SEXP(*Ptr_multiBed3)(SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP);
        static Ptr_multiBed3 p_multiBed3 = NULL;
        if (p_multiBed3 == NULL) {
            validateSignature("arma::mat(*multiBed3)(const std::string,int,int,const arma::mat,arma::Col<int>,arma::Col<int>,arma::Col<int>,arma::Col<int>,const int)");
            p_multiBed3 = (Ptr_multiBed3)R_GetCCallable("lassosum", "_lassosum_multiBed3");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_multiBed3(Shield<SEXP>(Rcpp::wrap(fileName)), Shield<SEXP>(Rcpp::wrap(N)), Shield<SEXP>(Rcpp::wrap(P)), Shield<SEXP>(Rcpp::wrap(input)), Shield<SEXP>(Rcpp::wrap(col_skip_pos)), Shield<SEXP>(Rcpp::wrap(col_skip)), Shield<SEXP>(Rcpp::wrap(keepbytes)), Shield<SEXP>(Rcpp::wrap(keepoffset)), Shield<SEXP>(Rcpp::wrap(trace)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<arma::mat >(rcpp_result_gen);
    }

    inline arma::mat multiBed3sp(const std::string fileName, int N, int P, const arma::vec beta, const arma::Col<int> nonzeros, const arma::Col<int> colpos, const int ncol, arma::Col<int> col_skip_pos, arma::Col<int> col_skip, arma::Col<int> keepbytes, arma::Col<int> keepoffset, const int trace) {
        typedef SEXP(*Ptr_multiBed3sp)(SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP);
        static Ptr_multiBed3sp p_multiBed3sp = NULL;
        if (p_multiBed3sp == NULL) {
            validateSignature("arma::mat(*multiBed3sp)(const std::string,int,int,const arma::vec,const arma::Col<int>,const arma::Col<int>,const int,arma::Col<int>,arma::Col<int>,arma::Col<int>,arma::Col<int>,const int)");
            p_multiBed3sp = (Ptr_multiBed3sp)R_GetCCallable("lassosum", "_lassosum_multiBed3sp");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_multiBed3sp(Shield<SEXP>(Rcpp::wrap(fileName)), Shield<SEXP>(Rcpp::wrap(N)), Shield<SEXP>(Rcpp::wrap(P)), Shield<SEXP>(Rcpp::wrap(beta)), Shield<SEXP>(Rcpp::wrap(nonzeros)), Shield<SEXP>(Rcpp::wrap(colpos)), Shield<SEXP>(Rcpp::wrap(ncol)), Shield<SEXP>(Rcpp::wrap(col_skip_pos)), Shield<SEXP>(Rcpp::wrap(col_skip)), Shield<SEXP>(Rcpp::wrap(keepbytes)), Shield<SEXP>(Rcpp::wrap(keepoffset)), Shield<SEXP>(Rcpp::wrap(trace)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<arma::mat >(rcpp_result_gen);
    }

    inline int elnet(double lambda1, double lambda2, const arma::vec& diag, const arma::mat& X, const arma::vec& r, double thr, arma::vec& x, arma::vec& yhat, int trace, int maxiter) {
        typedef SEXP(*Ptr_elnet)(SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP);
        static Ptr_elnet p_elnet = NULL;
        if (p_elnet == NULL) {
            validateSignature("int(*elnet)(double,double,const arma::vec&,const arma::mat&,const arma::vec&,double,arma::vec&,arma::vec&,int,int)");
            p_elnet = (Ptr_elnet)R_GetCCallable("lassosum", "_lassosum_elnet");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_elnet(Shield<SEXP>(Rcpp::wrap(lambda1)), Shield<SEXP>(Rcpp::wrap(lambda2)), Shield<SEXP>(Rcpp::wrap(diag)), Shield<SEXP>(Rcpp::wrap(X)), Shield<SEXP>(Rcpp::wrap(r)), Shield<SEXP>(Rcpp::wrap(thr)), Shield<SEXP>(Rcpp::wrap(x)), Shield<SEXP>(Rcpp::wrap(yhat)), Shield<SEXP>(Rcpp::wrap(trace)), Shield<SEXP>(Rcpp::wrap(maxiter)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<int >(rcpp_result_gen);
    }

    inline int repelnet(double lambda1, double lambda2, arma::vec& diag, arma::mat& X, arma::vec& r, double thr, arma::vec& x, arma::vec& yhat, int trace, int maxiter, arma::Col<int>& startvec, arma::Col<int>& endvec) {
        typedef SEXP(*Ptr_repelnet)(SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP);
        static Ptr_repelnet p_repelnet = NULL;
        if (p_repelnet == NULL) {
            validateSignature("int(*repelnet)(double,double,arma::vec&,arma::mat&,arma::vec&,double,arma::vec&,arma::vec&,int,int,arma::Col<int>&,arma::Col<int>&)");
            p_repelnet = (Ptr_repelnet)R_GetCCallable("lassosum", "_lassosum_repelnet");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_repelnet(Shield<SEXP>(Rcpp::wrap(lambda1)), Shield<SEXP>(Rcpp::wrap(lambda2)), Shield<SEXP>(Rcpp::wrap(diag)), Shield<SEXP>(Rcpp::wrap(X)), Shield<SEXP>(Rcpp::wrap(r)), Shield<SEXP>(Rcpp::wrap(thr)), Shield<SEXP>(Rcpp::wrap(x)), Shield<SEXP>(Rcpp::wrap(yhat)), Shield<SEXP>(Rcpp::wrap(trace)), Shield<SEXP>(Rcpp::wrap(maxiter)), Shield<SEXP>(Rcpp::wrap(startvec)), Shield<SEXP>(Rcpp::wrap(endvec)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<int >(rcpp_result_gen);
    }

    inline arma::mat genotypeMatrix(const std::string fileName, int N, int P, arma::Col<int> col_skip_pos, arma::Col<int> col_skip, arma::Col<int> keepbytes, arma::Col<int> keepoffset, const int fillmissing) {
        typedef SEXP(*Ptr_genotypeMatrix)(SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP);
        static Ptr_genotypeMatrix p_genotypeMatrix = NULL;
        if (p_genotypeMatrix == NULL) {
            validateSignature("arma::mat(*genotypeMatrix)(const std::string,int,int,arma::Col<int>,arma::Col<int>,arma::Col<int>,arma::Col<int>,const int)");
            p_genotypeMatrix = (Ptr_genotypeMatrix)R_GetCCallable("lassosum", "_lassosum_genotypeMatrix");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_genotypeMatrix(Shield<SEXP>(Rcpp::wrap(fileName)), Shield<SEXP>(Rcpp::wrap(N)), Shield<SEXP>(Rcpp::wrap(P)), Shield<SEXP>(Rcpp::wrap(col_skip_pos)), Shield<SEXP>(Rcpp::wrap(col_skip)), Shield<SEXP>(Rcpp::wrap(keepbytes)), Shield<SEXP>(Rcpp::wrap(keepoffset)), Shield<SEXP>(Rcpp::wrap(fillmissing)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<arma::mat >(rcpp_result_gen);
    }

    inline arma::vec normalize(arma::mat& genotypes) {
        typedef SEXP(*Ptr_normalize)(SEXP);
        static Ptr_normalize p_normalize = NULL;
        if (p_normalize == NULL) {
            validateSignature("arma::vec(*normalize)(arma::mat&)");
            p_normalize = (Ptr_normalize)R_GetCCallable("lassosum", "_lassosum_normalize");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_normalize(Shield<SEXP>(Rcpp::wrap(genotypes)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<arma::vec >(rcpp_result_gen);
    }

    inline List runElnet(arma::vec& lambda, double shrink, const std::string fileName, arma::vec& r, int N, int P, arma::Col<int>& col_skip_pos, arma::Col<int>& col_skip, arma::Col<int>& keepbytes, arma::Col<int>& keepoffset, double thr, arma::vec& x, int trace, int maxiter, arma::Col<int>& startvec, arma::Col<int>& endvec) {
        typedef SEXP(*Ptr_runElnet)(SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP);
        static Ptr_runElnet p_runElnet = NULL;
        if (p_runElnet == NULL) {
            validateSignature("List(*runElnet)(arma::vec&,double,const std::string,arma::vec&,int,int,arma::Col<int>&,arma::Col<int>&,arma::Col<int>&,arma::Col<int>&,double,arma::vec&,int,int,arma::Col<int>&,arma::Col<int>&)");
            p_runElnet = (Ptr_runElnet)R_GetCCallable("lassosum", "_lassosum_runElnet");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_runElnet(Shield<SEXP>(Rcpp::wrap(lambda)), Shield<SEXP>(Rcpp::wrap(shrink)), Shield<SEXP>(Rcpp::wrap(fileName)), Shield<SEXP>(Rcpp::wrap(r)), Shield<SEXP>(Rcpp::wrap(N)), Shield<SEXP>(Rcpp::wrap(P)), Shield<SEXP>(Rcpp::wrap(col_skip_pos)), Shield<SEXP>(Rcpp::wrap(col_skip)), Shield<SEXP>(Rcpp::wrap(keepbytes)), Shield<SEXP>(Rcpp::wrap(keepoffset)), Shield<SEXP>(Rcpp::wrap(thr)), Shield<SEXP>(Rcpp::wrap(x)), Shield<SEXP>(Rcpp::wrap(trace)), Shield<SEXP>(Rcpp::wrap(maxiter)), Shield<SEXP>(Rcpp::wrap(startvec)), Shield<SEXP>(Rcpp::wrap(endvec)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<List >(rcpp_result_gen);
    }

}

#endif // RCPP_lassosum_RCPPEXPORTS_H_GEN_
